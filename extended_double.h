#ifndef EXTFLOAT_H
#define EXTFLOAT_H

#include <stdint.h>
#include <cmath>
#include <ostream>
#include <iostream>
#include <limits>
#include <csignal>
#include <cassert>

#define ED_ALWAYS_INLINE inline __attribute__((__always_inline__))
#define ED_UNLIKELY(x) __builtin_expect(!!(x), 0)
#define ED_LIKELY(x) __builtin_expect(!!(x), 1)

#ifndef ED_ENABLE_NAN_WARNING
#   define ED_ENABLE_NAN_WARNING 0
#endif

#ifndef ED_ENABLE_ASSERTS_STATIC
#   define ED_ENABLE_ASSERTS_STATIC 0
#endif
#if ED_ENABLE_ASSERTS_STATIC
#   include <boost/static_assert.hpp>
#   define ED_ASSERT_STATIC(x) BOOST_STATIC_ASSERT(x)
#else
#	define ED_ASSERT_STATIC(x)
#endif

#ifndef ED_ENABLE_ASSERTS_NORMALIZATION
#   define ED_ENABLE_ASSERTS_NORMALIZATION 0
#endif
#if ED_ENABLE_ASSERTS_NORMALIZATION
#   define ED_ASSERT_NORMALIZATION(x) assert(x)
#else
#   define ED_ASSERT_NORMALIZATION(x)
#endif

struct extended_double {
	/**
	 * Default constructor for extended_double, sets the value to 0.
	 */
	ED_ALWAYS_INLINE
	extended_double()
	:m_fraction(0.0)
	,m_exponent_raw(0)
	{}

	/**
	 * Constructor for double -> extended_double conversions.
	 */
	ED_ALWAYS_INLINE
	extended_double(double v)
	:m_fraction(v)
	,m_exponent_raw(EXPONENT_EXCESS)
	{
        const double f = std::fabs(m_fraction);
        if ((f < 1.0) ||
            (f >= FRACTION_RESCALING_THRESHOLD))
            normalize_slowpath();
		check_consistency();
	}
    
    /**
     * Compute 2^i for integral i.
     */
    static extended_double pow2(int64_t exponent);

	/**
	 * Returns the fractional part of an extended_double.
	 *
	 * The absolute value of the fractional part is alwazs >= FRACTION_RESCALING_THRESHOLD
	 */
	ED_ALWAYS_INLINE
	double fraction() const {
		return m_fraction;
	}

	/**
	 * Returns the exponent of an extended_double.
	 *
	 * For non-zero values, the absolute value is always less than EXPONENT_EXCESS.
	 * For zero, the exponent is -EXPONENT_EXCESS.
	 */
	ED_ALWAYS_INLINE
	int64_t exponent() const {
		return m_exponent_raw - EXPONENT_EXCESS;
	}

	ED_ALWAYS_INLINE
	extended_double& operator=(double v) {
		*this = extended_double(v);
		return *this;
	}

	ED_ALWAYS_INLINE
	extended_double& operator+=(extended_double v) {
		make_exponents_uniform(*this, v);
		m_fraction += v.m_fraction;
        const double f = std::fabs(m_fraction);
        if ((f < 1.0) ||
            (f >= FRACTION_RESCALING_THRESHOLD))
            normalize_slowpath();
		check_consistency();
		return *this;
	}

	ED_ALWAYS_INLINE
	extended_double& operator-=(extended_double v) {
        *this += -v;
		return *this;
	}

	ED_ALWAYS_INLINE
	extended_double& operator*=(const extended_double& v) {
        /* Multiply fractions and compute new explicit raw exponent (e_raw_p) */
		m_fraction *= v.m_fraction;
        const int64_t e_raw_p = m_exponent_raw + v.exponent();
        
        /* If we computed 0.0 * x, for |x| < Inf, e_raw_p will now be greater
         * than 0 - more specifically, it will lie within the range
         *   [ 0, EXPONENT_MAX ],
         * when it should actually be zero. On the other hand, for the product
         * of two finite non-zero values e_raw_p will lie within
         *   [ 2*EXPONENT_MIN + EXPONENT_EXCESS , 2*EXPONENT_MAX + EXPONENT_EXCESS ].
         * Since the constants were selected such that
         *   EXPONENT_MAX < 2*EXPONENT_MIN + EXPONENT_EXCESS, the two cases can
         * be distinguished, and we reduce all new raw exponents e_raw_p to
         * zero if they are don't exceed EXPONENT_MAX.
         */
        ED_ASSERT_STATIC(EXPONENT_MAX < 2*EXPONENT_MIN + EXPONENT_EXCESS);
        m_exponent_raw = ((e_raw_p <= EXPONENT_MAX) ? 0 : e_raw_p);
        
        /* XXX: Check for exponent overflow / underflow */
        
        /* During multiplication, the fraction can grew beyond the scaling
         * threshold, but it cannot drop below one if it was previously
         * normalized (since then both factors where >= 1.0). To re-normalize,
         * it thus suffices to test if the fraction's absolute value exceeds
         * the threshold.
         */
		if (!(fabs(m_fraction) < FRACTION_RESCALING_THRESHOLD))
			normalize_after_multiply_slowpath();
        
        /* Result should be normalized now */
		check_consistency();
		return *this;
	}

	ED_ALWAYS_INLINE
	extended_double& operator/=(const extended_double& v) {
		m_fraction /= v.m_fraction;
		m_exponent_raw -= v.exponent();
        const double f = std::fabs(m_fraction);
        if ((f < 1.0) ||
            (f >= FRACTION_RESCALING_THRESHOLD))
            normalize_slowpath();
		check_consistency();
		return *this;
	}

	ED_ALWAYS_INLINE
	extended_double operator-() const {
		return extended_double(-fraction(), exponent());
	}

	friend extended_double fabs(const extended_double& v);
    
	friend double log(const extended_double& v);

	friend bool operator<(extended_double a, extended_double b);
	friend bool operator<=(extended_double a, extended_double b);
	friend bool operator==(extended_double a, extended_double b);
	friend bool operator!=(extended_double a, extended_double b);
	friend bool operator>=(extended_double a, extended_double b);
	friend bool operator>(extended_double a, extended_double b);

	template<typename T>
	friend T extended_double_cast(const extended_double& v);

	friend std::ostream& operator<<(std::ostream& dst, const extended_double& v);

private:
	/**
	 * Excess added to exponents when stored in m_exponent_raw.
	 */
	static const int64_t EXPONENT_EXCESS = 0x2000000000000000; // = 2^61 = -INT64_MIN/2

    /**
     * Maximum logical non-zero value of the exponent field.
     *
     * Maximum physical value is thus EXPONENT_MIN + EXPONENT_EXCESS.
     */
    static const int64_t EXPONENT_MIN = -0x0400000000000000; // = 2^58

    /**
     * Maximum logical non-infinite value of the exponent field.
     *
     * Maximum physical value is thus EXPONENT_MAX + EXPONENT_EXCESS.
     */
    static const int64_t EXPONENT_MAX = 0x0400000000000000; // = 2^58

    /**
     * Logical value for Infinity/NaN values of the exponent field.
     *
     * Physical value is thus EXPONENT_MAX + EXPONENT_EXCESS.
     */
    static const int64_t EXPONENT_INF = 0x2000000000000000; // = 2^61

	/**
	 * Base-2 Logarithm of rescaling threshold.
	 */
	static const int32_t FRACTION_RESCALING_THRESHOLD_LOG2 = 256;

    /**
     * Base-2 double Logarithm of rescaling threshold.
     */
    static const int32_t FRACTION_RESCALING_THRESHOLD_LOG2_LOG2 = 8;
    ED_ASSERT_STATIC((int32_t(1) << FRACTION_RESCALING_THRESHOLD_LOG2_LOG2)
                     == FRACTION_RESCALING_THRESHOLD_LOG2);

	/**
	 * Rescaling threshold, i.e. 2^FRACTION_RESCALING_THRESHOLD_LOG2
	 */
	static const double FRACTION_RESCALING_THRESHOLD;

	/**
	 * Natural logarithm of 2, used to convert from natural logarithms to base-2 logarithms.
	 */
	static const double LOG2;

	/**
	 * Excess used for exponents of IEEE754 double values.
	 */
	static const uint32_t IEEE754_DOUBLE_EXP_EXCESS = 1023;

	/**
	 * Raw exponent of IEEE754 double values +/- Infinity.
 	 */
	static const uint32_t IEEE754_DOUBLE_EXP_INF_RAW = 2047;

	/**
	 * Mask for the exponent field of IEEE754 double values.
	 */
	static const uint32_t IEEE754_DOUBLE_EXP_MASK = 2047;

	static const uint32_t IEEE754_DOUBLE_EXP_BITS = 11;

	static const uint32_t IEEE754_DOUBLE_MAN_BITS = 52;

	/**
	 * Union type used to access the fields of IEEE754 double values.
	 * Note that this is NON-PORTABLE! Its correctness depends on the
	 * bit-field layout chosen by the compiler! It does work, howerver,
	 * for GCC (and compatible compilers) on AMD64.
	 */
	typedef union {
		struct {
			uint64_t mantissa: 52;
			uint64_t exponent: 11;
			uint64_t sign: 1;
		} as_fields ;
		uint64_t as_uint64;
		double as_double;
	} ieee754_double_t;

	/**
	 * Fractional part of an extended_double.
	 *
	 * For normalizes values (i.e., every user-visible instance of extended_double
	 * that isn't zero, +/- Infinity or NaN), m_fraction lie within
	 *   ( FRACTION_RESCALING_THRESHOLD , FRACTION_RESCALING_THRESHOLD_INV ).
	 */
	double m_fraction;

	/**
	 * Exponent of an extended_double.
	 *
	 * The exponent is stored with an excess of EXPONENT_EXCESS, which ensures
	 * that the smallest allowed exponent (-EXPONENT_EXCESS) has bit pattern 0.
	 *
	 * The exponent is handled as in IEEE754. For zero values, it is set to
	 * the smallest allowed value (i.e. -EXPONENT_EXCESS, meaning raw value 0).
	 * For Infinity and +/- NaN, it is set to the largest allowed value,
	 * meaning EXPONENT_MAX (i.e. raw value EXPONENT_EXCESS + EXPONENT_MAX).
	 */
	int64_t m_exponent_raw;

	/**
	 * Constructs an extended_double from a fraction and an exponent.
	 *
	 * No normalization is performed!
	 */
	ED_ALWAYS_INLINE extended_double(double _fraction, int64_t _exponent)
			:m_fraction(_fraction), m_exponent_raw(_exponent + EXPONENT_EXCESS)
	{
		check_consistency();
	}

	void check_consistency() {
#if ED_ENABLE_ASSERTS_NORMALIZATION
		assert((m_exponent_raw == 0) || (m_exponent_raw == EXPONENT_INF + EXPONENT_EXCESS)
			   || ((m_exponent_raw >= EXPONENT_MIN + EXPONENT_EXCESS)
				   && (m_exponent_raw <= EXPONENT_MAX + EXPONENT_EXCESS)));
		assert((m_exponent_raw - EXPONENT_EXCESS) % FRACTION_RESCALING_THRESHOLD_LOG2 == 0);
		assert((m_fraction == 0.0) || (!std::isfinite(m_fraction)) ||
			   ((fabs(m_fraction) >= 1.0) &&
				(fabs(m_fraction) < FRACTION_RESCALING_THRESHOLD)));
		assert((m_exponent_raw == 0) == (m_fraction == 0.0));
		assert((m_exponent_raw == EXPONENT_INF + EXPONENT_EXCESS) == !std::isfinite(m_fraction));
#endif
#if ED_ENABLE_NAN_WARNING
		if (ED_UNLIKELY(!std::isfinite(m_fraction))) {
			std::cerr << "WARNING: NaN produced!" << std::endl;
		}
#endif
	}

    /**
     * Normalizes the fractional part to lie within
     *   [ 1 , FRACTION_RESCALING_THRESHOLD ),
     * and updates the exponent field accordingly.
     *
     * If the fractional value is sub-normal (inkl. zero), the instance will
     * afterwards have numerical value 0, and the exponent will be set to
     * -EXPONENT_EXCESS. If the fractional part is non-finite, the exponent
     * will be set to EXPONENT_MAX.
     */
	void normalize_slowpath();

    /**
     * Faster version of normalize_slowpath() to be used multiplying the
     * fractions and summing the exponents of two previously normalized values.
     */
    void normalize_after_multiply_slowpath();

	void exponent_overflowed();

	/**
	 * Rescales the fractional part of argument with the smaller absolute value such
	 * that both arguments afterwards have the same exponent.
	 *
	 * The argument with the smaller absolute value will afterwards be generally
	 * NOT normalized, i.e. its m_fraction will lie outside of
	 *   [ 1 , FRACTION_RESCALING_THRESHOLD ).
	 * The smaller absolute value may become zero during the rescaling, in which
	 * case the exponent will NOT necessarily be set to the smallest possible
	 * value! Any public function that calls make_exponents_uniform must thus
	 * ensure to eventually normalize its result before returning it.
	 */
	ED_ALWAYS_INLINE
	static void make_exponents_uniform(extended_double& a, extended_double& b) {
		if (a.m_exponent_raw != b.m_exponent_raw)
			make_exponents_uniform_slowpath(a, b);
	}

	static void make_exponents_uniform_slowpath(extended_double& a, extended_double& b);

	struct uniformity_factor {
		uniformity_factor(double _a_fraction_f, double _b_fraction_f,
		                   int64_t _a_exponent_mask, int64_t _b_exponent_mask)
			:a_fraction_f(_a_fraction_f)
			,b_fraction_f(_b_fraction_f)
			,a_exponent_mask(_a_exponent_mask)
			,b_exponent_mask(_b_exponent_mask)
		{}

		double a_fraction_f, b_fraction_f;
		int64_t a_exponent_mask, b_exponent_mask;
	};
	static const uniformity_factor s_uniformity_factors[5];
    
	double convert_to_double() const;
};

ED_ALWAYS_INLINE
extended_double operator+(extended_double a, const extended_double& b) {
	a += b;
	return a;
}

ED_ALWAYS_INLINE
extended_double operator-(extended_double a, const extended_double& b) {
	a -= b;
	return a;
}

ED_ALWAYS_INLINE
extended_double operator*(extended_double a, const extended_double& b) {
	a *= b;
	return a;
}

ED_ALWAYS_INLINE
extended_double operator/(extended_double a, const extended_double& b) {
	a /= b;
	return a;
}

ED_ALWAYS_INLINE
bool operator<(extended_double a, extended_double b) {
	extended_double::make_exponents_uniform(a, b);
	return a.fraction() < b.fraction();
}

ED_ALWAYS_INLINE
bool operator<=(extended_double a, extended_double b) {
	extended_double::make_exponents_uniform(a, b);
	return a.fraction() <= b.fraction();
}

ED_ALWAYS_INLINE
bool operator==(extended_double a, extended_double b) {
	extended_double::make_exponents_uniform(a, b);
	return a.fraction() == b.fraction();
}

ED_ALWAYS_INLINE
bool operator!=(extended_double a, extended_double b) {
	extended_double::make_exponents_uniform(a, b);
	return a.fraction() != b.fraction();
}

ED_ALWAYS_INLINE
bool operator>=(extended_double a, extended_double b) {
	extended_double::make_exponents_uniform(a, b);
	return a.fraction() >= b.fraction();
}

ED_ALWAYS_INLINE
bool operator>(extended_double a, extended_double b) {
	extended_double::make_exponents_uniform(a, b);
	return a.fraction() > b.fraction();
}

ED_ALWAYS_INLINE
extended_double fabs(const extended_double& v) {
	return extended_double(fabs(v.fraction()), v.exponent());
}

ED_ALWAYS_INLINE
double log(const extended_double& v) {
	return std::log(v.fraction()) + (static_cast<double>(v.exponent())
                                     * extended_double::LOG2);
}

ED_ALWAYS_INLINE
double log2(const extended_double& v) {
	return log2(v.fraction()) + static_cast<double>(v.exponent());
}

template<typename T>
T extended_double_cast(const extended_double& v);

template<>
ED_ALWAYS_INLINE
double extended_double_cast<double>(const extended_double& v) {
	return v.convert_to_double();
}

ED_ALWAYS_INLINE
bool isfinite(const extended_double& v) {
	return std::isfinite(v.fraction());
}

#endif
