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
	{
		set_exponent(0);
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
	 */
	ED_ALWAYS_INLINE
	double fraction() const {
		return m_fraction;
	}

	/**
	 * Returns the exponent of an extended_double.
	 */
	ED_ALWAYS_INLINE
	double exponent() const {
		ieee754_double_t v;
		v.as_double = m_exponent_raw;
		v.as_uint64 ^= EXPONENT_MASK;
		return v.as_double;
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
		/* Multiply fractions and add exponents.
		 * Since a zero value's exponent is -Infinity, +/- Infinity's exponent
		 * +Infinity and NaN's exponent NaN, summing the exponent behaves
		 * correctly for the cases
		 *   0 * 0 = 0 (-Inf + -Inf = -Inf),
		 *   +/-Infinity * +/-Infinity = Infinity ( +Inf + +Inf = +Inf),
		 *   0 * Infinity = NaN (-Inf + +Inf = NaN)
		 *   NaN * ... = NaN (NaN + ...= NaN)
		 */
		double ep = exponent() + v.exponent();
		double fp = fraction() * v.fraction();

		/* Re-normalize if necessary. Since we keep the fractions >= 1.0,
		 * the product of the fractions is surely >= 1.0, so we only need
		 * to check the upper bound.
		 */
		if (fabs(fp) >= FRACTION_RESCALING_THRESHOLD) {
			ep += FRACTION_RESCALING_THRESHOLD_LOG2;
			fp *= FRACTION_RESCALING_THRESHOLD_INV;
		}
		set_exponent(ep);
		set_fraction(fp);

        /* Result should be normalized now */
		check_consistency();
		return *this;
	}

	ED_ALWAYS_INLINE
	extended_double& operator/=(const extended_double& v) {
		/* Divide fractions and subtract exponents.
		 * Since a zero value's exponent is -Infinity, +/- Infinity's exponent
		 * +Infinity and NaN's exponent NaN, summing the exponent behaves
		 * correctly for the cases
		 *   0 / 0 = NaN (-Inf - -Inf = NaN),
		 *   +/-Infinity / +/-Infinity = Infinity ( +Inf - +Inf = NaN),
		 *   0 / Infinity = 0 (-Inf - +Inf = -Inf)
		 *   Infinity / 0 = Infinity (+Inf - -Inf = +Inf)
		 *   NaN / ... = NaN (NaN - ...= NaN)
		 *   ... / NaN = NaN (... - NaN = NaN)
		 */
		double ep = exponent() - v.exponent();
		double fp = fraction() / v.fraction();

		/* Re-normalize if necessary. Since we keep the fraction >= 1.0,
		 * the quotient of the fractions is surely smaller than the dividend's
		 * fraction, so we only need to check the lower bound.
		 */
		if (fabs(fp) < 1.0) {
			ep -= FRACTION_RESCALING_THRESHOLD_LOG2;
			fp *= FRACTION_RESCALING_THRESHOLD;
		}
		set_exponent(ep);
		set_fraction(fp);

		/* Result should be normalized now */
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
	 * Mask exponent is XORed with when stored in m_exponent_raw.
	 */
	static const uint64_t EXPONENT_MASK = 0xfff0000000000000;

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
 	 * Rescaling threshold, i.e. 2^-FRACTION_RESCALING_THRESHOLD_LOG2
 	 */
	static const double FRACTION_RESCALING_THRESHOLD_INV;

	/**
	 * Natural logarithm of 2, used to convert from natural logarithms to base-2 logarithms.
	 */
	static const double LN2;

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
	 *   [ 1.0 , FRACTION_RESCALING_THRESHOLD ).
	 */
	double m_fraction;

	/**
	 * Exponent of an extended_double.
	 *
	 * The exponent is stored XORed with EXPONENT_MASK, which ensures that the
	 * exponent for a zero fraction (-Infinity) has bit pattern 0.
	 *
	 * The exponent is handled similary to IEEE754. For zero values, it is set to
	 * -Infinity, for infinite values it is set to +Infinity, and for NaN
	 * it is set to NaN.
	 */
	double m_exponent_raw;

	/**
	 * Constructs an extended_double from a fraction and an exponent.
	 *
	 * No normalization is performed!
	 */
	ED_ALWAYS_INLINE
	extended_double(double _fraction, double _exponent)
	{
		set_fraction(_fraction);
		set_exponent(_exponent);
		check_consistency();
	}

	ED_ALWAYS_INLINE
	void set_exponent(double exponent) {
		ieee754_double_t v;
		v.as_double = exponent;
		v.as_uint64 ^= EXPONENT_MASK;
		m_exponent_raw = v.as_double;
	}

	ED_ALWAYS_INLINE
	void set_fraction(double fraction) {
		m_fraction = fraction;
	}

	/* XXX: For debugging only */
	double get_exponent();

	void check_consistency() {
#if ED_ENABLE_ASSERTS_NORMALIZATION
		const double INF = std::numeric_limits<double>::infinity();
		assert((fraction() == 0.0) == (exponent() == -INF));
		assert(std::isinf(fraction()) == (exponent() == INF));
		assert(std::isnan(fraction()) == (std::isnan(exponent())));
		if (std::isfinite(exponent())) {
			assert(double(int64_t(exponent())) == exponent());
			assert((int64_t(exponent())
					% extended_double::FRACTION_RESCALING_THRESHOLD_LOG2)
				   == 0);
			assert(std::fabs(fraction()) >= 1.0);
			assert(std::fabs(fraction()) < extended_double::FRACTION_RESCALING_THRESHOLD);
		}
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
     */
	void normalize_slowpath();

	/**
	 * Rescales the fractional part of argument with the smaller absolute value
	 * such that both arguments afterwards have the same exponent.
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
		if (a.exponent() != b.exponent())
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
                                     * extended_double::LN2);
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
