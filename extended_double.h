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
#define ED_UNLIKELY(x) __builtin_expect((x),0)

#ifndef ED_ENABLE_NAN_WARNING
#define ED_ENABLE_NAN_WARNING 0
#endif

#ifndef ED_ENABLE_ASSERTS_NORMALIZATION
#define ED_ENABLE_ASSERTS_NORMALIZATION 0
#endif

#ifndef ED_ENABLE_ASSERTS_OVERFLOW
#define ED_ENABLE_ASSERTS_OVERFLOW 0
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
		normalize();
		check_consistency();
	}

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
		normalize();
		check_consistency();
		return *this;
	}

	ED_ALWAYS_INLINE
	extended_double& operator-=(extended_double v) {
		make_exponents_uniform(*this, v);
		m_fraction -= v.m_fraction;
		normalize();
		check_consistency();
		return *this;
	}

	ED_ALWAYS_INLINE
	extended_double& operator*=(const extended_double& v) {
		m_fraction *= v.m_fraction;
		m_exponent_raw += v.exponent();
		normalize();
		check_consistency();
		return *this;
	}

	ED_ALWAYS_INLINE
	extended_double& operator/=(const extended_double& v) {
		m_fraction /= v.m_fraction;
		m_exponent_raw -= v.exponent();
		normalize();
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

	static const int64_t EXPONENT_MAX = 0x2000000000000000; // = 2^61

	/**
	 * Base-2 Logarithm of rescaling threshold.
	 */
	static const int64_t FRACTION_RESCALING_THRESHOLD_LOG2 = -256;

	/**
	 * Rescaling threshold, i.e. 2^FRACTION_RESCALING_THRESHOLD_LOG2
	 */
	static const double FRACTION_RESCALING_THRESHOLD;

	/**
	 * Multiplicative inverse of FRACTION_RESCALING_THRESHOLD, i.e. 2^-FRACTION_RESCALING_THRESHOLD_LOG2
	 */
	static const double FRACTION_RESCALING_THRESHOLD_INV;

	/**
	 * Natural logarithm of 2, used to convert from natural logarithms to base-2 logarithms.
	 */
	static const double LOG2;

	/**
	 * Excess used for exponents of IEEE754 double values.
	 */
	static const uint32_t IEEE754_DOUBLE_EXP_EXCESS = 1023;

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
		assert((m_fraction == 0.0) || (!std::isfinite(m_fraction)) ||
			   ((fabs(m_fraction) > FRACTION_RESCALING_THRESHOLD) &&
				(fabs(m_fraction) < FRACTION_RESCALING_THRESHOLD_INV)));
#endif
#if ED_ENABLE_ASSERTS_OVERFLOW
		assert((m_exponent_raw - EXPONENT_EXCESS) % -FRACTION_RESCALING_THRESHOLD_LOG2 == 0);
		assert((m_exponent_raw >= 0));
		assert((m_exponent_raw <= EXPONENT_EXCESS + EXPONENT_MAX));
#endif
#if ED_ENABLE_NAN_WARNING
		if (ED_UNLIKELY(!std::isfinite(m_fraction))) {
			std::cerr << "WARNING: NaN produced!" << std::endl;
		}
#endif
	}

	/**
	 * Normalizes the fractional part to lie within
	 *   ( FRACTION_RESCALING_THRESHOLD , FRACTION_RESCALING_THRESHOLD_INV ),
	 * and updates the exponent field accordingly.
	 *
	 * If the fractional value is sub-normal, the instance will afterwards
	 * have numerical value 0.
	 */
	ED_ALWAYS_INLINE
	void normalize() {
#if 0
		const double f = fabs(m_fraction);
		if ((f <= FRACTION_RESCALING_THRESHOLD) ||
            (f >= FRACTION_RESCALING_THRESHOLD_INV))
			normalize_slowpath();
#else
        ieee754_double_t v;
        v.as_double = m_fraction;
        const uint32_t e = (v.as_fields.exponent
                            - IEEE754_DOUBLE_EXP_EXCESS
                            - FRACTION_RESCALING_THRESHOLD_LOG2 - 1);
        if (e >= (-2*FRACTION_RESCALING_THRESHOLD_LOG2 - 1))
            normalize_slowpath();
#endif
	}

	void normalize_slowpath();

	/**
	 * Rescales the fractional part of argument with the smaller absolute value such
	 * that both arguments afterwards have the same exponent.
	 *
	 * The argument with the smaller absolute value will afterwards be generally
	 * NOT normalized, i.e. its m_fraction will lie outside of
	 *   ( FRACTION_RESCALING_THRESHOLD , FRACTION_RESCALING_THRESHOLD_INV ).
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

	ED_ALWAYS_INLINE
	double convert_to_double() const {
		const int32_t e_sat = std::max(std::min(exponent(),
												int64_t(IEEE754_DOUBLE_EXP_EXCESS + 1)),
									   -int64_t(IEEE754_DOUBLE_EXP_EXCESS));

		ieee754_double_t factor;
		factor.as_uint64 = 0;
		factor.as_fields.exponent = (e_sat + int32_t(IEEE754_DOUBLE_EXP_EXCESS));

		return fraction() * factor.as_double;
	}
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
	return log(v.fraction()) + static_cast<double>(v.exponent())*extended_double::LOG2;
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
