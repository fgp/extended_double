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
#define ED_ENABLE_NAN_WARNING 0
#define ED_ENABLE_ASSERTS_NORMALIZATION 1
#define ED_ENABLE_FPE 1

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
	static const int64_t EXPONENT_EXCESS = 0x4000000000000000; // = 2^62 = -INT64_MIN/2

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
	 * that isn't zero, +/- Infinity or NaN), m_fraction lie within [ 1, 2 ).
	 */
	double m_fraction;

	/**
	 * Exponent of an extended_double.
	 *
	 * The exponent is stored with an excess of EXPONENT_EXCESS, which ensures
	 * that the smallest allowed exponent (-EXPONENT_EXCESS) has bit pattern 0.
	 *
	 * For non-normal values (i.e. Zero, +/- Infinity, Nan) the exponent's raw
	 * value is zero. Note that therfore, contrary to how things are done in
	 * IEEE754, the logical value of the exponent is -EXPONENT_EXCESS for
	 * +/- Infinity and NaN!
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
			   ((fabs(m_fraction) >= 1.0) && (fabs(m_fraction) < 2.0)));
#endif
#if ED_ENABLE_NAN_WARNING
		if (ED_UNLIKELY(!std::isfinite(m_fraction))) {
			std::cerr << "WARNING: NaN produced!" << std::endl;
		}
#endif
	}

	/**
	 * Normalizes the fractional part to lie within [ 1, 2 ), and updates the
	 * exponent field accordingly to keep the numerical value the same.
	 *
	 * If the fractional value is sub-normal, the instance will afterwards
	 * have numerical value 0.
	 */
	ED_ALWAYS_INLINE
	void normalize() {
		/* Make fields of native double "m_fraction" available */
		ieee754_double_t v;
		v.as_double = m_fraction;

		/* Extract native exponent and determine if v is non-finite or sub-normal.
		 * For sub-normal values, e_type is 0, while for inifinities its 1.
		 * Note that zero is counted as being sub-normal here - as other sub-normal
		 * numbers, its native exponent is the smallest possible value.
		 */
		const uint32_t e = v.as_fields.exponent;
		const uint32_t e_type = (e + 1) & IEEE754_DOUBLE_EXP_MASK;

		/* Update exponent field. For non-finite or subnormal values, its
		 * set to the smallest possible value (even for +/- Inf!)
		 */
		m_exponent_raw += e - int64_t(IEEE754_DOUBLE_EXP_EXCESS);
		m_exponent_raw = (e_type <= 1) ? 0 : m_exponent_raw;

		/* Update m_fraction's native exponent and mantissa. For non-finite
		 * and sub-normal values, the original exponent is kept - this preserves
		 * Zero, +/- Infinity and NaN. For all sub-normal values, the mantissa
		 * is set to 0, i.e. they are converted to zero here! Note that doing
		 * this for all non-finite values would convert NaN into +/- Inf...
		 * For finite (and non-zero) values, the exponent is set to 0 (i.e. raw
		 * value IEEE754_DOUBLE_EXP_EXCESS).
		 */
		const uint32_t ep = (e_type <= 1) ? e : IEEE754_DOUBLE_EXP_EXCESS;
		const uint64_t mp = e_type ? v.as_fields.mantissa : 0 ;
		v.as_fields.exponent = ep;
		v.as_fields.mantissa = mp;
		m_fraction = v.as_double;
	}

	/**
	 * Rescales the fractional part of argument with the smaller absolute value such
	 * that both arguments afterwards have the same exponent.
	 *
	 * The argument with the smaller absolute value will afterwards be generally
	 * NOT normalized, i.e. its m_fraction will lie outside of [ 1, 2 ). The
	 * smaller absolute value may become zero during the rescaling, in which
	 * case the exponent will NOT necessarily be set to the smallest possible
	 * value! Any public function that calls make_exponents_uniform must thus
	 * ensure to eventually normalize its result before returning it.
	 */
	ED_ALWAYS_INLINE
	static void make_exponents_uniform(extended_double& a, extended_double& b) {
		/* Compute the difference between the two exponents (d_exp), and determine by
		 * which power of two the fractional parts of a resp. b must be multiplied
		 * (a_shift and b_shift). Note that the shift for the value with the larger
		 * exponent is always zero. Also compute saturated versions of these shifts
		 * which fit into a native double's exponent field.
		 */
		const int64_t d_exp = a.m_exponent_raw - b.m_exponent_raw;
		const int64_t a_shift = std::min(d_exp, int64_t(0));
		const int64_t b_shift = std::min(-d_exp, int64_t(0));
		const int64_t a_shift_sat = std::max(a_shift, -int64_t(IEEE754_DOUBLE_EXP_EXCESS));
		const int64_t b_shift_sat = std::max(b_shift, -int64_t(IEEE754_DOUBLE_EXP_EXCESS));

		/* Compute 2^shift for a_shift and b_shift. Instead of using std::pow, the
		 * saturated versions of the shifts are simply written into the exponent of
		 * a native double. If one of the saturated values is IEEE754_DOUBLE_EXPONENT_EXCESS,
		 * the resulting native double is zero, since in IEEE floating-point numbers
		 * zero is encoded as minimal exponent plus zero mantissa.
		 */
		ieee754_double_t a_factor, b_factor;
		a_factor.as_uint64 = b_factor.as_uint64 = 0;
		a_factor.as_fields.exponent = a_shift_sat + int64_t(IEEE754_DOUBLE_EXP_EXCESS);
		b_factor.as_fields.exponent = b_shift_sat + int64_t(IEEE754_DOUBLE_EXP_EXCESS);

		/* Multiple a and b by 2^shift_a resp. 2^shift_b, and adjust the exponents accordingly.
		 * Since either a_shift = a.exp - b.exp, or b_shift = b.exp - a.exp, the exponents
		 * of a and b will afterwards agree.
		 */
		a.m_fraction *= a_factor.as_double;
		b.m_fraction *= b_factor.as_double;
		a.m_exponent_raw -= a_shift;
		b.m_exponent_raw -= b_shift;
	}

	ED_ALWAYS_INLINE
	double convert_to_double() const {
		const int32_t e_sat = std::max(std::min(exponent(),
												int64_t(IEEE754_DOUBLE_EXP_EXCESS + 1)),
									   -int64_t(IEEE754_DOUBLE_EXP_EXCESS - 1));

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
