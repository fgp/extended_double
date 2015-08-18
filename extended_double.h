#ifndef EXTFLOAT_H
#define EXTFLOAT_H

#include <stdint.h>
#include <cmath>
#include <ostream>
#include <iostream>
#include <limits>
#include <csignal>

#define ED_ALWAYS_INLINE inline __attribute__((__always_inline__))
#define ED_UNLIKELY(x) __builtin_expect((x),0)
#define ED_ENABLE_NAN_WARNING 0
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
		normalize_upwards_full();
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

private:
	/**
	 * Rescales the fractional part of argument with the smaller absolute value such
	 * that both arguments afterwards have the same exponent.
	 *
	 * The argument with the smaller absolute value will afterwards be generally
	 * NOT normalized, i.e. its m_fraction may be <= FRACTION_RESCALING_THRESHOLD.
	 * Also, if the value becomes zero during the rescaling, its exponent won't
	 * be set to -EXPONENT_EXCESS! Any public function that calls make_exponents_uniform
	 * must thus ensure to call one of the normalization functions afterwards, to
	 * avoid leaking non-normalized values.
	 */
	ED_ALWAYS_INLINE
	static void make_exponents_uniform(extended_double& a, extended_double& b) {
#if 1
		/* Compute the difference between the two exponents (d_exp), and determine by
		 * which power of two the fractional parts of a resp. b must be multiplied
		 * (a_shift and b_shift). Note that the shift for the value with the larger
		 * exponent is always zero. Also compute saturated versions of these shifts
		 * which fit into a native double's exponent field.
		 */
		const int64_t d_exp = a.m_exponent_raw - b.m_exponent_raw;
		const int64_t a_shift = std::min(d_exp, int64_t(0));
		const int64_t b_shift = std::min(-d_exp, int64_t(0));
		const int64_t a_shift_sat = std::max(a_shift, -int64_t(IEEE754_DOUBLE_EXPONENT_EXCESS));
		const int64_t b_shift_sat = std::max(b_shift, -int64_t(IEEE754_DOUBLE_EXPONENT_EXCESS));

		/* Compute 2^shift for a_shift and b_shift. Instead of using std::pow, the
		 * saturated versions of the shifts are simply written into the exponent of
		 * a native double. If one of the saturated values is IEEE754_DOUBLE_EXPONENT_EXCESS,
		 * the resulting native double is zero, since in IEEE floating-point numbers
		 * zero is encoded as minimal exponent plus zero mantissa.
		 */
		typedef union {
			uint64_t bits;
			double factor;
		} factor_t;
		const factor_t a_factor = { (a_shift_sat + int64_t(IEEE754_DOUBLE_EXPONENT_EXCESS)) << IEEE754_DOUBLE_MANTISSA_BITS };
		const factor_t b_factor = { (b_shift_sat + int64_t(IEEE754_DOUBLE_EXPONENT_EXCESS)) << IEEE754_DOUBLE_MANTISSA_BITS };

		/* Multiple a and b by 2^shift_a resp. 2^shift_b, and adjust the exponents accordingly.
		 * Since either a_shift = a.exp - b.exp, or b_shift = b.exp - a.exp, the exponents
		 * of a and b will aftewards agree.
		 */
		a.m_fraction *= a_factor.factor;
		b.m_fraction *= b_factor.factor;
		a.m_exponent_raw -= a_shift;
		b.m_exponent_raw -= b_shift;
#else
		const int64_t a_e = a.exponent();
		const int64_t b_e = b.exponent();
		const int64_t e = std::max(a_e, b_e);
		a = extended_double(a.m_fraction * std::pow(2.0, a_e - e), e);
		b = extended_double(b.m_fraction * std::pow(2.0, b_e - e), e);
#endif
	}

	/**
	 * Performs a single rescaling step if necessary.
	 *
	 * If the fraction does not exceed FRACTION_RESCALING_THRESHOLD...
	 */
	void normalize_upwards_once() {
		/* If the value is zero, the exponent is set to logical value -EXPONENT_EXCESS,
		 * which corresponds to raw value 0. See make_exponent_uniform as to why.
		 * Otherwise, if the fractional part does not exceed FRACTION_RESCALING_THRESHOLD,
		 * it is rescaled and the exponent is increased accordingly. Since the
		 * raw and logical value of the exponent differes only by an additive constant,
		 * its not necessary to convert from raw to logical representation and back
		 * to do the addition.
		 */
		if (ED_UNLIKELY(m_fraction == 0.0))
			m_exponent_raw = 0;
		else if (ED_UNLIKELY(fabs(m_fraction) <= FRACTION_RESCALING_THRESHOLD)) {
			m_fraction *= FRACTION_RESCALING_THRESHOLD_INV;
			m_exponent_raw += FRACTION_RESCALING_THRESHOLD_LOG2;
		}

#if ED_ENABLE_FPE
		if (ED_UNLIKELY(m_exponent_raw < 0))
			raise(SIGFPE);
#endif
	}

	/**
	 * Performs as many rescaling steps as necessary to normalize the value
	 *
	 * If the fraction does not exceed FRACTION_RESCALING_THRESHOLD...
	 */
	void normalize_upwards_full() {
		/* If the value is zero, the exponent is set to logical value -EXPONENT_EXCESS,
		 * which corresponds to raw value 0. See make_exponent_uniform as to why.
		 * Otherwise, as long as the fractional part does not exceed FRACTION_RESCALING_THRESHOLD,
		 * it is rescaled and the exponent is increased accordingly. Since the
		 * raw and logical value of the exponent differes only by an additive constant,
		 * its not necessary to convert from raw to logical representation and back
		 * to do the addition.
		 */
		if (ED_UNLIKELY(m_fraction == 0.0))
			m_exponent_raw = 0;
		else {
			while (ED_UNLIKELY(fabs(m_fraction) <= FRACTION_RESCALING_THRESHOLD)) {
				m_fraction *= FRACTION_RESCALING_THRESHOLD_INV;
				m_exponent_raw += FRACTION_RESCALING_THRESHOLD_LOG2;
			}
		}

#if ED_ENABLE_FPE
		if (ED_UNLIKELY(m_exponent_raw < 0))
			raise(SIGFPE);
#endif
	}

public:
	ED_ALWAYS_INLINE
	extended_double& operator=(double v) {
		*this = extended_double(v);
		check_consistency();
		return *this;
	}

	ED_ALWAYS_INLINE
	extended_double& operator+=(extended_double v) {
		make_exponents_uniform(*this, v);
		m_fraction += v.m_fraction;
		normalize_upwards_full();
		check_consistency();
		return *this;
	}

	ED_ALWAYS_INLINE
	extended_double& operator-=(const extended_double& v) {
		return operator+=(-v);
	}
	
	ED_ALWAYS_INLINE
	extended_double& operator*=(const extended_double& v) {
		m_fraction *= v.m_fraction;
		m_exponent_raw += v.exponent();
		normalize_upwards_once();
		check_consistency();
		return *this;
	}
	
	ED_ALWAYS_INLINE
	extended_double& operator/=(const extended_double& v) {
		m_fraction /= v.m_fraction;
		m_exponent_raw -= v.exponent();
		normalize_upwards_full();
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

	static const int64_t IEEE754_DOUBLE_EXPONENT_EXCESS = 1023;

	static const unsigned IEEE754_DOUBLE_MANTISSA_BITS = 52;

	/**
	 * Fractional part of an extended_double.
	 *
	 * For normalizes values (i.e., every user-visible instance of extended_double),
	 * the absolute value of m_fraction is greater than FRACTION_RESCALING_THRESHOLD.
	 * FRACTION_RESCALING_THRESHOLD is chosen such that the product of two fractional
	 * parts (as the product of two doubles) doesn't underflow to zero.
	 *
	 * Currently, no upper bound is imposed on m_fraction, meaning that products MAY
	 * overflow to +/- infinity!
	 */
	double m_fraction;
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
#if ED_ENABLE_NAN_WARNING
		if (ED_UNLIKELY(!std::isfinite(m_fraction))) {
			std::cerr << "WARNING: NaN produced!" << std::endl;
		}
#endif
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
#if 1
	const uint64_t e_sat = std::max(std::min(v.exponent(),
											 extended_double::IEEE754_DOUBLE_EXPONENT_EXCESS + 1),
									-extended_double::IEEE754_DOUBLE_EXPONENT_EXCESS);

	typedef union {
		uint64_t bits;
		double factor;
	} factor_t;
	const factor_t factor = { ((e_sat + extended_double::IEEE754_DOUBLE_EXPONENT_EXCESS)
							   << extended_double::IEEE754_DOUBLE_MANTISSA_BITS) };

	return v.fraction() * factor.factor;
#else
	return v.fraction() * std::pow(2.0, v.exponent());
#endif
}

ED_ALWAYS_INLINE
bool isfinite(const extended_double& v) {
	return std::isfinite(v.fraction());
}

#endif
