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
#   ifdef NDEBUG
#       warning NDEBUG overrides ED_ENABLE_ASSERTS_NORMALIZATION
#   endif
#   define ED_ASSERT_NORMALIZATION(x) assert(x)
#else
#   define ED_ASSERT_NORMALIZATION(x)
#endif

#ifndef ED_ENABLE_SSE
#	define ED_ENABLE_SSE 1
#endif
#if ED_ENABLE_SSE
#	include <xmmintrin.h>
#	include <emmintrin.h>
#	include <pmmintrin.h>
#	include <tmmintrin.h>
#	include <smmintrin.h>
#endif

struct extended_double {
	/**
	 * Base-2 Logarithm of rescaling threshold.
	 */
	static const int32_t FRACTION_RESCALING_THRESHOLD_LOG2 = 512;

	/**
	 * Base-2 double Logarithm of rescaling threshold.
	 */
	static const int32_t FRACTION_RESCALING_THRESHOLD_LOG2_LOG2 = 9;
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
		/* The weird double-negation ensures that the test succeeds for NaN */
		if (!((f >= 1.0) &&
			  (f < FRACTION_RESCALING_THRESHOLD)))
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
#if ED_ENABLE_SSE
		const __m128d m =  _mm_load_sd(reinterpret_cast<const double*>(&EXPONENT_MASK));
		return _mm_cvtsd_f64(_mm_xor_pd(_mm_set_sd(m_exponent_raw), m));
#else
		ieee754_double_t v;
		v.as_double = m_exponent_raw;
		v.as_uint64 ^= EXPONENT_MASK;
		return v.as_double;
#endif
	}

	/**
	 * Assigns a double value to an extended_double
	 */
	ED_ALWAYS_INLINE
	extended_double& operator=(double v) {
		*this = extended_double(v);
		return *this;
	}

	ED_ALWAYS_INLINE
	extended_double& operator+=(const extended_double& v) {
		/* If exponents are equal, fractions are added and the sum renormalized,
		 * Otherwise, add_nonuniform_exponents_slowpath() rescales the fractions
		 * as necessary, adds, and normalizes.
		 *
		 * Exponents are tested for bitwise equality. This means that for NaN,
		 * both code paths could be taken, since multiple bit-patters mean NaN!
		 * Note that the exponent musn't be -0.0 (+0.0 is of course allowed!),
		 * otherwise we might call normalize_sum_uniform_exponents_slowpath()
		 * when the exponents are actually the same, which is doesn't handle!
		 */
		if (are_exponents_uniform(*this, v)) {
#if ED_ENABLE_SSE
			/* Sum fractions */
			const __m128d f_a = _mm_set_sd(m_fraction);
			const __m128d f_b = _mm_set_sd(v.m_fraction);
			const __m128d f_r = _mm_add_sd(f_a, f_b);
			set_fraction<0>(f_r);

			/* Handle over- or underflowed fraction.
			 * Adding 0x401 to the fraction's native exponent maps native
			 * exponents within [0, 511] to bit patterns 00xxxxxxxxx. The
			 * mask of 0x600 = 11000000000 ensures that only the leading bits
			 * are tested.
			 */
			const __m128i IEEE754_CMP_ADJ =
			_mm_set_epi64x(0, (int64_t(0x401) << IEEE754_DOUBLE_MAN_BITS));
			const __m128i IEEE754_CMP_MASK =
			_mm_set_epi64x(0, (int64_t(0x600) << IEEE754_DOUBLE_MAN_BITS));
			const __m128i f_r_cmp = _mm_add_epi16(_mm_castpd_si128(f_r),
												  IEEE754_CMP_ADJ);
			if (!_mm_test_all_zeros(f_r_cmp, IEEE754_CMP_MASK))
				normalize_sum_uniform_exponents_slowpath();
#else
			/* Sum fractions */
			set_fraction(fraction() + v.fraction());
			const double f_abs = std::fabs(fraction());

			/* Handle over- or underflowed fraction */
			if ((f_abs >= FRACTION_RESCALING_THRESHOLD)
				|| (f_abs < 1.0))
				normalize_sum_uniform_exponents_slowpath();
#endif
		}
		else
			add_nonuniform_exponents_slowpath(v);

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
#if ED_ENABLE_SSE
		m_exponent_raw = _mm_cvtsd_f64(mm_exponent_mask_pd(_mm_set_sd(exponent)));
#else
		ieee754_double_t v;
		v.as_double = exponent;
		v.as_uint64 ^= EXPONENT_MASK;
		m_exponent_raw = v.as_double;
#endif
	}

	ED_ALWAYS_INLINE
	void set_exponent_raw(double exponent) {
		m_exponent_raw = exponent;
	}

	ED_ALWAYS_INLINE
	void set_fraction(double fraction) {
		m_fraction = fraction;
	}

#if ED_ENABLE_SSE
	/**
	 * Returns exponent in the low double of a __m128d SSE vector
	 */
	ED_ALWAYS_INLINE __m128d exponent_m128d() const {
		return mm_exponent_mask_pd(_mm_load_sd(&m_exponent_raw));
	}

	/**
	 * Returns fraction in the low double of a __m128d SSE vector
	 */
	ED_ALWAYS_INLINE __m128d fraction_m128d() const {
		return _mm_load_sd(&m_fraction);
	}

	/**
	 * Sets exponent to low (I==0) or high (I==1) double in __m128d SSE vector.
	 */
	template<int I> ED_ALWAYS_INLINE
	void set_exponent(__m128d exponent) {
		ED_ASSERT_STATIC((I == 0) || (I == 1));
		const __m128d exponent_raw = mm_exponent_mask_pd(exponent);
		if (I == 0)
			_mm_store_sd(&m_exponent_raw, exponent_raw);
		else
			_mm_storeh_pd(&m_exponent_raw, exponent_raw);
	}

	/**
	 * Sets raw exponent to low (I==0) or high (I==1) double in __m128d SSE vector.
	 */
	template<int I> ED_ALWAYS_INLINE
	void set_exponent_raw(__m128d exponent_raw) {
		ED_ASSERT_STATIC((I == 0) || (I == 1));
		if (I == 0)
			_mm_store_sd(&m_exponent_raw, exponent_raw);
		else
			_mm_storeh_pd(&m_exponent_raw, exponent_raw);
	}

	/**
	 * Sets fraction to low (I==0) or high (I==1) double in __m128d SSE vector.
	 */
	template<int I> ED_ALWAYS_INLINE
	void set_fraction(__m128d fraction) {
		ED_ASSERT_STATIC((I == 0) || (I == 1));
		if (I == 0)
			_mm_store_sd(&m_fraction, fraction);
		else
			_mm_storeh_pd(&m_fraction, fraction);
	}

	/**
	 * Transforms raw to logical exponent and vice versa.
	 *
	 * This is used to ensure that an all-zero extended_double represents
	 * value 0, which requires the logical exponent to be -infinity.
	 */
	ED_ALWAYS_INLINE static __m128d mm_exponent_mask_pd(const __m128d e) {
		const __m128d M = _mm_set1_pd(-std::numeric_limits<double>::infinity());
		return _mm_xor_pd(e, M);
	}

	/**
	 * Returns the element-wise absolute value of an __m128d SSE vector
	 */
	ED_ALWAYS_INLINE static __m128d mm_abs_sd(const __m128d v) {
		const __m128d M = _mm_set_sd(-0.0);
		return _mm_andnot_pd(M, v);
	}
#endif

#ifndef NDEBUG
	/**
	 * Allows the logical exponent to be easy retrieved in GDB.
	 */
	double get_exponent();
#endif

	/**
	 * Returns true if the two values have the same exponent.
	 *
	 * Note that for NaN values, the result may either be true or false,
	 * depending on whether SSE is enabled, and on whether both NaN are
	 * represented by exactly the same bit pattern.
	 */
	ED_ALWAYS_INLINE
	static bool
	are_exponents_uniform(const extended_double& a, const extended_double& b) {
#if ED_ENABLE_SSE
		/* The bitwise comparison has indeterminate result for NaN exponents! */
		const __m128d a_e = _mm_set_sd(a.m_exponent_raw);
		const __m128d b_e = _mm_set_sd(b.m_exponent_raw);
		const __m128i eq = _mm_castpd_si128(_mm_xor_pd(a_e, b_e));
		return _mm_test_all_zeros(eq, eq);
#else
		/* It's tempting to compare the raw exponents here, but that is dangerous.
		 * Due to the XOR mask, the raw exponents might appear to be NaN, which
		 * have strange comparison semantics. Note that it *might* still be safe
		 * to use the raw exponents here, but as long as there's no proof that is
		 * is, let's rather be safe than sorry.
		 */
		return (a.exponent() == b.exponent());
#endif
	}

	/**
	 * Make sure the two values have the same exponent while retaining their
	 * numerical value (as far as possible).
	 *
	 * The smaller (absolute) value may be rounded down to zero, if the
	 * exponent difference is large.
	 */
	ED_ALWAYS_INLINE
	static void make_exponents_uniform(extended_double& a, extended_double& b) {
		if (!are_exponents_uniform(a, b))
			make_exponents_uniform_slowpath(a, b);
	}

	void check_consistency() const {
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

	void normalize_slowpath();

	void normalize_sum_uniform_exponents_slowpath();

	void add_nonuniform_exponents_slowpath(const extended_double& v);

	static void make_exponents_uniform_slowpath(extended_double& a, extended_double& b);

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

ED_ALWAYS_INLINE
bool isnan(const extended_double& v) {
	return std::isnan(v.fraction());
}

#endif
