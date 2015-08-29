#include <iomanip>

#include "extended_double.h"

const uint64_t extended_double::EXPONENT_MASK;

const int32_t extended_double::FRACTION_RESCALING_THRESHOLD_LOG2;

const int32_t extended_double::FRACTION_RESCALING_THRESHOLD_LOG2_LOG2;

const double extended_double::FRACTION_RESCALING_THRESHOLD =
std::ldexp(1.0, extended_double::FRACTION_RESCALING_THRESHOLD_LOG2);

const double extended_double::FRACTION_RESCALING_THRESHOLD_INV =
std::ldexp(1.0, -extended_double::FRACTION_RESCALING_THRESHOLD_LOG2);

const double extended_double::LN2 = std::log(2.0);

const uint32_t extended_double::IEEE754_DOUBLE_EXP_EXCESS;

const uint32_t extended_double::IEEE754_DOUBLE_EXP_INF_RAW;

const uint32_t extended_double::IEEE754_DOUBLE_EXP_MASK;

const uint32_t extended_double::IEEE754_DOUBLE_EXP_BITS;

const uint32_t extended_double::IEEE754_DOUBLE_MAN_BITS;

#ifndef NDEBUG
double extended_double::get_exponent()
{
	return exponent();
}
#endif

extended_double
extended_double::pow2(int64_t exponent)
{
	const int32_t e_mod = (exponent % FRACTION_RESCALING_THRESHOLD_LOG2);
	const int32_t e_nat = (e_mod >= 0) ? e_mod : e_mod + FRACTION_RESCALING_THRESHOLD_LOG2;

	ieee754_double_t f;
	f.as_uint64 = 0;
	f.as_fields.exponent = IEEE754_DOUBLE_EXP_EXCESS + e_nat;
	return extended_double(f.as_double, exponent - e_nat);
}

double
extended_double::convert_to_double() const
{
	ieee754_double_t r;
	r.as_double = m_fraction;

	/* Get fraction's native exponent */
	const int32_t f_e = int32_t(r.as_fields.exponent) -int32_t(IEEE754_DOUBLE_EXP_EXCESS);

	/* Add explicit exponent to fraction's native exponent.
	 * Order of ops carefully chosen so that NaN gets turned into largest value
	 */
	const double f_ep = std::min(double(IEEE754_DOUBLE_EXP_EXCESS+1),
								 std::max(exponent() + double(f_e),
										  double(-int32_t(IEEE754_DOUBLE_EXP_EXCESS))));
	ED_ASSERT_NORMALIZATION(!std::isnan(f_ep));
	ED_ASSERT_NORMALIZATION(!std::isnan(m_fraction) ||
							(f_ep == IEEE754_DOUBLE_EXP_EXCESS+1));

	/* Update fraction's native exponent */
	r.as_fields.exponent = uint32_t(int32_t(f_ep) + int32_t(IEEE754_DOUBLE_EXP_EXCESS));
	return r.as_double;
}


void
extended_double::normalize_slowpath()
{
	/* Make fields of native double "m_fraction" available */
	ieee754_double_t f;
	f.as_double = m_fraction;
	const uint32_t f_e_raw = f.as_fields.exponent;

	if (f_e_raw == 0) {
		/* Fraction is zero (or denormalized).
		 * Make sure the fraction is zero, but keep the sign. Exponent is
		 * set to -infinity (raw value 0)
		 */
		set_fraction(fraction() * 0.0);
		set_exponent_raw(0.0);
	}
	else if (f_e_raw != IEEE754_DOUBLE_EXP_INF_RAW) {
		/* Fraction is non-zero and finite */

		/* Compute native exponent mod FRACTION_RESCALING_THRESHOLD_LOG2.
		 * Result lies within [ 0, FRACTION_RESCALING_THRESHOLD_LOG2 ). This
		 * will become the new native exponent
		 */
		const int32_t f_e_n = (f_e_raw + 1) % FRACTION_RESCALING_THRESHOLD_LOG2;

		/* Compute difference e_delta between original and normalized native
		 * exponent.
		 */
		const int32_t f_e = int32_t(f_e_raw) - int32_t(IEEE754_DOUBLE_EXP_EXCESS);
		const int32_t e_delta = f_e - f_e_n;

		/* Update exponent and fraction. */
		set_exponent(exponent() + double(e_delta));
		f.as_fields.exponent = uint32_t(f_e_n + int32_t(IEEE754_DOUBLE_EXP_EXCESS));
		m_fraction = f.as_double;
	}
	else {
		/* Fraction has value +/- Infinity or NaN. Adjust exponent accordingly */
		set_exponent(std::fabs(m_fraction));
	}
}

#if ED_ENABLE_SSE

void
extended_double::normalize_sum_uniform_exponents_slowpath()
{
	/* Abbreviations for various constants */
	const __m128d TH = _mm_set_sd(FRACTION_RESCALING_THRESHOLD);
	const __m128d TH_INV = _mm_set_sd(FRACTION_RESCALING_THRESHOLD_INV);
	const __m128d TH_LOG2 = _mm_set_sd(FRACTION_RESCALING_THRESHOLD_LOG2);
	const __m128d TH_INV_LOG2 = _mm_set_sd(-FRACTION_RESCALING_THRESHOLD_LOG2);

	/* We assume that the fraction is *not* already normalized, and that it
	 * is the result of adding two normalized fractions. The fraction will
	 * thus either be zero, +/- Infinity, NaN, or lie within
	 *    [ 2^-IEEE754_DOUBLE_MAN_BITS , 1.0)
	 *  U [ FRACTION_RESCALING_THRESHOLD, 2*FRACTION_RESCALING_THRESHOLD )
	 *
	 * It thus suffices to set s to 1/FRACTION_RESCALING_THRESHOLD if the
	 * fraction exceeds the threshold, and to FRACTION_RESCALING_THRESHOLD
	 * otherwise.
	 */
	ED_ASSERT_NORMALIZATION(!((std::fabs(fraction()) >= 1.0)
							  && (std::fabs(fraction()) < FRACTION_RESCALING_THRESHOLD)));
	const __m128d f = fraction_m128d();
	const __m128d f_abs = mm_abs_sd(f);
	const __m128d s_up = _mm_cmplt_sd(f_abs, TH);
	const __m128d f_neq0 = _mm_cmpneq_sd(f, _mm_setzero_pd());
	const __m128d s = _mm_blendv_pd(TH_INV, TH, s_up);
	const __m128d d = _mm_blendv_pd(TH_LOG2, TH_INV_LOG2, s_up);

	/* Multiply fraction with s and add d = log2(s) to exponent. */
	const __m128d fp = _mm_mul_sd(f, s);
	const __m128d ep = _mm_add_sd(exponent_m128d(), d);
	set_fraction<0>(fp);

	/* The above yields the correct result for a fractions other than zero, since
	 * for +/-Infinity and NaN, the sum of exponents computed in operator+=
	 * will have the correct value.
	 *
	 * If the fraction is zero, however, the exponent must be set to -Infinity,
	 * which corresponds to raw value 0. The bitwise AND of the resulting raw
	 * exponent (after adding d, and logical-> raw transformation) ensures this.
	 */
	const __m128d ep_raw = _mm_and_pd(mm_exponent_mask_pd(ep), f_neq0);
	set_exponent_raw<0>(ep_raw);
}

void
extended_double::make_exponents_uniform_slowpath(extended_double& a, extended_double& b)
{
	const double TH_LOG2 = extended_double::FRACTION_RESCALING_THRESHOLD_LOG2;
	const double TH_INV = extended_double::FRACTION_RESCALING_THRESHOLD_INV;

	/* Construct two vectors e_ab = ( a, b ) and e_ba = ( b, a ) */
	const __m128d e_a = _mm_set_sd(a.exponent());
	const __m128d e_b = _mm_set_sd(b.exponent());
	const __m128d e_ab = _mm_unpacklo_pd(e_a, e_b);
	const __m128d e_ba = _mm_unpacklo_pd(e_b, e_a);

	/* Compute delta vector ( a - b, b - a) */
	const __m128d d = _mm_sub_pd(e_ab, e_ba);
	const __m128d d_lt0 = _mm_cmplt_pd(e_ab, e_ba);
	const __m128d d_eqth = _mm_cmpeq_pd(d, _mm_set1_pd(-TH_LOG2));

	/* Compute output exponent vector ( max(a,b), max(a,b) ) and store back */
	const __m128d res_e_ab = _mm_blendv_pd(e_ab, e_ba, d_lt0);
	a.set_exponent<0>(res_e_ab);
	b.set_exponent<0>(res_e_ab);

	/* Compute s =~ min(2^d, 0) by distinguishing the cases d=0, d=-1, d<=-2.
	 * For d <= -2, it suffices to set s = 0, since then a + b = max(a, b).
	 */
	const __m128d s = _mm_or_pd(_mm_or_pd(_mm_setzero_pd(),
										  _mm_andnot_pd(d_lt0, _mm_set1_pd(1.0))),
								_mm_and_pd(d_eqth, _mm_set1_pd(TH_INV)));

	/* Compute output fractions by multiplying with s */
	const __m128d f_ab = _mm_set_pd(b.fraction(), a.fraction());
	const __m128d res_f_ab = _mm_mul_pd(f_ab, s);
	a.set_fraction<0>(res_f_ab);
	b.set_fraction<1>(res_f_ab);
}

void
extended_double::add_nonuniform_exponents_slowpath(const extended_double& v)
{
	/* Various Constants */
	const __m128i TH_LOG2_IEEE754_EXP = _mm_set_epi64x(
													   0, (int64_t(FRACTION_RESCALING_THRESHOLD_LOG2)
														   << IEEE754_DOUBLE_MAN_BITS));
	const __m128d TH_LOG2 = _mm_set_sd(FRACTION_RESCALING_THRESHOLD_LOG2);

	/* Fetch fractions and exponents of the two values and compute difference
	 * between exponents. The difference must be non-zero, otherwise we'd
	 * have taken the fast path!
	 * */
	const __m128d e_a = exponent_m128d();
	const __m128d f_a = fraction_m128d();
	const __m128d e_b = v.exponent_m128d();
	const __m128d f_b = v.fraction_m128d();
	const __m128d e_d = _mm_sub_sd(e_a, e_b);
	const __m128d e_d_abs = mm_abs_sd(e_d);
	ED_ASSERT_NORMALIZATION(!(_mm_cvtsd_f64(e_d_abs) < FRACTION_RESCALING_THRESHOLD_LOG2));

	/* Set e_r to the larger exponent and f_r to the corresponding fraction */
	__m128d e_r = _mm_max_sd(e_a, e_b);
	__m128d f_r = _mm_blendv_pd(f_a, f_b, e_d);

	/* If the smaller exponent is more than one rescaling threshold less
	 * than the larger one, the sum of the two values is (after rounding)
	 * the same as the larger value, and we're done.
	 *
	 * We want the if-branch to be taken if e_d_abs is NaN or not equal to
	 * TH_LOG2. In theory, !_mm_ucomieq(e_d_abs, TH_LOG2) ought to achieve
	 * that, but at least with GCC is doesn't. UCOMISD sets the ZF flags if
	 * either both operands are equal, or if one of them is NaN. GCC enters the
	 * if-branch with JNE, thus branching only if the values are not equal and
	 * neither is NaN. With the operands casted back to double, OTOH, GCC
	 * correctly enters the if-branch with both JNE and JP, thus branching also
	 * if either value is NaN (in which case UCOMISD sets the PF flag).
	 */
	if (!(_mm_cvtsd_f64(e_d_abs) == _mm_cvtsd_f64(TH_LOG2))) {
		/* Compensare for _mm_max_sd's weird NaN handling. If either e_a or
		 * e_b is NaN, we make sure that both fraction and exponent are set
		 * to (quiet) NaN.
		 */
		const __m128d e_r_isnan = _mm_cmpunord_sd(e_a, e_b);
		set_exponent<0>(_mm_or_pd(e_r, e_r_isnan));
		set_fraction<0>(_mm_or_pd(f_r, e_r_isnan));
		return;
	}

	/* Set f_2 to the fraction of the value with the smaller exponent.
	 * Since the absolute difference between the two exponents is *finite* here
	 * (FRACTION_RESCALING_THRESHOLD_LOG2 to be precise), both fractions are
	 * guaranteed to be finite and non-zero (and of course normalized).
	 */
	const __m128d f_2 = _mm_blendv_pd(f_b, f_a, e_d);

	/* Multiply f_2 with FRACTION_RESCALING_THRESHOLD_INV to bring it to the same
	 * scale as f_r (i.e, they'd now both have exponent e_r) and add the
	 * fractions. Multiplying f_r can be done by simply adding to it's native
	 * exponent, since f_r is surely non-zero (see above) and finite.
	 */
	f_r = _mm_add_sd(_mm_castsi128_pd(_mm_sub_epi16(_mm_castpd_si128(f_2),
													TH_LOG2_IEEE754_EXP)),
					 f_r);

	/* |f_r| lay in the range [ 1, FRACTION_RESCALING_THRESHOLD ) and the scaled
	 * |f_2| in the range [ FRACTION_RESCALING_THRESHOLD_INV, 1.0 ). Thus, for
	 * the sum to exceed FRACTION_RESCALING_THRESHOLD, |f_r| would need to
	 * exceed FRACTION_RESCALING_THRESHOLD - 1. But since the threshold is chosen
	 * to be larger than 2^IEEE754_DOUBLE_MAN_BITS, this is only possible if
	 * already |f_r| >= FRACTION_RESCALING_THRESHOLD, which is impossible.
	 * It is therefore not necessary to check for an overflowed fraction here.
	 */
	const __m128d f_r_abs = mm_abs_sd(f_r);
	ED_ASSERT_NORMALIZATION(_mm_cvtsd_f64(f_r_abs) < FRACTION_RESCALING_THRESHOLD);

	/* A check for | f_r + f_2 | < 1.0 *is* necessary, but for similar reasons
	 * as above, | f_r + f_2 | >= 2^-IEEE754_DOUBLE_MAN_BITS. Thus, a single
	 * upwards rescaling round suffices, and there's no need to handle a zero
	 * result.
	 */
	if (_mm_ucomilt_sd(f_r_abs, _mm_set_sd(1.0))) {
		/* Rescale one if necessary to normalize the result. See above. */
		e_r = _mm_sub_sd(e_r, TH_LOG2);
		f_r = _mm_castsi128_pd(_mm_add_epi16(_mm_castpd_si128(f_r),
											 TH_LOG2_IEEE754_EXP));
	}

	set_exponent<0>(e_r);
	set_fraction<0>(f_r);
}

#else // !ED_ENABLE_SSE

void
extended_double::normalize_sum_uniform_exponents_slowpath()
{
	// XXX: We can do better!
	normalize_slowpath();
}

namespace {
	const double RESCALE_TH_INV = std::ldexp(1.0, -extended_double::FRACTION_RESCALING_THRESHOLD_LOG2);

	struct rescale_factors {
		rescale_factors(double _a_fraction_f, double _b_fraction_f,
						int64_t _a_exponent_mask, int64_t _b_exponent_mask)
		:a_fraction_f(_a_fraction_f)
		,b_fraction_f(_b_fraction_f)
		,a_exponent_mask(_a_exponent_mask)
		,b_exponent_mask(_b_exponent_mask)
		{}

		double a_fraction_f, b_fraction_f;
		int64_t a_exponent_mask, b_exponent_mask;
	};

	const rescale_factors rescale_factors_table[5] = {
		rescale_factors(0.0            , 1.0            ,  0, -1),
		rescale_factors(RESCALE_TH_INV , 1.0            ,  0, -1),
		rescale_factors(1.0            , 1.0            , -1,  0),
		rescale_factors(1.0            , RESCALE_TH_INV , -1,  0),
		rescale_factors(1.0            , 0.0            , -1,  0)
	};
}

void
extended_double::make_exponents_uniform_slowpath(extended_double& a, extended_double& b)
{
	const double e_delta = a.exponent() - b.exponent();
	const double THRESHOLD = 2*FRACTION_RESCALING_THRESHOLD_LOG2;
	const int32_t e_delta_sat = int32_t(std::min(std::max(-THRESHOLD, e_delta),
												 THRESHOLD));
	const int32_t e_delta_idx = (e_delta_sat >> extended_double::FRACTION_RESCALING_THRESHOLD_LOG2_LOG2);
	ED_ASSERT_NORMALIZATION(e_delta_idx >= -2);
	ED_ASSERT_NORMALIZATION(e_delta_idx <= 2);
	const rescale_factors& f = rescale_factors_table[e_delta_idx + 2];

	a.set_fraction(a.fraction() * f.a_fraction_f);
	b.set_fraction(b.fraction() * f.b_fraction_f);

	ieee754_double_t e_a, e_b, e_max;
	e_a.as_double = a.m_exponent_raw;
	e_b.as_double = b.m_exponent_raw;
	e_max.as_uint64 = ((e_a.as_uint64 & f.a_exponent_mask)
					   | (e_b.as_uint64 & f.b_exponent_mask));
	a.m_exponent_raw = b.m_exponent_raw = e_max.as_double;
}

void
extended_double::add_nonuniform_exponents_slowpath(const extended_double& v)
{
	/* Fetch fractions and exponents of the two values and compute difference
	 * between exponents. The difference must be non-zero, otherwise we'd
	 * have taken the fast path!
	 * */
	const double e_a = exponent();
	const double f_a = fraction();
	const double e_b = v.exponent();
	const double f_b = v.fraction();
	const double e_d = e_a - e_b;
	ED_ASSERT_NORMALIZATION(e_d != 0.0);

	/* Set e_r to the larger exponent and f_r to the corresponding fraction
	 * If the smaller exponent is more than one rescaling threshold less
	 * than the larger one, the sum of the two values is (after rounding)
	 * the same as the larger value. There's also no need to rescale in this case
	 */
	double e_r = std::max(e_a, e_b);
	double f_r = (e_d < 0.0) ? f_b : f_a;
	if (!(std::fabs(e_d) == FRACTION_RESCALING_THRESHOLD_LOG2)) {
		/* std::max(a,b) does (a < b) ? b : a. Thus, it will propagate a NaN
		 * value of a, but not of b. Therefore, force exponent and fraction to
		 * NaN if e_b is NaN.
		 */
		const bool is_nan = !(e_b == e_b);
		set_exponent(is_nan ? std::numeric_limits<double>::quiet_NaN() : e_r);
		set_fraction(is_nan ? std::numeric_limits<double>::quiet_NaN() : f_r);
		return;
	}

	/* The exponents differ by exactly one rescaling threshold. Thus, the
	 * smaller value's fraction must be divided by FRACTION_RESCALING_THRESHOLD
	 * before adding it to the larger value's fraction.
	 *
	 * The value that is actually added will thus lie within
	 *   [ FRACTION_RESCALING_THRESHOLD_INV, 1 ).
	 * Therefore, for the sum to exceed FRACTION_RESCALING_THRESHOLD, f_r would
	 * need to be at greater than FRACTION_RESCALING_THRESHOLD - 1. Since the
	 * interval [ FRACTION_RESCALING_THRESHOLD - 1, FRACTION_RESCALING_THRESHOLD)
	 * contains no representable numbers if FRACTION_RESCALING_THRESHOLD is
	 * large enough, no check for | f_r + f_2 | >= FRACTION_RESCALING_THRESHOLD
	 * is necessary.
	 *
	 * A check for | f_r + f_2 | < 1.0 *is* necessary, but for similar reasons
	 * as above, | f_r + f_2 | >= 2^-IEEE754_DOUBLE_MAN_BITS. Thus, a single
	 * upwards rescaling round suffices, and there's no need to handle a zero
	 * result.
	 */
	const double f_2 = (e_d < 0.0) ? f_a : f_b;
	f_r += f_2 * FRACTION_RESCALING_THRESHOLD_INV;
	ED_ASSERT_NORMALIZATION(std::fabs(f_r) < FRACTION_RESCALING_THRESHOLD);
	if (std::fabs(f_r) >= 1.0) {
		set_exponent(e_r);
		set_fraction(f_r);
		return;
	}

	/* Rescale one if necessary to normalize the result. See above. */
	e_r -= FRACTION_RESCALING_THRESHOLD_LOG2;
	f_r *= FRACTION_RESCALING_THRESHOLD;
	ED_ASSERT_NORMALIZATION(std::fabs(f_r) >= 1.0);
	set_exponent(e_r);
	set_fraction(f_r);
}

#endif // !ED_ENABLE_SSE


std::ostream&
operator<<(std::ostream& dst, const extended_double& v)
{
	if (v.m_exponent_raw == 0)
		dst << '0';
	else {
		int e = 0;
		const double f = std::frexp(v.fraction(), &e);
		dst << std::setprecision(std::numeric_limits<double>::digits10) << 2.0*f;
		dst << "*2^[" << (v.exponent() + double(e) - 1) << "]";
	}

	return dst;
}
