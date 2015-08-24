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

/*** XXX: For debugging only */
double extended_double::get_exponent() {
    return exponent();
}

extended_double extended_double::pow2(int64_t exponent) {
    const int32_t e_mod = (exponent % FRACTION_RESCALING_THRESHOLD_LOG2);
    const int32_t e_nat = (e_mod >= 0) ? e_mod : e_mod + FRACTION_RESCALING_THRESHOLD_LOG2;

    ieee754_double_t f;
    f.as_uint64 = 0;
    f.as_fields.exponent = IEEE754_DOUBLE_EXP_EXCESS + e_nat;
    return extended_double(f.as_double, exponent - e_nat);
}

double extended_double::convert_to_double() const {
    ieee754_double_t r;
    r.as_double = m_fraction;
    const int32_t f_e = int32_t(r.as_fields.exponent) -int32_t(IEEE754_DOUBLE_EXP_EXCESS);
    const double f_ep = std::max(std::min(exponent() + double(f_e),
                                          double(IEEE754_DOUBLE_EXP_EXCESS+1)),
                                 double(-int32_t(IEEE754_DOUBLE_EXP_EXCESS)));
    r.as_fields.exponent = uint32_t(int32_t(f_ep) + int32_t(IEEE754_DOUBLE_EXP_EXCESS));
    return r.as_double;
}


void extended_double::normalize_slowpath() {
    /* Make fields of native double "m_fraction" available */
    ieee754_double_t f;
    f.as_double = m_fraction;
    const uint32_t f_e_raw = f.as_fields.exponent;

    if (f_e_raw == 0) {
        /* Fraction is zero (or denormalized) */
        *this = extended_double();
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

const extended_double::uniformity_factor extended_double::s_uniformity_factors[5] = {
	extended_double::uniformity_factor(0.0       , 1.0       ,  0, -1),
	extended_double::uniformity_factor(exp2(-256), 1.0       ,  0, -1),
	extended_double::uniformity_factor(1.0       , 1.0       , -1,  0),
	extended_double::uniformity_factor(1.0       , exp2(-256), -1,  0),
	extended_double::uniformity_factor(1.0       , 0.0       , -1,  0)
};

void
extended_double::make_exponents_uniform_slowpath(extended_double& a, extended_double& b)
{
#if ED_ENABLE_SSE
    const __m128d EXP_LB = _mm_set1_pd(double(-int32_t(IEEE754_DOUBLE_EXP_EXCESS)));
    const __m128i EXP_EXC = _mm_set_epi16(0, 0, 0, 0,
                                          0, IEEE754_DOUBLE_EXP_EXCESS,
                                          0, IEEE754_DOUBLE_EXP_EXCESS);

    const __m128d e_a = _mm_set_sd(a.exponent());
    const __m128d e_b = _mm_set_sd(b.exponent());
    const __m128d e_ab = _mm_unpacklo_pd(e_a, e_b);
    const __m128d e_ba = _mm_unpacklo_pd(e_b, e_a);
    const __m128d d = _mm_sub_pd(e_ab, e_ba);
    const __m128d lt_ab = _mm_cmplt_pd(e_ab, e_ba);

    const __m128d res_e_ab = _mm_blendv_pd(e_ab, e_ba, lt_ab);
    a.set_exponent(_mm_cvtsd_f64(res_e_ab));
    b.set_exponent(_mm_cvtsd_f64(res_e_ab));

    const __m128i d_sat = _mm_cvtpd_epi32(_mm_max_pd(EXP_LB, _mm_min_pd(d, _mm_setzero_pd())));
    const __m128i dbl_exp = _mm_add_epi16(d_sat, EXP_EXC);
    const __m128i dbl = _mm_slli_epi32(_mm_shuffle_epi32(dbl_exp,
                                                         _MM_SHUFFLE(1, 3, 0, 3)),
                                       31 - IEEE754_DOUBLE_EXP_BITS);

    const __m128d f_ab = _mm_set_pd(b.fraction(), a.fraction());
    const __m128d s = _mm_castsi128_pd(dbl);
    const __m128d res_f_ab = _mm_mul_pd(f_ab, s);
    a.set_fraction(_mm_cvtsd_f64(res_f_ab));
    b.set_fraction(_mm_cvtsd_f64(_mm_unpackhi_pd(res_f_ab, res_f_ab)));
#else
    const double e_delta = a.exponent() - b.exponent();
    const double THRESHOLD = 2*FRACTION_RESCALING_THRESHOLD_LOG2;
    const int32_t e_delta_sat = int32_t(std::min(std::max(-THRESHOLD, e_delta),
                                                 THRESHOLD));
    const int32_t e_delta_idx = (e_delta_sat >> FRACTION_RESCALING_THRESHOLD_LOG2_LOG2);
    ED_ASSERT_NORMALIZATION(e_delta_idx >= -2);
    ED_ASSERT_NORMALIZATION(e_delta_idx <= 2);
    const uniformity_factor& f = s_uniformity_factors[e_delta_idx + 2];

    a.m_fraction *= f.a_fraction_f;
    b.m_fraction *= f.b_fraction_f;

    ieee754_double_t e_a, e_b, e;
    e_a.as_double = a.m_exponent_raw;
    e_b.as_double = b.m_exponent_raw;
    e.as_uint64 = ((e_a.as_uint64 & f.a_exponent_mask)
                   | (e_b.as_uint64 & f.b_exponent_mask));
    a.m_exponent_raw = b.m_exponent_raw = e.as_double;
#endif
}

std::ostream& operator<<(std::ostream& dst, const extended_double& v) {
	if (v.m_exponent_raw == 0)
		dst << '0';
	else {
		int e = 0;
		const double f = std::frexp(v.fraction(), &e);
		dst << f << "*2^" << (v.exponent() + double(e));
	}
    
    return dst;
}
