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
double extended_double::get_exponent() {
    return exponent();
}
#endif

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

#if !ED_ENABLE_SSE
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
#endif

#if !ED_ENABLE_SSE

ED_ALWAYS_INLINE
void
extended_double::rescale_fractions(const extended_double& a, const extended_double& b,
                                   double& a_f, double& b_f, double& e_raw)
{
    const double e_delta = a.exponent() - b.exponent();
    const double THRESHOLD = 2*FRACTION_RESCALING_THRESHOLD_LOG2;
    const int32_t e_delta_sat = int32_t(std::min(std::max(-THRESHOLD, e_delta),
                                                 THRESHOLD));
    const int32_t e_delta_idx = (e_delta_sat >> extended_double::FRACTION_RESCALING_THRESHOLD_LOG2_LOG2);
    ED_ASSERT_NORMALIZATION(e_delta_idx >= -2);
    ED_ASSERT_NORMALIZATION(e_delta_idx <= 2);
    const rescale_factors& f = rescale_factors_table[e_delta_idx + 2];
    
    a_f = a.fraction() * f.a_fraction_f;
    b_f = b.fraction() * f.b_fraction_f;
    
    ieee754_double_t e_a, e_b, e_max;
    e_a.as_double = a.m_exponent_raw;
    e_b.as_double = b.m_exponent_raw;
    e_max.as_uint64 = ((e_a.as_uint64 & f.a_exponent_mask)
                       | (e_b.as_uint64 & f.b_exponent_mask));
    e_raw = e_max.as_double;
}

#else

ED_ALWAYS_INLINE
void
extended_double::rescale_fractions(const extended_double& a, const extended_double& b,
                                   double& a_f, double& b_f, double& e_raw)
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
    const __m128d m =  _mm_load_sd(reinterpret_cast<const double*>(&EXPONENT_MASK));
    e_raw = _mm_cvtsd_f64(_mm_xor_pd(res_e_ab, m));
    
    /* Compute s =~ min(2^d, 0) by distinguishing the cases d=0, d=-1, d<=-2.
     * For d <= -2, it suffices to set s = 0, since then a + b = max(a, b).
     */
    const __m128i s = _mm_or_pd(_mm_or_pd(_mm_setzero_pd(),
                                          _mm_andnot_pd(d_lt0, _mm_set1_pd(1.0))),
                                _mm_and_pd(d_eqth, _mm_set1_pd(TH_INV)));
    
    /* Compute output fractions by multiplying with s */
    const __m128d f_ab = _mm_set_pd(b.fraction(), a.fraction());
    const __m128d res_f_ab = _mm_mul_pd(f_ab, s);
    a_f = _mm_cvtsd_f64(res_f_ab);
    b_f = _mm_cvtsd_f64(_mm_unpackhi_pd(res_f_ab, res_f_ab));
}

#endif // ED_ENABLE_SSE

void
extended_double::make_exponents_uniform_slowpath(extended_double& a, extended_double& b)
{
    double a_f, b_f, e_raw;
    rescale_fractions(a, b, a_f, b_f, e_raw);
    a.set_fraction(a_f);
    b.set_fraction(b_f);
    a.m_exponent_raw = b.m_exponent_raw = e_raw;
}

void
extended_double::add_slowpath(const extended_double& v)
{
    double t_f, v_f, e_raw;
    rescale_fractions(*this, v, t_f, v_f, e_raw);
    m_exponent_raw = e_raw;
    set_fraction(t_f + v_f);
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
