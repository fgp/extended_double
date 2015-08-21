#include "extended_double.h"

const int64_t extended_double::EXPONENT_EXCESS;

const int32_t extended_double::FRACTION_RESCALING_THRESHOLD_LOG2;

const double extended_double::FRACTION_RESCALING_THRESHOLD =
std::ldexp(1.0, extended_double::FRACTION_RESCALING_THRESHOLD_LOG2);

const double extended_double::LOG2 = std::log(2.0);

const uint32_t extended_double::IEEE754_DOUBLE_EXP_EXCESS;

const uint32_t extended_double::IEEE754_DOUBLE_EXP_INF_RAW;

const uint32_t extended_double::IEEE754_DOUBLE_EXP_MASK;

const uint32_t extended_double::IEEE754_DOUBLE_EXP_BITS;

const uint32_t extended_double::IEEE754_DOUBLE_MAN_BITS;

double extended_double::convert_to_double() const {
	/* Make fields of native double "m_fraction" available */
	ieee754_double_t v;
	v.as_double = m_fraction;

	/* Add exponent to native exponent, and clamp to valid IEEE754 double range */
	const int32_t e_sat = int32_t(std::max(std::min(exponent()
													+ int64_t(v.as_fields.exponent)
													- int64_t(IEEE754_DOUBLE_EXP_EXCESS),
													int64_t(IEEE754_DOUBLE_EXP_EXCESS + 1)),
										   -int64_t(IEEE754_DOUBLE_EXP_EXCESS)));

	/* Update native's exponent and return */
	v.as_fields.exponent = uint32_t(e_sat + int32_t(IEEE754_DOUBLE_EXP_EXCESS));
	return v.as_double;
}

void extended_double::normalize_after_multiply_slowpath() {
    /* Make fields of native double "m_fraction" available */
    ieee754_double_t f;
    f.as_double = m_fraction;
    const uint32_t f_e = f.as_fields.exponent;
    
    /* We only get here if the native exponent after a multiplication of two
     * normalized numbers exeends FRACTION_RESCALING_THRESHOLD_LOG2. This,
     * in particular, means that the fraction cannot be zero. In can be
     * infinity, though, if one of the factors was infinity. The code below
     * depends on this by
     *   a) Not checking for underflows of m_exponent_raw, since we only
     *      ever *increase it*. It must check for overflows, though!
     *   b) Not checking whether the native exponent is -EXPONENT_EXCESS,
     *      i.e. whether m_fraction is zero.
     * It seems product to verify these assumptions with an assert, though.
     */
    ED_ASSERT_NORMALIZATION(m_exponent_raw > 0);
    ED_ASSERT_NORMALIZATION(f_e >= (IEEE754_DOUBLE_EXP_EXCESS
                                    + FRACTION_RESCALING_THRESHOLD_LOG2));
    
	if (ED_LIKELY(f_e != IEEE754_DOUBLE_EXP_INF_RAW)) {
		/* Fraction is neither +/- Infinity nor NaN */

        /* Update exponent */
		m_exponent_raw += FRACTION_RESCALING_THRESHOLD_LOG2;

		if (ED_LIKELY(m_exponent_raw < EXPONENT_EXCESS + EXPONENT_MAX)) {
			/* Exponent still valid. Update fraction to reduced exponent */
			f.as_fields.exponent = f_e - FRACTION_RESCALING_THRESHOLD_LOG2;
			m_fraction = f.as_double;
		}
		else
            exponent_overflowed();
	}
	else {
		/* Fraction has value +/- Infinity or NaN. Adjust exponent accordingly */
		m_exponent_raw = EXPONENT_EXCESS + EXPONENT_INF;
	}
}

void extended_double::normalize_slowpath() {
    /* Make fields of native double "m_fraction" available */
    ieee754_double_t f;
    f.as_double = m_fraction;
    const uint32_t f_e_raw = f.as_fields.exponent;
    
    /* The following is the "full" version of normalize_after_multiply_slowpath().
     * It only assumes that m_exponent_raw is not *already* overflowed. The code
     * treast all de-normalized fractions (in the IEEE754 sense - not to be
     * confused with de-normalized instances of extended_double) as zero.
     */
    
    if (ED_LIKELY(f_e_raw != IEEE754_DOUBLE_EXP_INF_RAW)) {
        /* Fraction is neither +/- Infinity nor NaN */

        /* Compute mask which is all-zero if m_fraction is zero, otherwise all-one */
        const int64_t f_e_mask = -int32_t(f_e_raw > 0);

        /* Compute native exponent mod FRACTION_RESCALING_THRESHOLD_LOG2.
         * Result lies within [ 0, FRACTION_RESCALING_THRESHOLD_LOG2 ). This
         * will become the new native exponent (after zero masking)
         */
        const int32_t f_e_n = (f_e_raw + 1) % FRACTION_RESCALING_THRESHOLD_LOG2;
        
        /* Compute updated explicit exponent. If the fraction is zero, the explicit
         * exponent becomes zero as well (via masking). Otherwise, simply add
         * the delta between the original (f_e) and the normalized (f_e_n)
         * native exponent.
         */
        const int32_t f_e = int32_t(f_e_raw) - int32_t(IEEE754_DOUBLE_EXP_EXCESS);
        const int32_t e_delta = f_e - f_e_n;
        const int64_t e_p = (m_exponent_raw + e_delta) & f_e_mask;

		const int64_t e_p_test = e_p | (EXPONENT_EXCESS & ~f_e_mask);
		if (ED_LIKELY((e_p_test >= (EXPONENT_EXCESS + EXPONENT_MIN))
				      && (e_p_test <= (EXPONENT_EXCESS + EXPONENT_MAX))))
		{
            /* Exponent still valid. Update exponent and fraction.
             * If the fraction was zero (or denormalized), masking ensures that
             * it is set to zero.
             */
            m_exponent_raw = e_p;
            f.as_fields.exponent = f_e_n + int32_t(IEEE754_DOUBLE_EXP_EXCESS);
            f.as_uint64 &= f_e_mask | 0x8000000000000000;
            m_fraction = f.as_double;
        }
        else
            exponent_overflowed();
    }
    else {
        /* Fraction has value +/- Infinity or NaN. Adjust exponent accordingly */
        m_exponent_raw = EXPONENT_EXCESS + EXPONENT_INF;
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
	a_factor.as_uint64 = (uint64_t(a_shift_sat + int64_t(IEEE754_DOUBLE_EXP_EXCESS))
						  << IEEE754_DOUBLE_MAN_BITS);
	b_factor.as_uint64 = (uint64_t(b_shift_sat + int64_t(IEEE754_DOUBLE_EXP_EXCESS))
						  << IEEE754_DOUBLE_MAN_BITS);

	/* Multiple a and b by 2^shift_a resp. 2^shift_b, and adjust the exponents accordingly.
     * Since either a_shift = a.exp - b.exp, or b_shift = b.exp - a.exp, the exponents
     * of a and b will afterwards agree.
     */
	a.m_fraction *= a_factor.as_double;
	b.m_fraction *= b_factor.as_double;
	a.m_exponent_raw -= a_shift;
	b.m_exponent_raw -= b_shift;
}

void extended_double::exponent_overflowed() {
	__builtin_trap();
}

std::ostream& operator<<(std::ostream& dst, const extended_double& v) {
	if (v.m_exponent_raw == 0)
		dst << '0';
	else {
		int e = 0;
		const double f = std::frexp(v.fraction(), &e);
		dst << f << "*2^" << (v.exponent() + static_cast<int64_t>(e));
	}
    
    return dst;
}
