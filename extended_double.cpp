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

	if (ED_LIKELY(f_e != IEEE754_DOUBLE_EXP_INF_RAW)) {
		/* Fraction is neither +/- Infinity nor NaN */

		/* If the fraction increased beyond the rescaling threshold, increase
		 * the exponent and rescale the fraction. Note after multiplication,
		 * the fraction is certainly less than FRACTION_RESCALING_THRESHOLD ^ 2,
		 * so a single rescaling step always suffices here.
		 */
		const bool rescale = (f_e >= (FRACTION_RESCALING_THRESHOLD_LOG2
									  + IEEE754_DOUBLE_EXP_EXCESS));
		m_exponent_raw += rescale * FRACTION_RESCALING_THRESHOLD_LOG2;

		if (ED_LIKELY(m_exponent_raw < EXPONENT_EXCESS + EXPONENT_MAX)) {
			/* Exponent still valid. Update fraction to reduced exponent */
			f.as_fields.exponent = f_e - FRACTION_RESCALING_THRESHOLD_LOG2;
			m_fraction = f.as_double;
		}
		else {
			/* Exponent overflowed. Set value to +/- Infinity */
			m_exponent_raw = EXPONENT_EXCESS + EXPONENT_MAX;
			f.as_fields.exponent = IEEE754_DOUBLE_EXP_INF_RAW;
			f.as_fields.mantissa = 0;
			m_fraction = f.as_double;
		}
	}
	else {
		/* Fraction has value +/- Infinity or NaN. Adjust exponent accordingly */
		m_exponent_raw = EXPONENT_EXCESS + EXPONENT_MAX;
	}
}

void extended_double::normalize_slowpath() {
	/* Make fields of native double "m_fraction" available */
	ieee754_double_t v;
	v.as_double = m_fraction;

	/* Extract native exponent and determine if v is non-finite or sub-normal.
     * For sub-normal values, e_type is 0, while for inifinities its 1.
     * Note that zero is counted as being sub-normal here - as other sub-normal
     * numbers, its native exponent is the smallest possible value.
     */
	const uint32_t e = uint32_t(v.as_fields.exponent);
	const bool e_sub = (e == 0);
	const bool e_inf = (e == IEEE754_DOUBLE_EXP_MASK);

	/* Update exponent field. For normals, the fraction's native exponent
     * is added, rounded to a multiple of  FRACTION_RESCALING_THRESHOLD_LOG2.
     * For sub-normals (incl. zeros) the exponent is set to the smallest value,
     * while for infinities it is set to the largest.
     */
	const int32_t e_keep = (e + 1) % FRACTION_RESCALING_THRESHOLD_LOG2;
	const int32_t e_val = int32_t(e) - int32_t(IEEE754_DOUBLE_EXP_EXCESS);
	const int32_t e_delta = e_val - e_keep;
	const int64_t e_adj = m_exponent_raw + e_delta;
	m_exponent_raw = e_sub ? 0 : e_inf ? EXPONENT_EXCESS + EXPONENT_MAX : e_adj;

	/* Update m_fraction's native exponent and mantissa. For non-finite
     * and sub-normal values, the original exponent is kept - this preserves
     * Zero, +/- Infinity and NaN. For normal values, the remainder of the
     * exponent mod FRACTION_RESCALING_THRESHOLD_LOG2 is subtracted, since
     * that part of the exponent was "moved" to the exponent field above.
     * Sub-normal values are rounded down to zero by zeroing the mantissa,
     * to guarantee that the result is normalized.  Note that doing this for
     * all non-finite values would convert NaN to +/- Inf!
     */
	const uint32_t sp = uint32_t(v.as_fields.sign);
	const uint32_t ep = (e_sub || e_inf) ? e : uint32_t(e_keep + int32_t(IEEE754_DOUBLE_EXP_EXCESS));
	const uint64_t mp = e_sub ? 0 : v.as_fields.mantissa;
	ieee754_double_t vp;
	vp.as_uint64 = ((uint64_t(sp) << (IEEE754_DOUBLE_MAN_BITS + IEEE754_DOUBLE_EXP_BITS)) |
	 			    (uint64_t(ep) << IEEE754_DOUBLE_MAN_BITS) |
			        mp);
	m_fraction = vp.as_double;
}

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
