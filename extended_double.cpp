#include "extended_double.h"

const int64_t extended_double::EXPONENT_EXCESS;

const int64_t extended_double::FRACTION_RESCALING_THRESHOLD_LOG2;

const double extended_double::FRACTION_RESCALING_THRESHOLD =
std::ldexp(1.0, extended_double::FRACTION_RESCALING_THRESHOLD_LOG2);

const double extended_double::FRACTION_RESCALING_THRESHOLD_INV =
std::ldexp(1.0, -extended_double::FRACTION_RESCALING_THRESHOLD_LOG2);

const double extended_double::LOG2 =
log(2.0);

std::ostream& operator<<(std::ostream& dst, const extended_double& v) {
	if (v.m_exponent_raw == 0)
		dst << '0';
	else {
		int e = 0;
		const double f = std::frexp(v.fraction(), &e);
		dst << f << "*2^" << (v.exponent() + static_cast<int64_t>(e));
	}
}
