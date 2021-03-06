/*
 *  tests.cpp
 *  extended_double
 *
 *  Created by Florian Pflug on 17.08.2015.
 *  Copyright (c) 2015 Florian Pflug. All rights reserved.
 */

#include <boost/test/unit_test.hpp>

#include "extended_double.h"

BOOST_AUTO_TEST_SUITE(extended_double_tests)

namespace {
	double nextafter(double v) {
		union {
			double as_double;
			uint64_t as_uint64;
		} t = { v };
		++t.as_uint64;
		return t.as_double;
	}

	double firstbefore(double v) {
		union {
			double as_double;
			uint64_t as_uint64;
		} t = { v };
		--t.as_uint64;
		return t.as_double;
	}

	template<typename T>
	bool is_identical(const T& a, const T& b) {
		return !std::memcmp(&a, &b, sizeof(T));
	}

	template<typename T>
	bool is_not_identical(const T& a, const T& b) {
		return std::memcmp(&a, &b, sizeof(T));
	}
}

#define BOOST_CHECK_ED_IDENTICAL(x, y) \
BOOST_CHECK_PREDICATE(is_identical<extended_double>, (x)(y))
#define BOOST_CHECK_ED_NOT_IDENTICAL(x, y) \
BOOST_CHECK_PREDICATE(is_not_identical<extended_double>, (x)(y))

BOOST_AUTO_TEST_CASE(asserts_enabled_normalization) {
	BOOST_CHECK_EQUAL(ED_ENABLE_ASSERTS_NORMALIZATION, 1);
#ifdef NDEBUG
	const bool ndebug = true;
#else
	const bool ndebug = false;
#endif
	BOOST_CHECK_EQUAL(ndebug, false);
}

BOOST_AUTO_TEST_CASE(asserts_enabled_static) {
	BOOST_CHECK_EQUAL(ED_ENABLE_ASSERTS_STATIC, 1);
}

BOOST_AUTO_TEST_CASE(print_sse_status) {
	std::cerr << "Using SSE: ";
#if ED_ENABLE_SSE
	std::cerr << "Yes";
#else
	std::cerr << "No";
#endif
	std::cerr << std::endl;
}

BOOST_AUTO_TEST_CASE(conversions) {
	const double fractions[] = { -firstbefore(2.0), -1.9, -1.0 - 1.0/M_PI,
		-1.0, nextafter(1.0), 1.1, 1.0 + 1.0/M_PI, 1.9,
		firstbefore(2.0) };

	/* [ -1022, 1023 ] is the range the exponents of IEEE754 non-zero and finite
	 * double-precision values. -1023 is the exponent of zero (and sub-normals),
	 * and 1024 the exponent of +/- infinity.
	 */
	for(int32_t e = -1022; e <= 1023; ++e) {
		for(int f = 0; f < sizeof(fractions) / sizeof(double); ++f) {
			const double d = std::ldexp(fractions[f], e);
			BOOST_CHECK(d != 0.0);
			BOOST_CHECK(std::isfinite(d));

			const extended_double v = d;
			BOOST_CHECK_EQUAL(extended_double_cast<double>(v), d);
			BOOST_CHECK_EQUAL(int64_t(v.exponent())
							  % extended_double::FRACTION_RESCALING_THRESHOLD_LOG2,
							  0);

			int v_f_e = 0;
			const double v_fp = std::frexp(v.fraction(), &v_f_e);
			BOOST_CHECK_GE(std::fabs(v_fp), 0.5);
			BOOST_CHECK_LT(std::fabs(v_fp), 1.0);
			v_f_e -= 1;

			BOOST_CHECK_GE(v_f_e, 0);
			BOOST_CHECK_LT(v_f_e, extended_double::FRACTION_RESCALING_THRESHOLD_LOG2);
			BOOST_CHECK_EQUAL(v.exponent() + v_f_e, e);
		}
	}
}

BOOST_AUTO_TEST_CASE(zeros) {
	extended_double v0;
	BOOST_CHECK_EQUAL(v0, extended_double(0.0));
	BOOST_CHECK_EQUAL(v0.fraction(), 0.0);

	std::memset(&v0, 0, sizeof(v0));
	BOOST_CHECK_EQUAL(v0, extended_double(0.0));
	BOOST_CHECK_EQUAL(v0.fraction(), 0.0);

	extended_double v0m(-0.0);
	BOOST_CHECK_EQUAL(v0m, v0);
	BOOST_CHECK_ED_NOT_IDENTICAL(v0, v0m);
	BOOST_CHECK_ED_IDENTICAL(-v0, v0m);
	BOOST_CHECK_ED_IDENTICAL(-v0m, v0);

	BOOST_CHECK_ED_IDENTICAL(v0 + v0, v0);
	BOOST_CHECK_ED_IDENTICAL(v0 + v0m, v0);
	BOOST_CHECK_ED_IDENTICAL(v0m + v0, v0);

	BOOST_CHECK_ED_IDENTICAL(v0m + v0m, v0m);
	BOOST_CHECK_ED_IDENTICAL(v0m - v0, v0m);

	const double fractions[] = { -firstbefore(2.0), -1.9, -1.0 - 1.0/M_PI,
		-1.0, nextafter(1.0), 1.1, 1.0 + 1.0/M_PI, 1.9,
		firstbefore(2.0) };

	for(int32_t e = -2048; e <= 2048; ++e) {
		for(int f = 0; f < sizeof(fractions) / sizeof(double); ++f) {
			const extended_double v = fractions[f] * extended_double::pow2(e);
			BOOST_CHECK_NE(v, v0);
			BOOST_CHECK_NE(v, v0m);
			BOOST_CHECK_EQUAL(v / fractions[f], extended_double::pow2(e));

			/* Operations with positive zero */

			BOOST_CHECK_EQUAL(v * v0, v0);
			BOOST_CHECK_EQUAL(v * v0, v0m);
			BOOST_CHECK_ED_IDENTICAL(v * v0, (v >= 0.0) ? v0 : v0m);

			BOOST_CHECK_EQUAL(v / v0, fractions[f] * std::numeric_limits<double>::infinity());
			BOOST_CHECK_EQUAL(v + v0, v);
			BOOST_CHECK_EQUAL(v - v0, v);

			BOOST_CHECK_EQUAL(v0 * v, v0);
			BOOST_CHECK_EQUAL(v0 * v, v0m);
			BOOST_CHECK_ED_IDENTICAL(v0 * v, (v >= 0.0) ? v0 : v0m);

			BOOST_CHECK_EQUAL(v0 / v, v0);
			BOOST_CHECK_EQUAL(v0 / v, v0m);
			BOOST_CHECK_ED_IDENTICAL(v0 / v, (v >= 0.0) ? v0 : v0m);

			BOOST_CHECK_EQUAL(v0 + v, v);
			BOOST_CHECK_EQUAL(v0 - v, -v);

			/* Operations with negative zero */

			BOOST_CHECK_EQUAL(v * v0m, v0);
			BOOST_CHECK_EQUAL(v * v0m, v0m);
			BOOST_CHECK_ED_IDENTICAL(v * v0m, (v >= 0.0) ? v0m : v0);

			BOOST_CHECK_EQUAL(v / v0m, -fractions[f] * std::numeric_limits<double>::infinity());
			BOOST_CHECK_EQUAL(v + v0m, v);
			BOOST_CHECK_EQUAL(v - v0m, v);

			BOOST_CHECK_EQUAL(v0m * v, v0);
			BOOST_CHECK_EQUAL(v0m * v, v0m);
			BOOST_CHECK_ED_IDENTICAL(v0m * v, (v >= 0.0) ? v0m : v0);

			BOOST_CHECK_EQUAL(v0m / v, v0);
			BOOST_CHECK_EQUAL(v0m / v, v0m);
			BOOST_CHECK_ED_IDENTICAL(v0m / v, (v >= 0.0) ? v0m : v0);

			BOOST_CHECK_EQUAL(v0m + v, v);
			BOOST_CHECK_EQUAL(v0m - v, -v);
		}
	}
}

BOOST_AUTO_TEST_CASE(infinities) {
	const extended_double v_pinf(std::numeric_limits<double>::infinity());
	BOOST_CHECK_EQUAL(extended_double_cast<double>(v_pinf),
					  std::numeric_limits<double>::infinity());
	BOOST_CHECK(!isfinite(v_pinf));

	const extended_double v_ninf(-std::numeric_limits<double>::infinity());
	BOOST_CHECK_EQUAL(extended_double_cast<double>(v_ninf),
					  -std::numeric_limits<double>::infinity());
	BOOST_CHECK(!isfinite(v_ninf));

	const extended_double v_small = extended_double::pow2(-1024);
	BOOST_CHECK_NE(v_small, 0.0);
	BOOST_CHECK_EQUAL(log2(v_small), -1024);

	const extended_double v_large = extended_double::pow2(1024);
	BOOST_CHECK_EQUAL(log2(v_large), 1024);


	BOOST_CHECK_EQUAL(v_pinf, v_pinf);
	BOOST_CHECK_EQUAL(v_ninf, v_ninf);
	BOOST_CHECK_EQUAL(v_pinf, -v_ninf);
	BOOST_CHECK_EQUAL(-v_pinf, v_ninf);

	BOOST_CHECK_LT(v_ninf, v_pinf);
	BOOST_CHECK_LE(v_ninf, v_pinf);
	BOOST_CHECK_NE(v_ninf, v_pinf);
	BOOST_CHECK_GE(v_pinf, v_ninf);
	BOOST_CHECK_GT(v_pinf, v_ninf);

	BOOST_CHECK_LT(v_ninf, -v_large);
	BOOST_CHECK_LT(v_ninf, -v_small);
	BOOST_CHECK_LT(v_ninf, v_small);
	BOOST_CHECK_LT(v_ninf, v_large);
	BOOST_CHECK_GT(v_pinf, -v_large);
	BOOST_CHECK_GT(v_pinf, -v_small);
	BOOST_CHECK_GT(v_pinf, v_small);
	BOOST_CHECK_GT(v_pinf, v_large);

	BOOST_CHECK_EQUAL(v_pinf, v_pinf + v_small);
	BOOST_CHECK_EQUAL(v_pinf, v_pinf - v_small);
	BOOST_CHECK_EQUAL(v_pinf, v_pinf + 1e-100);
	BOOST_CHECK_EQUAL(v_pinf, v_pinf - 1e-100);
	BOOST_CHECK_EQUAL(v_pinf, v_pinf + 1.0);
	BOOST_CHECK_EQUAL(v_pinf, v_pinf - 1.0);
	BOOST_CHECK_EQUAL(v_pinf, v_pinf + 1e100);
	BOOST_CHECK_EQUAL(v_pinf, v_pinf - 1e100);
	BOOST_CHECK_EQUAL(v_pinf, v_pinf + v_large);
	BOOST_CHECK_EQUAL(v_pinf, v_pinf - v_large);
	BOOST_CHECK_EQUAL(v_pinf, v_pinf + v_pinf);
	BOOST_CHECK_EQUAL(v_pinf, v_pinf - v_ninf);

	BOOST_CHECK_EQUAL(v_pinf, v_pinf * v_small);
	BOOST_CHECK_EQUAL(v_pinf, v_pinf * 1e-100);
	BOOST_CHECK_EQUAL(v_pinf, v_pinf * 1.0);
	BOOST_CHECK_EQUAL(v_pinf, v_pinf * 1e100);
	BOOST_CHECK_EQUAL(v_pinf, v_pinf * v_large);
	BOOST_CHECK_EQUAL(v_pinf, v_pinf * v_pinf);
	BOOST_CHECK_EQUAL(v_ninf, v_ninf * v_small);
	BOOST_CHECK_EQUAL(v_ninf, v_ninf * 1e-100);
	BOOST_CHECK_EQUAL(v_ninf, v_ninf * 1.0);
	BOOST_CHECK_EQUAL(v_ninf, v_ninf * 1e100);
	BOOST_CHECK_EQUAL(v_ninf, v_ninf * v_large);
	BOOST_CHECK_EQUAL(v_ninf, v_ninf * v_pinf);

	BOOST_CHECK_EQUAL(v_ninf, v_pinf * -v_small);
	BOOST_CHECK_EQUAL(v_ninf, v_pinf * -1e-100);
	BOOST_CHECK_EQUAL(v_ninf, v_pinf * -1.0);
	BOOST_CHECK_EQUAL(v_ninf, v_pinf * -1e100);
	BOOST_CHECK_EQUAL(v_ninf, v_pinf * -v_large);
	BOOST_CHECK_EQUAL(v_ninf, v_pinf * v_ninf);
	BOOST_CHECK_EQUAL(v_pinf, v_ninf * v_ninf);
	BOOST_CHECK_EQUAL(v_pinf, v_ninf * -v_small);
	BOOST_CHECK_EQUAL(v_pinf, v_ninf * -1e-100);
	BOOST_CHECK_EQUAL(v_pinf, v_ninf * -1.0);
	BOOST_CHECK_EQUAL(v_pinf, v_ninf * -1e100);
	BOOST_CHECK_EQUAL(v_pinf, v_ninf * -v_small);
	BOOST_CHECK_EQUAL(v_pinf, v_ninf * v_ninf);

	BOOST_CHECK_EQUAL(v_pinf, v_pinf / v_small);
	BOOST_CHECK_EQUAL(v_pinf, v_pinf / 1e-100);
	BOOST_CHECK_EQUAL(v_pinf, v_pinf / 1.0);
	BOOST_CHECK_EQUAL(v_pinf, v_pinf / 1e100);
	BOOST_CHECK_EQUAL(v_pinf, v_pinf / v_large);
	BOOST_CHECK_EQUAL(v_ninf, v_ninf / 1e-100);
	BOOST_CHECK_EQUAL(v_ninf, v_ninf / v_small);
	BOOST_CHECK_EQUAL(v_ninf, v_ninf / 1.0);
	BOOST_CHECK_EQUAL(v_ninf, v_ninf / 1e100);
	BOOST_CHECK_EQUAL(v_ninf, v_ninf / v_small);

	BOOST_CHECK_EQUAL(v_ninf, v_pinf / -v_small);
	BOOST_CHECK_EQUAL(v_ninf, v_pinf / -1e-100);
	BOOST_CHECK_EQUAL(v_ninf, v_pinf / -1.0);
	BOOST_CHECK_EQUAL(v_ninf, v_pinf / -1e100);
	BOOST_CHECK_EQUAL(v_ninf, v_pinf / -v_large);
	BOOST_CHECK_EQUAL(v_pinf, v_ninf / -v_small);
	BOOST_CHECK_EQUAL(v_pinf, v_ninf / -1e-100);
	BOOST_CHECK_EQUAL(v_pinf, v_ninf / -1.0);
	BOOST_CHECK_EQUAL(v_pinf, v_ninf / -1e100);
	BOOST_CHECK_EQUAL(v_pinf, v_ninf / -v_large);

	BOOST_CHECK_EQUAL(v_pinf, v_small + v_pinf);
	BOOST_CHECK_EQUAL(-v_pinf, v_small - v_pinf);
	BOOST_CHECK_EQUAL(v_pinf, 1e-100 + v_pinf);
	BOOST_CHECK_EQUAL(-v_pinf, 1e-100 - v_pinf);
	BOOST_CHECK_EQUAL(v_pinf, 1.0 + v_pinf);
	BOOST_CHECK_EQUAL(-v_pinf, 1.0 - v_pinf);
	BOOST_CHECK_EQUAL(v_pinf, 1e100 + v_pinf);
	BOOST_CHECK_EQUAL(-v_pinf, 1e100 - v_pinf);
	BOOST_CHECK_EQUAL(v_pinf, v_large + v_pinf);
	BOOST_CHECK_EQUAL(-v_pinf, v_large - v_pinf);
	BOOST_CHECK_EQUAL(v_pinf, v_pinf + v_pinf);
	BOOST_CHECK_EQUAL(-v_pinf, v_ninf - v_pinf);

	BOOST_CHECK_EQUAL(v_pinf, v_small * v_pinf);
	BOOST_CHECK_EQUAL(v_pinf, 1e-100 * v_pinf);
	BOOST_CHECK_EQUAL(v_pinf, 1.0 * v_pinf);
	BOOST_CHECK_EQUAL(v_pinf, 1e100 * v_pinf);
	BOOST_CHECK_EQUAL(v_pinf, v_large * v_pinf);
	BOOST_CHECK_EQUAL(v_pinf, v_pinf * v_pinf);
	BOOST_CHECK_EQUAL(v_ninf, v_small * v_ninf);
	BOOST_CHECK_EQUAL(v_ninf, 1e-100 * v_ninf);
	BOOST_CHECK_EQUAL(v_ninf, 1.0 * v_ninf);
	BOOST_CHECK_EQUAL(v_ninf, 1e100 * v_ninf);
	BOOST_CHECK_EQUAL(v_ninf, v_large * v_ninf);
	BOOST_CHECK_EQUAL(v_ninf, v_pinf * v_ninf);

	BOOST_CHECK_EQUAL(v_ninf, -v_small * v_pinf);
	BOOST_CHECK_EQUAL(v_ninf, -1e-100 * v_pinf);
	BOOST_CHECK_EQUAL(v_ninf, -1.0 * v_pinf);
	BOOST_CHECK_EQUAL(v_ninf, -1e100 * v_pinf);
	BOOST_CHECK_EQUAL(v_ninf, -v_large * v_pinf);
	BOOST_CHECK_EQUAL(v_ninf, v_ninf * v_pinf);
	BOOST_CHECK_EQUAL(v_pinf, v_ninf * v_ninf);
	BOOST_CHECK_EQUAL(v_pinf, -v_small * v_ninf);
	BOOST_CHECK_EQUAL(v_pinf, -1e-100 * v_ninf);
	BOOST_CHECK_EQUAL(v_pinf, -1.0 * v_ninf);
	BOOST_CHECK_EQUAL(v_pinf, -1e100 * v_ninf);
	BOOST_CHECK_EQUAL(v_pinf, -v_small * v_ninf);
	BOOST_CHECK_EQUAL(v_pinf, v_ninf * v_ninf);

	BOOST_CHECK_EQUAL(0.0, v_small / v_pinf);
	BOOST_CHECK_EQUAL(0.0, 1e-100 / v_pinf);
	BOOST_CHECK_EQUAL(0.0, 1.0 / v_pinf);
	BOOST_CHECK_EQUAL(0.0, 1e100 / v_pinf);
	BOOST_CHECK_EQUAL(0.0, v_large / v_pinf);
	BOOST_CHECK_EQUAL(0.0, 1e-100 / v_ninf);
	BOOST_CHECK_EQUAL(0.0, v_small / v_ninf);
	BOOST_CHECK_EQUAL(0.0, 1.0 / v_ninf);
	BOOST_CHECK_EQUAL(0.0, 1e100 / v_ninf);
	BOOST_CHECK_EQUAL(0.0, v_small / v_ninf);

	BOOST_CHECK_EQUAL(0.0, -v_small / v_pinf);
	BOOST_CHECK_EQUAL(0.0, -1e-100 / v_pinf);
	BOOST_CHECK_EQUAL(0.0, -1.0 / v_pinf);
	BOOST_CHECK_EQUAL(0.0, -1e100 / v_pinf);
	BOOST_CHECK_EQUAL(0.0, -v_large / v_pinf);
	BOOST_CHECK_EQUAL(0.0, -v_small / v_ninf);
	BOOST_CHECK_EQUAL(0.0, -1e-100 / v_ninf);
	BOOST_CHECK_EQUAL(0.0, -1.0 / v_ninf);
	BOOST_CHECK_EQUAL(0.0, -1e100 / v_ninf);
	BOOST_CHECK_EQUAL(0.0, -v_large / v_ninf);
}

BOOST_AUTO_TEST_CASE(nans) {
	const extended_double v_qnan(std::numeric_limits<double>::quiet_NaN());
	const extended_double v_snan(std::numeric_limits<double>::signaling_NaN());
	BOOST_CHECK(std::isnan(extended_double_cast<double>(v_qnan)));
	BOOST_CHECK(std::isnan(extended_double_cast<double>(v_snan)));

	/* Mixing NaNs */

	BOOST_CHECK(isnan(v_qnan + v_qnan));
	BOOST_CHECK(isnan(v_qnan - v_qnan));
	BOOST_CHECK(isnan(v_qnan * v_qnan));
	BOOST_CHECK(isnan(v_qnan / v_qnan));

	BOOST_CHECK(isnan(v_snan + v_qnan));
	BOOST_CHECK(isnan(v_snan - v_qnan));
	BOOST_CHECK(isnan(v_snan * v_qnan));
	BOOST_CHECK(isnan(v_snan / v_qnan));

	BOOST_CHECK(isnan(v_qnan + v_snan));
	BOOST_CHECK(isnan(v_qnan - v_snan));
	BOOST_CHECK(isnan(v_qnan * v_snan));
	BOOST_CHECK(isnan(v_qnan / v_snan));

	BOOST_CHECK(isnan(v_snan + v_snan));
	BOOST_CHECK(isnan(v_snan - v_snan));
	BOOST_CHECK(isnan(v_snan * v_snan));
	BOOST_CHECK(isnan(v_snan / v_snan));

	const extended_double v_small = extended_double::pow2(-1024);
	BOOST_CHECK_NE(v_small, 0.0);
	BOOST_CHECK_EQUAL(log2(v_small), -1024);

	const extended_double v_large = extended_double::pow2(1024);
	BOOST_CHECK_EQUAL(log2(v_large), 1024);

	/* Quiet NaN */

	BOOST_CHECK(isnan(v_qnan + v_small));
	BOOST_CHECK(isnan(v_qnan - v_small));
	BOOST_CHECK(isnan(v_qnan + 1e-100));
	BOOST_CHECK(isnan(v_qnan - 1e-100));
	BOOST_CHECK(isnan(v_qnan + 1.0));
	BOOST_CHECK(isnan(v_qnan - 1.0));
	BOOST_CHECK(isnan(v_qnan + 1e100));
	BOOST_CHECK(isnan(v_qnan - 1e100));
	BOOST_CHECK(isnan(v_qnan + v_large));
	BOOST_CHECK(isnan(v_qnan - v_large));

	BOOST_CHECK(isnan(v_qnan * v_small));
	BOOST_CHECK(isnan(v_qnan * 1e-100));
	BOOST_CHECK(isnan(v_qnan * 1.0));
	BOOST_CHECK(isnan(v_qnan * 1e100));
	BOOST_CHECK(isnan(v_qnan * v_large));
	BOOST_CHECK(isnan(v_qnan * -v_small));
	BOOST_CHECK(isnan(v_qnan * -1e-100));
	BOOST_CHECK(isnan(v_qnan * -1.0));
	BOOST_CHECK(isnan(v_qnan * -1e100));
	BOOST_CHECK(isnan(v_qnan * -v_large));

	BOOST_CHECK(isnan(v_qnan / v_small));
	BOOST_CHECK(isnan(v_qnan / 1e-100));
	BOOST_CHECK(isnan(v_qnan / 1.0));
	BOOST_CHECK(isnan(v_qnan / 1e100));
	BOOST_CHECK(isnan(v_qnan / v_large));
	BOOST_CHECK(isnan(v_qnan / -v_small));
	BOOST_CHECK(isnan(v_qnan / -1e-100));
	BOOST_CHECK(isnan(v_qnan / -1.0));
	BOOST_CHECK(isnan(v_qnan / -1e100));
	BOOST_CHECK(isnan(v_qnan / -v_large));

	BOOST_CHECK(isnan(v_small + v_qnan));
	BOOST_CHECK(isnan(v_small - v_qnan));
	BOOST_CHECK(isnan(1e-100 + v_qnan));
	BOOST_CHECK(isnan(1e-100 - v_qnan));
	BOOST_CHECK(isnan(1.0 + v_qnan));
	BOOST_CHECK(isnan(1.0 - v_qnan));
	BOOST_CHECK(isnan(1e100 + v_qnan));
	BOOST_CHECK(isnan(1e100 - v_qnan));
	BOOST_CHECK(isnan(v_large + v_qnan));
	BOOST_CHECK(isnan(v_large - v_qnan));

	BOOST_CHECK(isnan(v_small * v_qnan));
	BOOST_CHECK(isnan(1e-100 * v_qnan));
	BOOST_CHECK(isnan(1.0 * v_qnan));
	BOOST_CHECK(isnan(1e100 * v_qnan));
	BOOST_CHECK(isnan(v_large * v_qnan));
	BOOST_CHECK(isnan(-v_small * v_qnan));
	BOOST_CHECK(isnan(-1e-100 * v_qnan));
	BOOST_CHECK(isnan(-1.0 * v_qnan));
	BOOST_CHECK(isnan(-1e100 * v_qnan));
	BOOST_CHECK(isnan(-v_large * v_qnan));

	BOOST_CHECK(isnan(v_small / v_qnan));
	BOOST_CHECK(isnan(1e-100 / v_qnan));
	BOOST_CHECK(isnan(1.0 / v_qnan));
	BOOST_CHECK(isnan(1e100 / v_qnan));
	BOOST_CHECK(isnan(v_large / v_qnan));
	BOOST_CHECK(isnan(-v_small / v_qnan));
	BOOST_CHECK(isnan(-1e-100 / v_qnan));
	BOOST_CHECK(isnan(-1.0 / v_qnan));
	BOOST_CHECK(isnan(-1e100 / v_qnan));
	BOOST_CHECK(isnan(-v_large / v_qnan));

	/* Signalling NaN */

	BOOST_CHECK(isnan(v_snan + v_small));
	BOOST_CHECK(isnan(v_snan - v_small));
	BOOST_CHECK(isnan(v_snan + 1e-100));
	BOOST_CHECK(isnan(v_snan - 1e-100));
	BOOST_CHECK(isnan(v_snan + 1.0));
	BOOST_CHECK(isnan(v_snan - 1.0));
	BOOST_CHECK(isnan(v_snan + 1e100));
	BOOST_CHECK(isnan(v_snan - 1e100));
	BOOST_CHECK(isnan(v_snan + v_large));
	BOOST_CHECK(isnan(v_snan - v_large));

	BOOST_CHECK(isnan(v_snan * v_small));
	BOOST_CHECK(isnan(v_snan * 1e-100));
	BOOST_CHECK(isnan(v_snan * 1.0));
	BOOST_CHECK(isnan(v_snan * 1e100));
	BOOST_CHECK(isnan(v_snan * v_large));
	BOOST_CHECK(isnan(v_snan * -v_small));
	BOOST_CHECK(isnan(v_snan * -1e-100));
	BOOST_CHECK(isnan(v_snan * -1.0));
	BOOST_CHECK(isnan(v_snan * -1e100));
	BOOST_CHECK(isnan(v_snan * -v_large));

	BOOST_CHECK(isnan(v_snan / v_small));
	BOOST_CHECK(isnan(v_snan / 1e-100));
	BOOST_CHECK(isnan(v_snan / 1.0));
	BOOST_CHECK(isnan(v_snan / 1e100));
	BOOST_CHECK(isnan(v_snan / v_large));
	BOOST_CHECK(isnan(v_snan / -v_small));
	BOOST_CHECK(isnan(v_snan / -1e-100));
	BOOST_CHECK(isnan(v_snan / -1.0));
	BOOST_CHECK(isnan(v_snan / -1e100));
	BOOST_CHECK(isnan(v_snan / -v_large));

	BOOST_CHECK(isnan(v_small + v_snan));
	BOOST_CHECK(isnan(v_small - v_snan));
	BOOST_CHECK(isnan(1e-100 + v_snan));
	BOOST_CHECK(isnan(1e-100 - v_snan));
	BOOST_CHECK(isnan(1.0 + v_snan));
	BOOST_CHECK(isnan(1.0 - v_snan));
	BOOST_CHECK(isnan(1e100 + v_snan));
	BOOST_CHECK(isnan(1e100 - v_snan));
	BOOST_CHECK(isnan(v_large + v_snan));
	BOOST_CHECK(isnan(v_large - v_snan));

	BOOST_CHECK(isnan(v_small * v_snan));
	BOOST_CHECK(isnan(1e-100 * v_snan));
	BOOST_CHECK(isnan(1.0 * v_snan));
	BOOST_CHECK(isnan(1e100 * v_snan));
	BOOST_CHECK(isnan(v_large * v_snan));
	BOOST_CHECK(isnan(-v_small * v_snan));
	BOOST_CHECK(isnan(-1e-100 * v_snan));
	BOOST_CHECK(isnan(-1.0 * v_snan));
	BOOST_CHECK(isnan(-1e100 * v_snan));
	BOOST_CHECK(isnan(-v_large * v_snan));

	BOOST_CHECK(isnan(v_small / v_snan));
	BOOST_CHECK(isnan(1e-100 / v_snan));
	BOOST_CHECK(isnan(1.0 / v_snan));
	BOOST_CHECK(isnan(1e100 / v_snan));
	BOOST_CHECK(isnan(v_large / v_snan));
	BOOST_CHECK(isnan(-v_small / v_snan));
	BOOST_CHECK(isnan(-1e-100 / v_snan));
	BOOST_CHECK(isnan(-1.0 / v_snan));
	BOOST_CHECK(isnan(-1e100 / v_snan));
	BOOST_CHECK(isnan(-v_large / v_snan));
}

BOOST_AUTO_TEST_CASE(arithmetic) {
	const int64_t exponents[] = { -0x100000000, -1500, -512, -511,
		-1, 0, 511, 512, 1500, 0x100000000 };
	const double fractions[] = {
		-firstbefore(extended_double::FRACTION_RESCALING_THRESHOLD),
		-extended_double::FRACTION_RESCALING_THRESHOLD/2.0,
		nextafter(2.0), -firstbefore(2.0), -1.0 - 1.0/M_PI, nextafter(-1.0),
		-1.0, firstbefore(-1.0), firstbefore(1.0), 1.0,
		nextafter(1.0), 1.0 + 1.0/M_PI, firstbefore(2.0), nextafter(2.0),
		extended_double::FRACTION_RESCALING_THRESHOLD/2.0,
		firstbefore(extended_double::FRACTION_RESCALING_THRESHOLD),
	};

	for(int e1_i = 0; e1_i < sizeof(exponents) / sizeof(int64_t); ++e1_i) {
		for(int e2_i = 0; e2_i < sizeof(exponents) / sizeof(int64_t); ++e2_i) {
			const int64_t e1_1 = exponents[e1_i];
			const int64_t e2_1 = exponents[e1_i];
			for(int e1_2 = 0; e1_2 <= 513; e1_2 += 57) {
				for(int e2_2 = 0; e2_2 <= 513; e2_2 += 57) {
					const int e1 = e1_1 + e1_2;
					const int e2 = e2_1 + e2_2;
					for(int f1 = 0; f1 < sizeof(fractions) / sizeof(double); ++f1) {
						for(int f2 = 0; f2 < sizeof(fractions) / sizeof(double); ++f2) {
							const extended_double v1 = fractions[f1] * extended_double::pow2(e1);
							const extended_double v2 = fractions[f2] * extended_double::pow2(e2);

							/* Check multiplication */
							const extended_double p = v1 * v2;
							const double p_q = p.fraction() / (fractions[f1] * fractions[f2]);
							BOOST_CHECK_EQUAL(std::ceil(log2(p_q)), log2(p_q));
							BOOST_CHECK_EQUAL(log2(p_q) + p.exponent(), e1 + e2);

							/* Check division */
							const extended_double q = v1 / v2;
							const double q_q = q.fraction() / (fractions[f1] / fractions[f2]);
							BOOST_CHECK_EQUAL(std::ceil(log2(q_q)), log2(q_q));
							BOOST_CHECK_EQUAL(log2(q_q) + q.exponent(), e1 - e2);

							/* Check addition */
							const extended_double s = v1 + v2;
							const double s1 = fractions[f1] * exp2(std::min(e1 - e2, 0));
							const double s2 = fractions[f2] * exp2(std::min(e2 - e1, 0));
							if ((s1 + s2) == 0.0)
								BOOST_CHECK_EQUAL(s, 0.0);
							else {
								const double s_q = s.fraction() / (s1 + s2);
								BOOST_CHECK_EQUAL(std::ceil(log2(s_q)), log2(s_q));
								BOOST_CHECK_EQUAL(log2(s_q) + s.exponent(),
												  std::max(e1, e2));
							}

							/* Check subtraction */
							const extended_double d = v1 - v2;
							const double d1 = fractions[f1] * exp2(std::min(e1 - e2, 0));
							const double d2 = fractions[f2] * exp2(std::min(e2 - e1, 0));
							if ((d1 - d2) == 0.0)
								BOOST_CHECK_EQUAL(d, 0.0);
							else {
								const double d_q = d.fraction() / (d1 - d2);
								BOOST_CHECK_EQUAL(std::ceil(log2(d_q)), log2(d_q));
								BOOST_CHECK_EQUAL(log2(d_q) + d.exponent(),
												  std::max(e1, e2));
							}
						}
					}
				}
			}
		}
	}
}

BOOST_AUTO_TEST_CASE(basic) {
	extended_double v0;
	BOOST_CHECK(v0 == extended_double(0.0));
	std::memset(&v0, 0, sizeof(v0));
	BOOST_CHECK(v0 == extended_double(0.0));
	const extended_double vs = extended_double(std::pow(2.0, -1000)) * extended_double(std::pow(2.0, -1000));
	BOOST_CHECK(extended_double_cast<double>(vs) == 0.0);
	const extended_double vh(0.5);
	BOOST_CHECK(extended_double_cast<double>(vh) == 0.5);
	const extended_double vx(0.13829348290834290348239);
	const extended_double v1(1.0);
	BOOST_CHECK(extended_double_cast<double>(v1) == 1.0);
	const extended_double v2(2.0);
	BOOST_CHECK(extended_double_cast<double>(v2) == 2.0);

	BOOST_CHECK((v0 < vs) && (vs < vx) && (vx < vh) && (vh < v1) && (v1 < v2));
	BOOST_CHECK((v0 <= vs) && (vs <= vx) && (vx <= vh) && (vh <= v1) && (v1 <= v2));
	BOOST_CHECK((v2 >= v1) && (v1 >= vh) && (vh >= vx) && (vx >= vs) && (vs >= v0));
	BOOST_CHECK((v2 > v1) && (v1 > vh) && (vh > vx) && (vx > vs) && (vs > v0));

	BOOST_CHECK(v0 * v0 == v0);
	BOOST_CHECK(v0 * vs == v0);
	BOOST_CHECK(v0 * vh == v0);
	BOOST_CHECK(v0 * vx == v0);
	BOOST_CHECK(v0 * v1 == v0);
	BOOST_CHECK(v0 * v2 == v0);

	BOOST_CHECK(v0 + v0 == v0);
	BOOST_CHECK(v0 + vs == vs);
	BOOST_CHECK(v0 + vh == vh);
	BOOST_CHECK(v0 + vx == vx);
	BOOST_CHECK(v0 + v1 == v1);
	BOOST_CHECK(v0 + v2 == v2);

	BOOST_CHECK(v1 / (v0 + v1 / vs) == vs);
	BOOST_CHECK(v1 / (v0 + v1 / vh) == vh);
	BOOST_CHECK(v1 / (v0 + v1 / vx) == vx);
	BOOST_CHECK(v1 / (v0 + v1 / v1) == v1);
	BOOST_CHECK(v1 / (v0 + v1 / v2) == v2);

	BOOST_CHECK((v0 * vs) / vs == v0);
	BOOST_CHECK((vs * vs) / vs == vs);
	BOOST_CHECK((vh * vs) / vs == vh);


	BOOST_CHECK((vx * vs) / vs == vx);
	BOOST_CHECK((v1 * vs) / vs == v1);
	BOOST_CHECK((v2 * vs) / vs == v2);

	BOOST_CHECK(vs + (-vs) == v0);
	BOOST_CHECK(vh + (-vh) == v0);
	BOOST_CHECK(vx + (-vx) == v0);
	BOOST_CHECK(v1 + (-v1) == v0);
	BOOST_CHECK(v2 + (-v2) == v0);

	BOOST_CHECK(vs - vs == v0);
	BOOST_CHECK(vh - vh == v0);
	BOOST_CHECK(vx - vx == v0);
	BOOST_CHECK(v1 - v1 == v0);
	BOOST_CHECK(v2 - v2 == v0);

	BOOST_CHECK(vs / vs == v1);
	BOOST_CHECK(vh / vh == v1);
	BOOST_CHECK(vx / vx == v1);
	BOOST_CHECK(v1 / v1 == v1);
	BOOST_CHECK(v2 / v2 == v1);

	const extended_double vexpm300 = std::exp(-300);
	BOOST_CHECK_CLOSE(log(vexpm300), -300, 1e-13);
	BOOST_CHECK_CLOSE(log(vexpm300*vexpm300), -600, 1e-13);
	BOOST_CHECK_CLOSE(log(vexpm300*vexpm300*vexpm300), -900, 1e-13);
	BOOST_CHECK_CLOSE(log(vexpm300*vexpm300*vexpm300*vexpm300), -1200, 1e-13);
	BOOST_CHECK_LT(vexpm300*vexpm300, vexpm300);
	BOOST_CHECK_LT(vexpm300*vexpm300*vexpm300, vexpm300*vexpm300);
	BOOST_CHECK_LT(vexpm300*vexpm300*vexpm300*vexpm300, vexpm300*vexpm300*vexpm300);
	BOOST_CHECK_EQUAL(vexpm300*vexpm300 + vexpm300, vexpm300);

	BOOST_CHECK_EQUAL(vs*vs*vs*vs, vs*vs*vs*vs);
	BOOST_CHECK_NE(vs*vs*vs*vs, vs*vs*vs);
	BOOST_CHECK_EQUAL(extended_double_cast<double>(vs*vs*vs*vs), 0.0);
}

BOOST_AUTO_TEST_SUITE_END()
