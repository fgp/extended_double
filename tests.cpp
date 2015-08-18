/*
 *  test/benchmark.cpp
 *  hoppe
 *
 *  Created by Florian Pflug on 03.12.2014.
 *  Copyright (c) 2014 Florian Pflug. All rights reserved.
 */

#include <boost/test/unit_test.hpp>

#include "extended_double.h"

BOOST_AUTO_TEST_SUITE(extended_double_tests)

BOOST_AUTO_TEST_CASE(infinities) {
	const extended_double v_pinf(std::numeric_limits<double>::infinity());
	const extended_double v_ninf(-std::numeric_limits<double>::infinity());

	BOOST_CHECK_EQUAL(extended_double_cast<double>(v_pinf),
					  std::numeric_limits<double>::infinity());
	BOOST_CHECK_EQUAL(extended_double_cast<double>(v_ninf),
					  -std::numeric_limits<double>::infinity());
	BOOST_CHECK(!isfinite(v_pinf));
	BOOST_CHECK(!isfinite(v_ninf));

	BOOST_CHECK_EQUAL(v_pinf, v_pinf);
	BOOST_CHECK_EQUAL(v_ninf, v_ninf);
	BOOST_CHECK_EQUAL(v_pinf, -v_ninf);
	BOOST_CHECK_EQUAL(-v_pinf, v_ninf);

	BOOST_CHECK(v_ninf < v_pinf);
	BOOST_CHECK(v_ninf <= v_pinf);
	BOOST_CHECK(v_ninf != v_pinf);
	BOOST_CHECK(v_pinf >= v_ninf);
	BOOST_CHECK(v_pinf > v_ninf);

	BOOST_CHECK_EQUAL(v_pinf, v_pinf + 1.0);
	BOOST_CHECK_EQUAL(v_pinf, v_pinf - 1.0);
	BOOST_CHECK_EQUAL(v_pinf, v_pinf + v_pinf);
	BOOST_CHECK_EQUAL(v_pinf, v_pinf - v_ninf);

	BOOST_CHECK_EQUAL(v_pinf, v_pinf * 0.5);
	BOOST_CHECK_EQUAL(v_pinf, v_pinf * 1.0);
	BOOST_CHECK_EQUAL(v_pinf, v_pinf * 2.0);
	BOOST_CHECK_EQUAL(v_pinf, v_pinf * v_pinf);
	BOOST_CHECK_EQUAL(v_ninf, v_ninf * 0.5);
	BOOST_CHECK_EQUAL(v_ninf, v_ninf * 1.0);
	BOOST_CHECK_EQUAL(v_ninf, v_ninf * 2.0);

	BOOST_CHECK_EQUAL(v_ninf, v_pinf * -0.5);
	BOOST_CHECK_EQUAL(v_ninf, v_pinf * -1.0);
	BOOST_CHECK_EQUAL(v_ninf, v_pinf * -2.0);
	BOOST_CHECK_EQUAL(v_pinf, v_ninf * v_ninf);
	BOOST_CHECK_EQUAL(v_pinf, v_ninf * -0.5);
	BOOST_CHECK_EQUAL(v_pinf, v_ninf * -1.0);
	BOOST_CHECK_EQUAL(v_pinf, v_ninf * -2.0);
	BOOST_CHECK_EQUAL(v_ninf, v_ninf * v_pinf);

	BOOST_CHECK_EQUAL(v_pinf, v_pinf / 0.5);
	BOOST_CHECK_EQUAL(v_pinf, v_pinf / 1.0);
	BOOST_CHECK_EQUAL(v_pinf, v_pinf / 2.0);
	BOOST_CHECK_EQUAL(v_ninf, v_ninf / 0.5);
	BOOST_CHECK_EQUAL(v_ninf, v_ninf / 1.0);
	BOOST_CHECK_EQUAL(v_ninf, v_ninf / 2.0);

	BOOST_CHECK_EQUAL(v_ninf, v_pinf / -0.5);
	BOOST_CHECK_EQUAL(v_ninf, v_pinf / -1.0);
	BOOST_CHECK_EQUAL(v_ninf, v_pinf / -2.0);
	BOOST_CHECK_EQUAL(v_pinf, v_ninf / -0.5);
	BOOST_CHECK_EQUAL(v_pinf, v_ninf / -1.0);
	BOOST_CHECK_EQUAL(v_pinf, v_ninf / -2.0);
}

BOOST_AUTO_TEST_CASE(basic)
{
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
	BOOST_CHECK(log(vexpm300) == -300);
	BOOST_CHECK(log(vexpm300*vexpm300) == -600);
	BOOST_CHECK(log(vexpm300*vexpm300*vexpm300) == -900);
	BOOST_CHECK(log(vexpm300*vexpm300*vexpm300*vexpm300) == -1200);
	BOOST_CHECK(vexpm300*vexpm300 < vexpm300);
	BOOST_CHECK(vexpm300*vexpm300*vexpm300 < vexpm300*vexpm300);
	BOOST_CHECK(vexpm300*vexpm300*vexpm300*vexpm300 < vexpm300*vexpm300*vexpm300);
	BOOST_CHECK(vexpm300*vexpm300 + vexpm300 == vexpm300);

	BOOST_CHECK(vs*vs*vs*vs == vs*vs*vs*vs);
	BOOST_CHECK(vs*vs*vs*vs != vs*vs*vs);
	BOOST_CHECK(extended_double_cast<double>(vs*vs*vs*vs) == 0.0);
}

BOOST_AUTO_TEST_SUITE_END()
