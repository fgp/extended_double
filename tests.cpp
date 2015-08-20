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

	BOOST_CHECK_EQUAL(v_pinf, v_pinf * 1e-100);
	BOOST_CHECK_EQUAL(v_pinf, v_pinf * 1.0);
	BOOST_CHECK_EQUAL(v_pinf, v_pinf * 1e100);
	BOOST_CHECK_EQUAL(v_pinf, v_pinf * v_pinf);
	BOOST_CHECK_EQUAL(v_ninf, v_ninf * 1e-100);
	BOOST_CHECK_EQUAL(v_ninf, v_ninf * 1.0);
	BOOST_CHECK_EQUAL(v_ninf, v_ninf * 1e100);

	BOOST_CHECK_EQUAL(v_ninf, v_pinf * -1e-100);
	BOOST_CHECK_EQUAL(v_ninf, v_pinf * -1.0);
	BOOST_CHECK_EQUAL(v_ninf, v_pinf * -1e100);
	BOOST_CHECK_EQUAL(v_pinf, v_ninf * v_ninf);
	BOOST_CHECK_EQUAL(v_pinf, v_ninf * -1e-100);
	BOOST_CHECK_EQUAL(v_pinf, v_ninf * -1.0);
	BOOST_CHECK_EQUAL(v_pinf, v_ninf * -1e100);
	BOOST_CHECK_EQUAL(v_ninf, v_ninf * v_pinf);

	BOOST_CHECK_EQUAL(v_pinf, v_pinf / 1e-100);
	BOOST_CHECK_EQUAL(v_pinf, v_pinf / 1.0);
	BOOST_CHECK_EQUAL(v_pinf, v_pinf / 1e100);
	BOOST_CHECK_EQUAL(v_ninf, v_ninf / 1e-100);
	BOOST_CHECK_EQUAL(v_ninf, v_ninf / 1.0);
	BOOST_CHECK_EQUAL(v_ninf, v_ninf / 1e100);

	BOOST_CHECK_EQUAL(v_ninf, v_pinf / -1e-100);
	BOOST_CHECK_EQUAL(v_ninf, v_pinf / -1.0);
	BOOST_CHECK_EQUAL(v_ninf, v_pinf / -1e100);
	BOOST_CHECK_EQUAL(v_pinf, v_ninf / -1e-100);
	BOOST_CHECK_EQUAL(v_pinf, v_ninf / -1.0);
	BOOST_CHECK_EQUAL(v_pinf, v_ninf / -1e100);
}

BOOST_AUTO_TEST_CASE(normalization) {
    BOOST_CHECK_EQUAL(extended_double(std::ldexp(1.0,     0)).exponent(),   0);
    BOOST_CHECK_EQUAL(extended_double(std::ldexp(1.0,   255)).exponent(),   0);
	BOOST_CHECK_EQUAL(extended_double(std::ldexp(1.0,   256)).exponent(), 256);
	BOOST_CHECK_EQUAL(extended_double(std::ldexp(1.0,   511)).exponent(), 256);
	BOOST_CHECK_EQUAL(extended_double(std::ldexp(1.0,   512)).exponent(), 512);
	BOOST_CHECK_EQUAL(extended_double(std::ldexp(1.0,   767)).exponent(), 512);
	BOOST_CHECK_EQUAL(extended_double(std::ldexp(1.0,   768)).exponent(), 768);
	BOOST_CHECK_EQUAL(extended_double(std::ldexp(1.0,  1023)).exponent(), 768);

	BOOST_CHECK_EQUAL(extended_double(std::ldexp(1.0,    -1)).exponent(),  -256);
	BOOST_CHECK_EQUAL(extended_double(std::ldexp(1.0,  -256)).exponent(),  -256);
	BOOST_CHECK_EQUAL(extended_double(std::ldexp(1.0,  -257)).exponent(),  -512);
    BOOST_CHECK_EQUAL(extended_double(std::ldexp(1.0,  -512)).exponent(),  -512);
    BOOST_CHECK_EQUAL(extended_double(std::ldexp(1.0,  -513)).exponent(),  -768);
    BOOST_CHECK_EQUAL(extended_double(std::ldexp(1.0,  -768)).exponent(),  -768);
    BOOST_CHECK_EQUAL(extended_double(std::ldexp(1.0,  -769)).exponent(), -1024);
	BOOST_CHECK_EQUAL(extended_double(std::ldexp(1.0, -1022)).exponent(), -1024);

    BOOST_CHECK_EQUAL(extended_double_cast<double>(extended_double(std::ldexp(1.0, 0))),
                      1.0);
    BOOST_CHECK_EQUAL(extended_double_cast<double>(extended_double(std::ldexp(1.0, 255))),
                      std::ldexp(1.0, 255));
    BOOST_CHECK_EQUAL(extended_double_cast<double>(extended_double(std::ldexp(1.0, 256))),
                      std::ldexp(1.0, 256));
    BOOST_CHECK_EQUAL(extended_double_cast<double>(extended_double(std::ldexp(1.0, 511))),
                      std::ldexp(1.0, 511));
    BOOST_CHECK_EQUAL(extended_double_cast<double>(extended_double(std::ldexp(1.0, 512))),
                      std::ldexp(1.0, 512));
    BOOST_CHECK_EQUAL(extended_double_cast<double>(extended_double(std::ldexp(1.0, 767))),
                      std::ldexp(1.0, 767));
    BOOST_CHECK_EQUAL(extended_double_cast<double>(extended_double(std::ldexp(1.0, 768))),
                      std::ldexp(1.0, 768));
    BOOST_CHECK_EQUAL(extended_double_cast<double>(extended_double(std::ldexp(1.0, 1023))),
                      std::ldexp(1.0, 1023));
    BOOST_CHECK_EQUAL(extended_double_cast<double>(extended_double(std::ldexp(1.0,  1024))),
                      std::numeric_limits<double>::infinity());
    
    BOOST_CHECK_EQUAL(extended_double_cast<double>(extended_double(std::ldexp(-1.0, 0))),
                      -1.0);
    BOOST_CHECK_EQUAL(extended_double_cast<double>(extended_double(std::ldexp(-1.0, 255))),
                      std::ldexp(-1.0, 255));
    BOOST_CHECK_EQUAL(extended_double_cast<double>(extended_double(std::ldexp(-1.0, 256))),
                      std::ldexp(-1.0, 256));
    BOOST_CHECK_EQUAL(extended_double_cast<double>(extended_double(std::ldexp(-1.0, 511))),
                      std::ldexp(-1.0, 511));
    BOOST_CHECK_EQUAL(extended_double_cast<double>(extended_double(std::ldexp(-1.0, 512))),
                      std::ldexp(-1.0, 512));
    BOOST_CHECK_EQUAL(extended_double_cast<double>(extended_double(std::ldexp(-1.0, 767))),
                      std::ldexp(-1.0, 767));
    BOOST_CHECK_EQUAL(extended_double_cast<double>(extended_double(std::ldexp(-1.0, 768))),
                      std::ldexp(-1.0, 768));
    BOOST_CHECK_EQUAL(extended_double_cast<double>(extended_double(std::ldexp(-1.0, 1023))),
                      std::ldexp(-1.0, 1023));
    BOOST_CHECK_EQUAL(extended_double_cast<double>(extended_double(std::ldexp(-1.0,  1024))),
                      -std::numeric_limits<double>::infinity());
    
    BOOST_CHECK_EQUAL(extended_double_cast<double>(extended_double(std::ldexp(1.0, -255))),
                      std::ldexp(1.0, -255));
    BOOST_CHECK_EQUAL(extended_double_cast<double>(extended_double(std::ldexp(1.0, -256))),
                      std::ldexp(1.0, -256));
    BOOST_CHECK_EQUAL(extended_double_cast<double>(extended_double(std::ldexp(1.0, -511))),
                      std::ldexp(1.0, -511));
    BOOST_CHECK_EQUAL(extended_double_cast<double>(extended_double(std::ldexp(1.0, -512))),
                      std::ldexp(1.0, -512));
    BOOST_CHECK_EQUAL(extended_double_cast<double>(extended_double(std::ldexp(1.0, -767))),
                      std::ldexp(1.0, -767));
    BOOST_CHECK_EQUAL(extended_double_cast<double>(extended_double(std::ldexp(1.0, -768))),
                      std::ldexp(1.0, -768));
    BOOST_CHECK_EQUAL(extended_double_cast<double>(extended_double(std::ldexp(1.0, -1022))),
                      std::ldexp(1.0, -1022));
    BOOST_CHECK_EQUAL(extended_double_cast<double>(extended_double(std::ldexp(1.0, -1023))),
                      0.0);
    
    BOOST_CHECK_EQUAL(extended_double_cast<double>(extended_double(std::ldexp(-1.0, -255))),
                      std::ldexp(-1.0, -255));
    BOOST_CHECK_EQUAL(extended_double_cast<double>(extended_double(std::ldexp(-1.0, -256))),
                      std::ldexp(-1.0, -256));
    BOOST_CHECK_EQUAL(extended_double_cast<double>(extended_double(std::ldexp(-1.0, -511))),
                      std::ldexp(-1.0, -511));
    BOOST_CHECK_EQUAL(extended_double_cast<double>(extended_double(std::ldexp(-1.0, -512))),
                      std::ldexp(-1.0, -512));
    BOOST_CHECK_EQUAL(extended_double_cast<double>(extended_double(std::ldexp(-1.0, -767))),
                      std::ldexp(-1.0, -767));
    BOOST_CHECK_EQUAL(extended_double_cast<double>(extended_double(std::ldexp(-1.0, -768))),
                      std::ldexp(-1.0, -768));
    BOOST_CHECK_EQUAL(extended_double_cast<double>(extended_double(std::ldexp(-1.0, -1022))),
                      std::ldexp(-1.0, -1022));
    BOOST_CHECK_EQUAL(extended_double_cast<double>(extended_double(std::ldexp(-1.0, -1023))),
                      0.0);
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
