/*
 *  tests.cpp
 *  extended_double
 *
 *  Created by Florian Pflug on 17.08.2015.
 *  Copyright (c) 2015 Florian Pflug. All rights reserved.
 */

#include "extended_double.h"

extended_double multiply(extended_double a, extended_double b);

extended_double add(extended_double a, extended_double b);

extended_double multiply(extended_double a, extended_double b) {
    return a * b;
}

extended_double add(extended_double a, extended_double b) {
    return a + b;
}
