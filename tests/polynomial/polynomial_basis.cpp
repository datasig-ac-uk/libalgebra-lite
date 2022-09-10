//
// Created by user on 08/09/22.
//

#include "polynomial_fixture.h"

#include <libalgebra_lite/coefficients.h>


TEST_F(PolynomialFixture, testDegreeEmptyKey) {
    monomial k;

    EXPECT_EQ(0, basis->degree(k));
}

TEST_F(PolynomialFixture, testDegreeSingleVariable) {

    monomial k {letter_type('x', 1)};

    EXPECT_EQ(1, basis->degree(k));
}

TEST_F(PolynomialFixture, testDegreeSingleVarPower2) {
    monomial k(letter_type('x', 1), 2);

    EXPECT_EQ(2, basis->degree(k));
}

TEST_F(PolynomialFixture, testDegreeTwoVars) {
    std::map<letter_type, deg_t> tmp;
    tmp[letter_type('x', 1)] = 1;
    tmp[letter_type('y', 1)] = 1;

    monomial k(tmp);

    EXPECT_EQ(2, basis->degree(k));
}

TEST_F(PolynomialFixture, testDegreeTwoVarPower2) {
    std::map<letter_type, deg_t> tmp;
    tmp[letter_type('x', 1)] = 2;
    tmp[letter_type('x', 2)] = 1;

    EXPECT_EQ(3, basis->degree(monomial(tmp)));
}

TEST_F(PolynomialFixture, testEvalSingleLetter) {
    letter_type letter('x', 1);

    std::map<letter_type, typename lal::rational_field::scalar_type> map;
    map[letter] = 2;

    monomial k(letter, 2);

//    EXPECT_EQ(typename lal::rational_field::scalar_type(4),
//            k.eval<lal::rational_field>(map));

}
