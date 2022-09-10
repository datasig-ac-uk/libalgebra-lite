//
// Created by user on 08/09/22.
//

#ifndef LIBALGEBRA_LITE_POLYNOMIAL_FIXTURE_H
#define LIBALGEBRA_LITE_POLYNOMIAL_FIXTURE_H

#include <libalgebra_lite/polynomial_basis.h>

#include <gtest/gtest.h>

struct PolynomialFixture : public ::testing::Test
{
    using deg_t = lal::deg_t;
    using dimn_t = lal::dimn_t;
    using let_t = lal::let_t;

    using monomial = lal::monomial;
    using letter_type = typename monomial::letter_type;


    std::shared_ptr<const lal::polynomial_basis> basis;

    PolynomialFixture() : basis(new lal::polynomial_basis)
    {}

};



#endif //LIBALGEBRA_LITE_POLYNOMIAL_FIXTURE_H
