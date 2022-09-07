//
// Created by sam on 07/09/22.
//

#include "tensor_fixture.h"




struct FreeTensorMultiplierFixture : public TensorFixture
{
    lal::free_tensor_multiplier multiplier;


    FreeTensorMultiplierFixture() : TensorFixture(), multiplier(TensorFixture::width)
    {}

};


TEST_F(FreeTensorMultiplierFixture, testMutliplyEmptyKey) {

    key_type lhs, rhs, expected;

    auto result = multiplier(*basis, lhs, rhs);

    ASSERT_EQ(1, result.size());
    EXPECT_EQ(expected, result[0].first);
    EXPECT_EQ(1, result[0].second);
}

TEST_F(FreeTensorMultiplierFixture, testMultiplierEmptyLetter) {

    key_type lhs1(1, 0), rhs1, lhs2, rhs2(1, 3);

    auto expected1(lhs1);
    auto expected2(rhs2);

    auto result1 = multiplier(*basis, lhs1, rhs1);
    auto result2 = multiplier(*basis, lhs2, rhs2);

    ASSERT_EQ(1, result1.size());
    EXPECT_EQ(expected1, result1[0].first);
    EXPECT_EQ(1, result1[0].second);

    ASSERT_EQ(1, result2.size());
    EXPECT_EQ(expected2, result2[0].first);
    EXPECT_EQ(1, result2[0].second);
}


TEST_F(FreeTensorMultiplierFixture, testMultiplierTwoLetters) {
    key_type lhs(1, 0), rhs(1, 2);

    key_type expected(2, 2);

    auto result = multiplier(*basis, lhs, rhs);

    ASSERT_EQ(1, result.size());
    EXPECT_EQ(expected, result[0].first);
    EXPECT_EQ(1, result[0].second);
}


TEST_F(FreeTensorMultiplierFixture, testMultiplierOutOfBounds) {
    key_type lhs(5, 0), rhs(5, 0);

    auto result = multiplier(*basis, lhs, rhs);

    ASSERT_EQ(0, result.size());
}
