//
// Created by sam on 07/09/22.
//


#include "tensor_fixture.h"



struct ShuffleMultiplierFixture : public TensorFixture
{
    lal::shuffle_tensor_multiplier multiplier;

    ShuffleMultiplierFixture() : TensorFixture(), multiplier(TensorFixture::width)
    {}

};


TEST_F(ShuffleMultiplierFixture, testMultiplierEmpty) {

    key_type lhs, rhs;

    auto result = multiplier(*basis, lhs, rhs);

    ASSERT_EQ(1, result.size());
    EXPECT_EQ(key_type(0, 0), result[0].first);
    EXPECT_EQ(1, result[0].second);
}


TEST_F(ShuffleMultiplierFixture, testMultiplierEmptyLetter) {
    key_type lhs1(1, 0), rhs1, lhs2, rhs2(1, 3);

    auto result1 = multiplier(*basis, lhs1, rhs1);
    auto result2 = multiplier(*basis, lhs2, rhs2);

    ASSERT_EQ(1, result1.size());
    EXPECT_EQ(lhs1, result1[0].first);
    EXPECT_EQ(1, result1[0].second);

    ASSERT_EQ(1, result2.size());
    EXPECT_EQ(rhs2, result2[0].first);
    EXPECT_EQ(1, result2[0].second);
}



TEST_F(ShuffleMultiplierFixture, testMultiplierLetterLetter) {
    key_type lhs(1, 0), rhs(1, 3);

    auto result = multiplier(*basis, lhs, rhs);

    ASSERT_EQ(2, result.size());
    EXPECT_EQ(key(4, 1), result[1].first);
    EXPECT_EQ(1, result[1].second);
    EXPECT_EQ(key(1, 4), result[0].first);
    EXPECT_EQ(1, result[0].second);

}


TEST_F(ShuffleMultiplierFixture, testMutliplierWordWord) {

    auto lhs = key(1, 2);
    auto rhs = key(3, 4);

    auto result = multiplier(*basis, lhs, rhs);

    ASSERT_EQ(6, result.size());

    EXPECT_EQ(key(1, 2, 3, 4), result[0].first);
    EXPECT_EQ(1, result[0].second);

    EXPECT_EQ(key(1, 3, 2, 4), result[1].first);
    EXPECT_EQ(1, result[1].second);

    EXPECT_EQ(key(1, 3, 4, 2), result[2].first);
    EXPECT_EQ(1, result[2].second);

    EXPECT_EQ(key(3, 1, 2, 4), result[3].first);
    EXPECT_EQ(1, result[3].second);

    EXPECT_EQ(key(3, 1, 4, 2), result[4].first);
    EXPECT_EQ(1, result[4].second);

    EXPECT_EQ(key(3, 4, 1, 2), result[5].first);
    EXPECT_EQ(1, result[5].second);
}
