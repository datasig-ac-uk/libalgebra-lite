//
// Created by sam on 07/09/22.
//


#include "tensor_fixture.h"


struct HalfShuffleMultiplierFixture : public TensorFixture
{

    lal::half_shuffle_tensor_multiplier multiplier;

    HalfShuffleMultiplierFixture() : TensorFixture(), multiplier(TensorFixture::width)
    {}

};



TEST_F(HalfShuffleMultiplierFixture, testHalfShuffleEmptyEmpty) {
    auto lhs = key();
    auto rhs = key();

    auto& result = multiplier(*basis, lhs, rhs);

    ASSERT_EQ(0, result.size());
//    EXPECT_EQ(key(), result[0].first);
//    EXPECT_EQ(1, result[0].second);
}

TEST_F(HalfShuffleMultiplierFixture, testHalfShuffleEmtpyLetter) {
    auto lhs = key();
    auto rhs = key(1);

    auto& result = multiplier(*basis, lhs, rhs);

    ASSERT_EQ(0, result.size());
//    EXPECT_EQ(key(1), result[0].first);
//    EXPECT_EQ(1, result[0].second);
}
