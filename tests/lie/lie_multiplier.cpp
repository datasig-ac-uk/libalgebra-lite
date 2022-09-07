//
// Created by user on 28/08/22.
//
#include <gtest/gtest.h>
#include <libalgebra_lite/hall_set.h>
#include <libalgebra_lite/lie.h>

using namespace lal;

struct LieMultiplierFixture : public ::testing::Test
{
    std::shared_ptr<const hall_basis> lie_basis;
    lie_multiplier multiplier;

    using key_type = typename hall_basis::key_type;
    using pref_type = const boost::container::small_vector_base<std::pair<key_type, int>>&;

    LieMultiplierFixture()
        : lie_basis(new hall_basis {5, 5}),
          multiplier(5)
    {}

};


TEST_F(LieMultiplierFixture, test_product_letters_in_order) {
    key_type k1(1, 0), k2(1, 1);

    const auto& result = multiplier(*lie_basis, k1, k2);

    ASSERT_EQ(result.size(), 1);
    EXPECT_EQ(result[0].first, key_type(2, 0));
    EXPECT_EQ(result[0].second, 1);
}

TEST_F(LieMultiplierFixture, test_product_letters_reverse_order) {
    key_type k1(1, 1), k2(1, 0);

    pref_type result = multiplier(*lie_basis, k1, k2);

    ASSERT_EQ(result.size(), 1);
    EXPECT_EQ(result[0].first, key_type(2, 0));
    EXPECT_EQ(result[0].second, -1);
}
