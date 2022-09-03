//
// Created by user on 15/08/22.
//

#include <gtest/gtest.h>
#define LAL_TESTING
#include "libalgebra_lite/packed_integer.h"
#undef LAL_TESTING

#include <cstdint>

namespace lal {
namespace dtl {

struct packed_integer_access
{

    template <typename I, typename P>
    static I raw_data(const packed_integer<I, P>& arg) {
        return arg.m_data;
    }

    template <typename I, typename P>
    static I integral_mask(const packed_integer<I, P>&) {
        return packed_integer<I, P>::integral_mask;
    }

    template <typename I, typename P>
    static I packed_mask(const packed_integer<I, P>&) {
        return packed_integer<I, P>::packed_mask;
    }

    template <typename I, typename P>
    static int remaining_bits(const packed_integer<I, P>&)
    {
        return packed_integer<I, P>::remaining_bits;
    }


};

}
}





struct PackedIntegerFixture : public ::testing::Test
{
    using integral_t = std::uint64_t;
    using packed_t = char;

    using pint_t = lal::dtl::packed_integer<std::uint64_t, char>;

};

TEST_F(PackedIntegerFixture, test_remaining_bits) {
    ASSERT_EQ(64-8, lal::dtl::packed_integer_access::remaining_bits(pint_t(0, 0)));
}


TEST_F(PackedIntegerFixture, test_create_packed) {
    packed_t item = 'x';
    pint_t packed(item);

    EXPECT_EQ(static_cast<packed_t>(packed), item);
    EXPECT_EQ(static_cast<integral_t>(packed), 0);
}

TEST_F(PackedIntegerFixture, test_create_integeral) {
    integral_t item = 1234567;
    pint_t packed(item);

    EXPECT_EQ(static_cast<packed_t>(packed), packed_t(0));
    EXPECT_EQ(static_cast<integral_t>(packed), item);
}

TEST_F(PackedIntegerFixture, test_create_both) {
    packed_t item1 = 'x';
    integral_t item2 = 1234567;
    pint_t packed(item1, item2);

    EXPECT_EQ(static_cast<packed_t>(packed), item1);
    EXPECT_EQ(static_cast<integral_t>(packed), item2);
}
