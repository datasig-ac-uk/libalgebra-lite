//
// Created by user on 08/09/22.
//

#ifndef LIBALGEBRA_LITE_LIE_FIXTURE_H
#define LIBALGEBRA_LITE_LIE_FIXTURE_H

#include <libalgebra_lite/lie.h>
#include <gtest/gtest.h>



struct LieFixture : public ::testing::Test
{
    using deg_t = lal::deg_t;
    deg_t width = 5;
    deg_t depth = 5;

    using key_type = typename lal::hall_basis::key_type;

    std::shared_ptr<const lal::hall_basis> basis;

    LieFixture() : basis(new lal::hall_basis(width, depth))
    {}

    template <typename Let>
    key_type key(Let let) const
    {
        return basis->key_of_letter(lal::let_t(let));
    }

    key_type key(key_type arg) const
    {
        return arg;
    }

    template <typename Let>
    key_type key(std::initializer_list<Let> args) const
    {
        assert(args.size() == 2);
        return key(args[0], args[1]);
    }

    template <typename Left, typename Right>
    key_type key(Left left, Right right) const
    {
        return basis->find({key(left), key(right)}).it->second;
    }


};


#endif //LIBALGEBRA_LITE_LIE_FIXTURE_H
