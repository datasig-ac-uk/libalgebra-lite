//
// Created by user on 09/09/22.
//

#ifndef LIBALGEBRA_LITE_MAPS_FIXTURE_H
#define LIBALGEBRA_LITE_MAPS_FIXTURE_H

#include <libalgebra_lite/lie.h>
#include <libalgebra_lite/free_tensor.h>
#include <libalgebra_lite/maps.h>
#include <libalgebra_lite/coefficients.h>

#include <gtest/gtest.h>

class MapsFixture : public ::testing::Test
{
public:
    using deg_t = lal::deg_t;
    deg_t width = 5;
    deg_t depth = 5;

    std::shared_ptr<const lal::tensor_basis> tbasis;
    std::shared_ptr<const lal::hall_basis> lbasis;
    lal::maps maps;


    MapsFixture()
        : tbasis(new lal::tensor_basis(width, depth)),
          lbasis(new lal::hall_basis(width, depth)),
          maps(tbasis, lbasis)
    {}

    using tkey_type = typename lal::tensor_basis::key_type;
    using lkey_type = typename lal::hall_basis::key_type;

    tkey_type tkey(std::initializer_list<lal::let_t> lets) const {
        typename tkey_type::index_type idx(0);
        for (auto let : lets) {
            idx *= width;
            idx += let - 1;
        }
        return tkey_type(lets.size(), idx);
    }

    template <typename... Lets>
    tkey_type tkey(Lets... lets) const
    {
        return tkey({lal::let_t(lets)...});
    }

    lkey_type lkey(lkey_type key) const { return key; }

    template <typename Let>
    lkey_type lkey(Let let) const { return lbasis->key_of_letter(let); }

    template <typename Left, typename Right>
    lkey_type lkey(Left left, Right right) const
    {
        return lbasis->find({lkey(left), lkey(right)}).it->second;
    }

};


#endif //LIBALGEBRA_LITE_MAPS_FIXTURE_H
