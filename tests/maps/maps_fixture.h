//
// Created by user on 09/09/22.
//

#ifndef LIBALGEBRA_LITE_MAPS_FIXTURE_H
#define LIBALGEBRA_LITE_MAPS_FIXTURE_H

#include <libalgebra_lite/lie.h>
#include <libalgebra_lite/free_tensor.h>
#include <libalgebra_lite/maps.h>
#include <libalgebra_lite/coefficients.h>
#include <libalgebra_lite/polynomial.h>

#include <gtest/gtest.h>

class MapsFixture : public ::testing::Test
{
public:
    using deg_t = lal::deg_t;
    deg_t width = 5;
    deg_t depth = 5;

    lal::basis_pointer<lal::tensor_basis> tbasis;
    lal::basis_pointer<lal::hall_basis> lbasis;
    lal::maps maps;

    using poly_t = lal::polynomial_ring::scalar_type;
    using rational_type = lal::polynomial_ring::rational_type;
    using dtensor_t = lal::free_tensor<lal::polynomial_ring, lal::dense_vector, lal::dtl::standard_storage>;
    using dlie_t = lal::lie<lal::polynomial_ring, lal::dense_vector, lal::dtl::standard_storage>;

    MapsFixture()
        : tbasis(new lal::tensor_basis(width, depth)),
          lbasis(new lal::hall_basis(width, depth)),
          maps(tbasis, lbasis)
    {}


    ~MapsFixture() {
        delete const_cast<lal::tensor_basis*>(static_cast<const lal::tensor_basis*>(tbasis));
        delete const_cast<lal::hall_basis*>(static_cast<const lal::hall_basis*>(lbasis));
    }

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

private:

    template <typename S, typename T>
    static constexpr S pow(S base, T exponent) noexcept
    {
        if (exponent == 0) return S(1);
        if (exponent == 1) return base;
        auto half_pow = pow(base, exponent / 2);
        auto result = half_pow * half_pow;
        if (exponent & 1) result *= base;
        return result;
    }

    lal::dimn_t lkey_to_word_fragment(lkey_type key) const
    {
        if (key.degree() == 1) {
            return lal::dimn_t(lbasis->to_letter(key));
        }
        auto parents = lbasis->parents(key);
        if (key.degree() == 2) {
            return 900 + lkey_to_word_fragment(parents.first)*10 + lkey_to_word_fragment(parents.second);
        }

        // insert a 0 between groups of brackets
        lal::dimn_t shift = pow(lal::dimn_t(10), parents.second.degree()+1);
        auto top_shift = pow(lal::dimn_t(10), parents.first.degree())*shift;
        return 9*top_shift + lkey_to_word_fragment(parents.first)*shift + lkey_to_word_fragment(parents.second);
    }


    lal::monomial lkey_to_poly(char prefix, lkey_type key) const
    {
        return lal::monomial(typename lal::monomial::letter_type(prefix, lkey_to_word_fragment(key)), 1);
    }

    lal::monomial tkey_to_poly(char prefix, tkey_type key) const
    {
        lal::dimn_t word = 0;
        auto deg = static_cast<deg_t>(key.degree());
        auto index = key.index();
        for (deg_t i=0; i<deg; ++i) {
            word *= 10;
            auto old = index;
            index /= width;
            word += 1 + (old - index*width);
        }
        return lal::monomial(typename lal::monomial::letter_type(prefix, word), 1);
    }

public:
    dtensor_t generic_dtensor(char prefix, deg_t degree=3) const
    {
        const auto& sizes = tbasis->sizes();

        dtensor_t result(tbasis);
        for (deg_t d=0; d<=degree; ++d) {
            tkey_type k(d, 0);
            for (lal::dimn_t i=0; i<sizes[d]; ++i, ++k) {
                result[k] = poly_t(tkey_to_poly(prefix, k), rational_type(1));
            }
        }
        return result;
    }

    dlie_t generic_lie(char prefix, lal::deg_t degree=3) const
    {
        dlie_t result(lbasis);
        for (lal::deg_t d=1; d<=degree; ++d) {
            lkey_type k(d, 0);
            auto size = lbasis->size_of_degree(d);
            for (lal::dimn_t i=0; i<size; ++i, ++k) {
                result[k] = poly_t(lkey_to_poly(prefix, k), rational_type(1));
            }
        }
        return result;
    }

};


#endif //LIBALGEBRA_LITE_MAPS_FIXTURE_H
