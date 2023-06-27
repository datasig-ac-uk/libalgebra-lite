//
// Created by sam on 07/09/22.
//

#ifndef LIBALGEBRA_LITE_TENSOR_FIXTURE_H
#define LIBALGEBRA_LITE_TENSOR_FIXTURE_H


#include <libalgebra_lite/free_tensor.h>
#include <libalgebra_lite/shuffle_tensor.h>

#include <libalgebra_lite/polynomial.h>

#include <gtest/gtest.h>

struct TensorFixture : public ::testing::Test
{
    lal::deg_t width = 5;
    lal::deg_t depth = 5;
    lal::basis_pointer<lal::tensor_basis> basis;

    using key_type = typename lal::tensor_basis::key_type;

    using scalar_type = typename lal::polynomial_ring::scalar_type;
    using rational_type = typename lal::polynomial_ring::rational_type;

    TensorFixture() : basis(new lal::tensor_basis(width, depth))
    {}

    ~TensorFixture() { delete const_cast<lal::tensor_basis*>(static_cast<const lal::tensor_basis*>(basis)); }

    key_type key(std::initializer_list<lal::let_t> arg) {
        typename key_type::index_type idx = 0;
        for (auto let : arg) {
            idx *= width;
            idx += let - 1;
        }
        return key_type(arg.size(), idx);
    }

    template <typename... Ints>
    key_type key(Ints... ints) {
        return key({lal::let_t(ints)...});
    }

    lal::monomial key_to_monomial(char letter, key_type arg) const {
        using letter_type = typename lal::monomial::letter_type;

        auto idx = arg.index();
        typename letter_type::integral_type tmp(0);
        decltype(idx) pow = 1;
        for (decltype(idx) i = 0; i < arg.degree(); ++i) {
            tmp += (1 + idx % width) * pow;
            pow *= 10;
            idx /= width;
        }
        return lal::monomial(letter_type(letter, tmp));
    }

    rational_type one() const { return rational_type(1); }

    template <typename TensorType>
    TensorType generic_tensor(char letter) const {
        const auto &powers = basis->powers();

        TensorType result(basis);

        for (lal::deg_t d = 0; d <= basis->depth(); ++d) {
            for (lal::dimn_t i = 0; i < powers[d]; ++i) {
                key_type key(d, i);
                result[key_type(d, i)] = scalar_type(key_to_monomial(letter, key), 1);
            }
        }
        return result;
    }

    template <typename Fn>
    void for_each_key(Fn op) const
    {
        const auto& powers = basis->powers();
        for (lal::deg_t d=0; d<=depth; ++d) {
            for (lal::dimn_t i=0; i<powers[d]; ++i) {
                key_type key(d, i);
                op(key);
            }
        }
    }

};




#endif //LIBALGEBRA_LITE_TENSOR_FIXTURE_H
