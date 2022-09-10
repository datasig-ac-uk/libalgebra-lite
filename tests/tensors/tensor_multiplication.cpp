//
// Created by sam on 07/09/22.
//

#include "tensor_fixture.h"
#include <libalgebra_lite/polynomial.h>

//
template <typename TensorType>
struct FreeTensorMultiplicationFixture : public TensorFixture
{

public:

    std::shared_ptr<const lal::free_tensor_multiplication> multiplication;

    using tensor_type = TensorType;
    using scalar_type = typename lal::polynomial_ring::scalar_type;
    using rational_type = typename lal::polynomial_ring::rational_type;

    rational_type one() const { return rational_type(1); }

    lal::monomial key_to_monomial(char letter, key_type arg) const
    {
        using letter_type = typename lal::monomial::letter_type;

        auto idx = arg.index();
        typename letter_type::integral_type tmp(0);
        auto pow = 1;
        for (auto i=0; i<arg.degree(); ++i) {
            tmp += (1 + idx % width)*pow;
            pow *= 10;
            idx /= width;
        }
        return lal::monomial(letter_type(letter, tmp));
    }

    tensor_type generic_tensor(char letter) const
    {
        const auto& powers = basis->powers();

        tensor_type result(basis, multiplication);

        for (lal::deg_t d=0; d<=basis->depth(); ++d) {
            for (lal::dimn_t i=0; i<powers[d]; ++i) {
                key_type key(d, i);
                result[key_type(d, i)] = scalar_type(key_to_monomial(letter, key), 1);
            }
        }
        return result;
    }

    std::vector<lal::monomial> deconstruct_key(char left, char right, key_type key)
    {
        using letter_type = typename lal::monomial::letter_type;

        std::vector<lal::monomial> result;
        result.reserve(key.degree()+1);

        const auto& powers = basis->powers();

        for (auto lh_deg = 0; lh_deg <= key.degree(); ++lh_deg) {
            auto rh_deg = key.degree() - lh_deg;
            auto idx = key.index();

            std::map<letter_type, lal::deg_t> tmp;
            tmp[letter_type(left, idx/powers[lh_deg])] = 1;
            tmp[letter_type(right, idx%powers[lh_deg])] = 1;
            result.emplace_back(tmp);
        }

        return result;
    }


};


TYPED_TEST_SUITE_P(FreeTensorMultiplicationFixture);

TYPED_TEST_P(FreeTensorMultiplicationFixture, testMultiplication) {

    auto lhs = this->generic_tensor('x');
    auto rhs = this->generic_tensor('y');

    auto result = lhs*rhs;


    const auto& powers = this->basis->powers();

    const auto& first = result[TensorFixture::key_type(0, 0)];
    EXPECT_EQ(first[this->key_to_monomial('x', TensorFixture::key_type(0, 0))], this->one());
    EXPECT_EQ(first[this->key_to_monomial('y', TensorFixture::key_type(0, 0))], this->one());

    for (auto d=1; d<=this->basis->depth(); ++d) {
        for (auto i=0; i<powers[d]; ++i) {
            TensorFixture::key_type key(d, i);
            const auto& val = result[key];
            auto components = this->deconstruct_key('x', 'y', key);

            ASSERT_EQ(val.size(), components.size());

            auto k = 0;
            for (const auto& item : val) {
                EXPECT_EQ(item.key(), components[k]);
            }

        }

    }

}

REGISTER_TYPED_TEST_SUITE_P(FreeTensorMultiplicationFixture,
        testMultiplication);

//REGISTER_TYPED_TEST_SUITE_P(
//        FreeTensorMultiplication,
//        FreeTensorMultiplicationFixture,
//        )

using FreeTensorMultiplicationCases = ::testing::Types<
        lal::free_tensor<lal::polynomial_ring, lal::dense_vector, lal::dtl::standard_storage>
        >;
INSTANTIATE_TYPED_TEST_SUITE_P(FTM, FreeTensorMultiplicationFixture, FreeTensorMultiplicationCases);
