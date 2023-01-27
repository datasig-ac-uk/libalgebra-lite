//
// Created by sam on 07/09/22.
//

#include "tensor_fixture.h"
#include <libalgebra_lite/polynomial.h>

#include <libalgebra_lite/unpacked_tensor_word.h>

//
template <typename TensorType>
struct FreeTensorMultiplicationFixture : public TensorFixture
{

public:

    std::shared_ptr<const lal::free_tensor_multiplication> multiplication;

    using tensor_type = TensorType;
    using scalar_type = typename lal::polynomial_ring::scalar_type;
    using rational_type = typename lal::polynomial_ring::rational_type;

    FreeTensorMultiplicationFixture()
        : multiplication(lal::multiplication_registry<lal::free_tensor_multiplication>::get(TensorFixture::width))
    {}

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

        auto deg = key.degree();
        if (deg == 0) {
            return { lal::monomial(letter_type(left, 0))*lal::monomial(letter_type(right, 0)) };
        }
        std::vector<lal::monomial> result;
        if (deg == 1) {
            result.reserve(2);
            result.emplace_back(lal::monomial(letter_type(left, 0))*lal::monomial(letter_type(right, key.index()+1)));
            result.emplace_back(lal::monomial(letter_type(left, key.index()+1))*lal::monomial(letter_type(right, 0)));
            return result;
        }

        result.reserve(key.degree()+1);


        const auto& powers = basis->powers();

        lal::unpacked_tensor_word tword(width, key);
        for (auto lh_deg = 0; lh_deg <= key.degree(); ++lh_deg) {

            std::map<letter_type, lal::deg_t> tmp;
            auto split = tword.split(lh_deg);
            tmp[letter_type(left, split.first.pack_with_base(10, 1))] = 1;
            tmp[letter_type(right, split.second.pack_with_base(10, 1))] = 1;
            result.emplace_back(tmp);
        }

        return result;
    }


};


TYPED_TEST_SUITE_P(FreeTensorMultiplicationFixture);

TYPED_TEST_P(FreeTensorMultiplicationFixture, testMultiplication) {

    auto lhs = this->generic_tensor('x');
    auto rhs = this->generic_tensor('y');

    const auto result = lhs*rhs;


    const auto& powers = this->basis->powers();

//    std::cout << result << '\n';

    const auto& first = result[TensorFixture::key_type(0, 0)];

    auto x0 = this->key_to_monomial('x', TensorFixture::key_type(0, 0));
    auto y0 = this->key_to_monomial('y', TensorFixture::key_type(0, 0));

    EXPECT_EQ(first[x0 * y0], this->one());

    for (auto d=1; d<=this->basis->depth(); ++d) {
        for (auto i=0; i<powers[d]; ++i) {
            TensorFixture::key_type key(d, i);
            const auto& val = result[key];
            this->basis->print_key(std::cout, key);
            auto components = this->deconstruct_key('x', 'y', key);

            ASSERT_EQ(val.size(), components.size());

            auto k = 0;
            for (const auto& item : val) {
                EXPECT_EQ(item.key(), components[k++]);
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
