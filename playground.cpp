//
// Created by user on 07/08/22.
//

#include "libalgebra_lite/algebra.h"
#include "libalgebra_lite/coefficients.h"
#include "libalgebra_lite/dense_vector.h"
#include "libalgebra_lite/vector.h"

#include <iostream>
#include <memory>

template <lal::dimn_t Dimension>
struct integer_basis
{
    using key_type = unsigned;
    using idx_type = lal::dimn_t;

    static constexpr idx_type key_to_index(key_type key) noexcept { return key; }
    static constexpr key_type index_to_key(idx_type idx) noexcept { return idx; }

    static constexpr lal::dimn_t size(int) { return Dimension; }
    static constexpr lal::deg_t degree(key_type) noexcept { return 0; }

};


using basis_type = integer_basis<5>;



namespace lal {

template <>
struct basic_multiplier<basis_type>
{

    using key_type = typename basis_type::key_type;

    std::vector<std::pair<key_type, int>> operator()(const key_type& lhs, const key_type& rhs) const
    {
        if (lhs == rhs) { return {{ lhs, 1 }}; }
        return {};
    }

};

} // namespace alg

template <typename B, typename C>
using dense_vector = lal::dense_vector<B, C>;

using multiplier = lal::basic_multiplier<basis_type>;
using vector_type = lal::vector<basis_type, lal::float_field, dense_vector>;


int main()
{
    using traits = lal::multiplication_traits<lal::base_multiplication<multiplier>>;
    auto basis = std::make_shared<const basis_type>();
    auto mult = std::make_shared<lal::base_multiplication<multiplier>>(basis);

    vector_type lhs (basis, { 1., 2., 3., 4., 5. });
    std::cout << "lhs\n";
    for (auto item : lhs) {
        std::cout << item.key() << ' ' << item.value() << '\n';
    }
    vector_type rhs (basis, { 2., 2., 2., 2., 2.});
    std::cout << "rhs\n";
    for (auto item : rhs) {
        std::cout << item.key() << ' ' << item.value() << '\n';
    }
    vector_type result (basis);

    traits::multiply_and_add(mult, result, lhs, rhs, [](auto t) { return t; });

    std::cout << "result\n";
    for (auto item : result) {
        std::cout << item.key() << ' ' << item.value() << '\n';
    }

    return 0;
}
