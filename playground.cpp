//
// Created by user on 07/08/22.
//

#include <libalgebra_lite/algebra.h>
#include <libalgebra_lite/coefficients.h>
#include <libalgebra_lite/vectors/dense_vector.h>
#include <libalgebra_lite/vector.h>

#include <iostream>

template <alg::dimn_t Dimension>
struct integer_basis
{
    using key_type = unsigned;
    using idx_type = alg::dimn_t;

    static constexpr idx_type key_to_index(key_type key) noexcept { return key; }
    static constexpr key_type index_to_key(idx_type idx) noexcept { return idx; }

    static constexpr alg::dimn_t size(int) { return Dimension; }
    static constexpr alg::deg_t degree(key_type) noexcept { return 0; }

};


using basis_type = integer_basis<5>;



namespace alg {

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
using multiplier = alg::basic_multiplier<basis_type>;
using vector_type = alg::vector<basis_type, alg::float_field, alg::dense_vector>;


int main()
{
    using traits = alg::multiplication_traits<alg::base_multiplication<multiplier>>;
    alg::base_multiplication<multiplier> mult{};
    auto basis = new basis_type;

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

    delete basis;
    return 0;
}
