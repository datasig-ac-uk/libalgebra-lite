//
// Created by user on 25/07/22.
//

#ifndef LIBALGEBRA_LITE_TRAITS_H
#define LIBALGEBRA_LITE_TRAITS_H

#include <libalgebra_lite/implementation_types.h>
#include <utility>

namespace alg {

template <typename Basis>
struct basis_trait {
    using key_type = typename Basis::key_type;

    static dimn_t max_dimension(const Basis& basis) noexcept { return basis.size(-1); };

    static dimn_t key_to_index(const Basis& basis, const key_type& k) noexcept { return basis.key_to_index(k); }
    static key_type index_to_key(const Basis& basis, dimn_t idx) noexcept { return basis.index_to_key(idx); }
    static deg_t degree(const Basis& basis, const key_type& key) noexcept
    { return basis.degree(key); }

    static dimn_t start_of_degree(const Basis& basis, deg_t deg) noexcept
    { return basis.start_of_degree(deg); }
    static dimn_t size(const Basis& basis, deg_t deg) noexcept
    { return basis.size(static_cast<int>(deg)); }
    std::pair<dimn_t, deg_t> get_dimension_degree(const Basis& basis, dimn_t dim)
    {
        auto key = index_to_key(basis, dim);
        auto deg = degree(basis, key);
        return {size(basis, deg), deg};
    }

};


} // namespace alg

#endif //LIBALGEBRA_LITE_TRAITS_H
