//
// Created by user on 25/07/22.
//

#ifndef LIBALGEBRA_LITE_TRAITS_H
#define LIBALGEBRA_LITE_TRAITS_H

#include <libalgebra_lite/implementation_types.h>

#include <functional>
#include <utility>
#include <memory>

namespace alg {

struct with_degree_tag {
    deg_t degree;
};
struct without_degree_tag { };

namespace dtl {

template <typename Key, typename KeyOrder>
struct key_value_ordering
{
    template <typename S>
    using pair = std::pair<Key, S>;

    KeyOrder order;

    template <typename S>
    bool operator()(const pair<S>& lhs, const pair<S>& rhs) const noexcept
    {
        return order(lhs.first, rhs.first);
    }

};

} // namespace dtl

template <typename Basis>
struct basis_trait {
    using key_type = typename Basis::key_type;
    using degree_tag = without_degree_tag;

    static dimn_t max_dimension(const Basis& basis) noexcept { return basis.size(-1); };

    static dimn_t key_to_index(const Basis& basis, const key_type& k) noexcept { return basis.key_to_index(k); }
    static key_type index_to_key(const Basis& basis, dimn_t idx) noexcept { return basis.index_to_key(idx); }
    static deg_t degree(const Basis& basis, const key_type& key) noexcept
    { return basis.degree(key); }
    static deg_t max_degree(const Basis& basis) noexcept
    { return basis.max_degree(); }

    static dimn_t start_of_degree(const Basis& basis, deg_t deg) noexcept
    { return basis.start_of_degree(deg); }
    static dimn_t size(const Basis& basis, deg_t deg) noexcept
    { return basis.size(static_cast<int>(deg)); }
    static std::pair<dimn_t, deg_t> get_next_dimension(const Basis& basis, dimn_t dim, deg_t hint=0)
    {
        auto key = index_to_key(basis, dim);
        auto deg = degree(basis, key);
        return {size(basis, deg), deg};
    }

    using key_ordering = std::less<key_type>;
    using kv_ordering = dtl::key_value_ordering<key_type, key_ordering>;

    template <typename... Args>
    static std::shared_ptr<const Basis> basis_factory(Args... args)
    {
        return { nullptr };
    }


};


template <typename Basis1, typename Basis2>
struct is_basis_compatible : std::false_type
{};



} // namespace alg

#endif //LIBALGEBRA_LITE_TRAITS_H
