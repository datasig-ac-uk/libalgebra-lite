//
// Created by user on 27/08/22.
//


#include "libalgebra_lite/lie.h"

#include <map>

namespace lal {

lie_multiplier::product_type
lie_multiplier::key_prod_impl(key_type lhs, key_type rhs) const
{
    if (lhs>rhs) {
        return generic_minus(operator()(rhs, lhs));
    }

    product_type result;
    if (p_basis->degree(lhs) + p_basis->degree(rhs) > p_basis->depth()) {
        return result;
    }

    auto found = p_basis->find(parent_type(lhs, rhs));
    if (found.found) {
        result.emplace_back(found.it->second, 1);
    } else {
        auto lparent = p_basis->lparent(rhs);
        auto rparent = p_basis->rparent(rhs);

        auto result_left = generic_by_key(operator()(lhs, lparent), rparent);
        auto result_right = generic_by_key(operator()(lhs, rparent), lparent);

        return generic_sub(result_left, result_right);
    }

    return result;
}

lie_multiplier::product_type lie_multiplier::generic_minus(
        product_ref_type arg)
{
    product_type result;
    result.reserve(arg.size());

    for (const auto& item : arg) {
        result.emplace_back(item.first, -item.second);
    }

    return result;
}
lie_multiplier::product_type lie_multiplier::generic_sub(
        const boost::container::small_vector_base<std::pair<key_type, int>>& lhs,
        const boost::container::small_vector_base<std::pair<key_type, int>>& rhs)
{
    std::map<key_type, int> tmp;
    tmp.insert(lhs.begin(), lhs.end());

    for (const auto& item : rhs) {
        tmp[item.first] -= item.second;
    }

    return {tmp.begin(), tmp.end()};
}

lie_multiplier::product_type lie_multiplier::generic_by_key(
        const boost::container::small_vector_base<std::pair<key_type, int>>& lhs, lie_multiplier::key_type rhs) const
{
    std::map<key_type, int> tmp;

    for (const auto& outer : lhs) {
        for (const auto& inner : operator()(outer.first, rhs)) {
            tmp[inner.first] += outer.second * inner.second;
        }
    }

    return {tmp.begin(), tmp.end()};
}

lie_multiplier::product_ref_type lie_multiplier::operator()(
        lie_multiplier::key_type lhs, lie_multiplier::key_type rhs) const
{
    static const product_type null;
    if (lhs == rhs) {
        return null;
    }

    std::lock_guard<std::recursive_mutex> access(m_lock);

    parent_type parents {lhs, rhs};
    auto& found = m_cache[parents];
    if (!found.empty()) {
        return found;
    }

    found = key_prod_impl(lhs, rhs);
    return found;
}


} // namespace alg
