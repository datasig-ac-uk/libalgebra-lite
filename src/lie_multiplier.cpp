//
// Created by user on 27/08/22.
//


#include "libalgebra_lite/lie.h"


namespace lal {

template class base_multiplier<lie_multiplier, hall_basis, 2>;

lie_multiplier::product_type
lie_multiplier::key_prod_impl(key_type lhs, key_type rhs) const
{
    if (lhs>rhs) {
        return lie_multiplier::uminus(operator()(rhs, lhs));
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

        auto result_left = mul(operator()(lhs, lparent), rparent);
        auto result_right = mul(operator()(lhs, rparent), lparent);

        return sub(result_left, result_right);
    }

    return result;
}


lie_multiplier::reference lie_multiplier::operator()(
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
