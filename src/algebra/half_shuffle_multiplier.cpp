//
// Created by sam on 03/09/22.
//

#include "libalgebra_lite/free_tensor.h"





using namespace lal;

typename half_shuffle_tensor_multiplier::product_type
half_shuffle_tensor_multiplier::shuffle(key_type lhs, key_type rhs) const
{
    typename base_type::product_type result;

    const auto lhs_deg = lhs.degree();
    const auto rhs_deg = rhs.degree();

    if (lhs_deg == 0) {
        product_type result;
        result[rhs] = 1;
        return result;
    }
    if (rhs_deg == 0) {
        product_type result;
        result[lhs] = 1;
        return result;
    }

    return base_type::add(operator()(lhs, rhs), operator()(rhs, lhs));
}

typename half_shuffle_tensor_multiplier::product_type
lal::half_shuffle_tensor_multiplier::key_prod_impl(key_type lhs, key_type rhs) const
{
    product_type result;

    const auto lhs_deg = lhs.degree();
    const auto rhs_deg = rhs.degree();

    if (lhs_deg + rhs_deg <= p_basis->depth()) {

        if (lhs_deg == 0) {
            goto finish_half_shuffle;
        }
        if (rhs_deg == 0) {
            result[lhs] = 1;
            goto finish_half_shuffle;
        }

        const auto lparent = p_basis->lparent(lhs);
        const auto& right_part = shuffle(p_basis->rparent(lhs), rhs);

        result.reserve(right_part.size());
        for (const auto& item : right_part) {
            result.emplace_back(concat_product(lparent, item.first), item.second);
        }
    }

finish_half_shuffle:
    return result;
}


typename half_shuffle_tensor_multiplier::reference
half_shuffle_tensor_multiplier::operator()(key_type lhs, key_type rhs) const
{
    static const boost::container::small_vector<typename base_type::pair_type, 0> null;

    if (lhs.degree() + rhs.degree() >= p_basis->depth()) {
        return null;
    }

    parent_type parents{lhs, rhs};
    auto found = m_cache.find(parents);
    if (found != m_cache.end()) {
        return found->second;
    }

    return m_cache[parents] = key_prod_impl(lhs, rhs);
}
