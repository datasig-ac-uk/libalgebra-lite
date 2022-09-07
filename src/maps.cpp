//
// Created by user on 05/09/22.
//

#include "libalgebra_lite/maps.h"

using namespace lal;

dtl::generic_commutator::tensor_type dtl::generic_commutator::operator()(
        const boost::container::small_vector_base<pair_type>& lhs,
        const boost::container::small_vector_base<pair_type>& rhs)
{
    std::map<key_type, int> tmp;

    for (const auto& lhs_val : lhs) {
        for (const auto& rhs_val : rhs) {
            for (const auto& inner : m_mul.multiply(m_basis, lhs_val.first, rhs_val.first)) {
                tmp[inner.first] += inner.second*lhs_val.second*rhs_val.second;
            }
            for (const auto& inner : m_mul.multiply(m_basis, rhs_val.first, lhs_val.first)) {
                tmp[inner.first] -= inner.second*rhs_val.second*lhs_val.second;
            }
        }
    }
    return {tmp.begin(), tmp.end()};
}

maps::generic_tensor maps::expand_letter(let_t letter)
{
    return {{tkey_type(1, letter), 1}};
}

maps::generic_lie maps::rbracketing_impl(maps::lkey_type lhs, maps::lkey_type rhs) const
{
}

typename maps::glie_ref maps::rbracketing(maps::tkey_type tkey) const
{
    static const boost::container::small_vector<lie_pair, 0> null;
    return null;
}
typename maps::gtensor_ref maps::expand(maps::lkey_type lkey) const
{
    return m_expand(lkey);
}
