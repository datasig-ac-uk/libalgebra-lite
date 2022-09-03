//
// Created by sam on 03/09/22.
//


#include "libalgebra_lite/free_tensor.h"


using namespace lal;


typename shuffle_tensor_multiplier::base_type::key_type
shuffle_tensor_multiplier::operator()(
        typename base_type::key_type lhs,
        typename base_type::key_type rhs) const
{
    if (lhs.degree() + rhs.degree() >= half_type::p_basis->depth()) {
        return typename half_type::product_type();
    }

    return half_type::shuffle(lhs, rhs);
}
