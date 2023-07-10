//
// Created by sam on 03/09/22.
//

#include "libalgebra_lite/shuffle_tensor.h"
#include "libalgebra_lite/free_tensor.h"


#include <mutex>
#include <unordered_map>



using namespace lal;


typename shuffle_tensor_multiplier::base_type::product_type
shuffle_tensor_multiplier::operator()(
        const tensor_basis& basis,
        typename base_type::key_type lhs,
        typename base_type::key_type rhs) const
{
    if (static_cast<deg_t>(lhs.degree() + rhs.degree()) > basis.depth()) {
        return {};
    }


    return half_type::shuffle(basis, lhs, rhs);
}



namespace lal {
template class multiplication_registry<right_half_shuffle_multiplication>;
template class multiplication_registry<left_half_shuffle_multiplication>;
template class multiplication_registry<shuffle_tensor_multiplication>;
}
