//
// Created by user on 01/09/22.
//


#include "libalgebra_lite/coefficients.h"
#include "libalgebra_lite/polynomial.h"


namespace lal {

template class coefficient_ring<polynomial<float_field>, float>;
template class coefficient_ring<polynomial<double_field>, double>;
template class coefficient_ring<polynomial<rational_field>, typename rational_field::scalar_type>;


}
