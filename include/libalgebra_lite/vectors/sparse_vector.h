//
// Created by user on 26/07/22.
//

#ifndef LIBALGEBRA_LITE_SPARSE_VECTOR_H
#define LIBALGEBRA_LITE_SPARSE_VECTOR_H

#include <libalgebra_lite/implementation_types.h>
#include <libalgebra_lite/basis/traits.h>
#include <libalgebra_lite/coefficients/coefficients.h>

#include <boost/container/flat_map.hpp>

namespace alg {

template <typename Basis, typename Coefficients>
class sparse_vector
{
    using basis_traits = basis_traits<Basis>;
    using coeff_traits = coefficient_trait<Coefficients>;

public:
    using basis_type = Basis;
    using key_type = typename basis_traits::key_type;

    using coefficient_ring = typename coeff_traits::coefficient_ring;
    using scalar_type = typename coeff_traits::scalar_type;
    using rational_type = typename coeff_traits::rational_type;

private:
    boost::container::flat_map<key_type, scalar_type> m_data;

public:



};


} // namespace alg

#endif //LIBALGEBRA_LITE_SPARSE_VECTOR_H
