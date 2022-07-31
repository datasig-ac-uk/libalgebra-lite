#include "library.h"
#include <libalgebra_lite/implementation_types.h>
#include <libalgebra_lite/vectors/dense_vector.h>
#include <libalgebra_lite/basis/index_key.h>
#include <libalgebra_lite/basis/tensor_basis.h>
#include <libalgebra_lite/basis/hall_set.h>
#include <libalgebra_lite/coefficients.h>
#include <libalgebra_lite/basis/traits.h>
#include <iostream>

void hello()
{
    alg::dense_vector<alg::tensor_basis, alg::double_field> vect;
    std::cout << "Hello, World!" << std::endl;
}
