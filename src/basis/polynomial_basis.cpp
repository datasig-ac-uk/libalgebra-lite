//
// Created by user on 05/09/22.
//

#include "libalgebra_lite/polynomial_basis.h"

#include <mutex>
#include <unordered_map>

namespace lal {



std::shared_ptr<const polynomial_basis> basis_registry<polynomial_basis>::get()
{
    static const std::shared_ptr<const polynomial_basis> basis(new polynomial_basis);
    return basis;
}

}
