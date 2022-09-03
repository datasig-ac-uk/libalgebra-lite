//
// Created by user on 01/09/22.
//


#include "include/libalgebra_lite/polynomial_basis.h"


using namespace lal;


deg_t lal::monomial::degree() const noexcept
{
    deg_t result{0};
    for (const auto& item : m_data) {
        result += item.second;
    }
    return result;
}
