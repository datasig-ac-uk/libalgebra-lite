//
// Created by user on 05/09/22.
//

#include "libalgebra_lite/polynomial_basis.h"

#include <mutex>
#include <unordered_map>
#include <ostream>

namespace lal {


std::ostream& operator<<(std::ostream& os, const monomial& arg) noexcept
{
    bool first = true;
    for (const auto& item : arg) {
        os << item.first << '^' << item.second;
        if (first) {
            first = false;
        } else {
            os << ' ';
        }
    }
    return os;
}
deg_t monomial::degree() const noexcept
{
    return std::accumulate(m_data.begin(), m_data.end(), 0,
            [](const auto& curr, const auto& key) { return curr + key.second; });
}
std::ostream& polynomial_basis::print_key(std::ostream& os, const polynomial_basis::key_type& key) const
{
    return os << key;
}

std::shared_ptr<const polynomial_basis> basis_registry<polynomial_basis>::get()
{
    static const std::shared_ptr<const polynomial_basis> basis(new polynomial_basis);
    return basis;
}

}
