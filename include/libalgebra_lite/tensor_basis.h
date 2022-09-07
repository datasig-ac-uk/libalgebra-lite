//
// Created by user on 25/07/22.
//

#ifndef LIBALGEBRA_LITE_TENSOR_BASIS_H
#define LIBALGEBRA_LITE_TENSOR_BASIS_H

#include "implementation_types.h"
#include "libalgebra_lite_export.h"

#include <vector>
#include <cassert>
#include <string>
#include <iosfwd>

#include "registry.h"
#include "index_key.h"
#include "basis_traits.h"

namespace lal {



class  tensor_basis
{
    deg_t m_width;
    deg_t m_depth;

    std::vector<dimn_t> m_powers;
    std::vector<dimn_t> m_sizes;

public:
    using key_type = index_key<>;

    tensor_basis(deg_t width, deg_t depth);

    deg_t width() const noexcept { return m_width; }
    deg_t depth() const noexcept { return m_depth; }

    static constexpr deg_t degree(const key_type& arg) noexcept
    { return deg_t(arg.degree()); }

    key_type lparent(const key_type& arg) const noexcept
    {
        auto degree = arg.degree();
        return key_type(1, arg.degree() / m_powers[degree]);
    }
    key_type rparent(const key_type& arg) const noexcept
    {
        auto degree = arg.degree();
        return key_type(degree-1, arg.index() % m_powers[degree]);
    }

    std::string key_to_string(const key_type&) const noexcept;
    std::ostream& stream_out(std::ostream&, const key_type&) const noexcept;

    static key_type key_of_letter(let_t letter) noexcept
    { return key_type(1, letter); }
    let_t first_letter(const key_type& arg) const noexcept
    { return let_t(lparent(arg).index()); }
    static bool letter(const key_type& arg) noexcept
    { return arg.degree() == 1; }

    dimn_t start_of_degree(deg_t deg) const noexcept
    {
        if (deg == 0) {
            return 0;
        } else {
            return m_sizes[deg-1];
        }
    }
    dimn_t size(int i) const noexcept
    {
        if (i >= 0) {
            return m_sizes[i];
        } else {
            return m_sizes[m_depth];
        }
    }

    const std::vector<dimn_t>& powers() const noexcept
    { return m_powers; }

};


extern template class LIBALGEBRA_LITE_EXPORT basis_registry<tensor_basis>;


} // namespace lal



#endif //LIBALGEBRA_LITE_TENSOR_BASIS_H
