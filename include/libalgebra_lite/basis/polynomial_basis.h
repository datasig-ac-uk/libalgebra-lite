//
// Created by user on 28/08/22.
//

#ifndef LIBALGEBRA_LITE_POLYNOMIAL_BASIS_H
#define LIBALGEBRA_LITE_POLYNOMIAL_BASIS_H

#include <libalgebra_lite/implementation_types.h>
#include <libalgebra_lite/utilities/packed_integer.h>

#include <functional>

#include <numeric>
#include <boost/container/small_vector.hpp>
#include <boost/container/flat_map.hpp>

namespace lal {

struct polynomial_basis
{
    using letter_type = dtl::packed_integer<dimn_t, char>;
    using small_vec = boost::container::small_vector<std::pair<letter_type, deg_t>, 1>;
    using key_type = boost::container::flat_map<letter_type, deg_t, std::less<>, small_vec>;

    static key_type key_of_letter(letter_type letter)
    {
        key_type result;
        result[letter] = 1;
        return result;
    }

    static deg_t degree(const key_type& key)
    {
        return std::accumulate(key.begin(), key.end(), 0,
                [](const auto& curr, const auto& key) { return curr + key.second; });
    }

    struct key_order
    {
        bool operator()(const key_type& lhs, const key_type& rhs) const;
    };

};


} // namespace lal

#endif //LIBALGEBRA_LITE_POLYNOMIAL_BASIS_H
