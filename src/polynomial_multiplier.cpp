//
// Created by user on 30/08/22.
//


#include "libalgebra_lite/polynomial.h"

#include <map>

namespace lal {


typename polynomial_multiplier::product_type
polynomial_multiplier::operator()(const key_type& lhs, const key_type& rhs) const
{
    std::map<letter_type, deg_t> tmp(lhs.begin(), lhs.end());

    for (const auto& item : rhs) {
        tmp[item.first] += item.second;
    }
    product_type result;
    result.emplace_back(key_type(tmp.begin(), tmp.end()), 1);

    return result;
}


}
