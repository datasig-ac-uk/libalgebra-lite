//
// Created by user on 26/07/22.
//

#include "include/libalgebra_lite/hall_set.h"

#include <sstream>
#include <algorithm>
#include <mutex>
#include <unordered_map>
#include <boost/functional/hash.hpp>

namespace lal {


hall_set::hall_set(hall_set::degree_type width, hall_set::degree_type depth)
        : width(width), current_degree(0)
{
    data.reserve(1 + width);
    letters.reserve(width);
    sizes.reserve(2);
    l2k.reserve(width);
    degree_ranges.reserve(2);

    key_type zero_key(0, 0);

    data.push_back({zero_key, zero_key});
    degree_ranges.push_back({0, 1});
    sizes.push_back(0);

    for (letter_type l = 1; l <= width; ++l) {
        key_type key(1, l-1);
        parent_type parents{zero_key, key};
        letters.push_back(l);
        data.push_back(parents);
        reverse_map.insert(std::make_pair(parents, key));
        l2k.push_back(key);
    }

    std::pair<size_type, size_type> range {
            degree_ranges[current_degree].second,
            data.size()
    };

    degree_ranges.push_back(range);
    sizes.push_back(width);
    ++current_degree;

    if (depth > 1) {
        grow_up(depth);
    }
}

void hall_set::grow_up(hall_set::degree_type new_depth)
{

    for (degree_type d = current_degree + 1; d <= new_depth; ++d) {
        size_type k_index = 0;
        for (degree_type e = 1; 2 * e <= d; ++e) {
            letter_type i_lower, i_upper, j_lower, j_upper;
            i_lower = degree_ranges[e].first;
            i_upper = degree_ranges[e].second;
            j_lower = degree_ranges[d - e].first;
            j_upper = degree_ranges[d - e].second;

            for (letter_type i = i_lower; i < i_upper; ++i) {
                key_type ik(e, i - i_lower);
                for (letter_type j = std::max(j_lower, i + 1); j < j_upper; ++j) {
                    key_type jk(d-e, j-j_lower);
                    if (data[j].first <= ik) {
                        parent_type parents(ik, jk);
                        data.push_back(parents);
                        reverse_map.insert(std::make_pair(parents, key_type(d, k_index++)));
                    }
                }
            }
        }

        std::pair<size_type, size_type> range;
        range.first = degree_ranges[current_degree].second;
        range.second = data.size();
        degree_ranges.push_back(range);
        // The hall set contains an entry for the "god element" 0,
        // so subtract one from the size.
        sizes.push_back(data.size() - 1);

        ++current_degree;
    }
}

hall_set::key_type hall_set::key_of_letter(let_t let) const noexcept
{
    return typename hall_set::key_type {1, letters[let]};
}
hall_set::size_type hall_set::size(deg_t deg) const noexcept
{
    if (deg < sizes.size()) {
        return sizes[deg];
    }
    return sizes.back();
}
hall_set::hall_set(const hall_set& existing, hall_set::degree_type deg)
{

}
hall_set::size_type hall_set::size_of_degree(deg_t) const noexcept
{
    return 0;
}
hall_set::letter_type hall_set::get_letter(dimn_t idx) const noexcept
{
    return 0;
}
typename hall_set::find_result hall_set::find(hall_set::parent_type parent) const noexcept
{
    find_result result;
    result.it = reverse_map.find(parent);
    result.found = (result.it != reverse_map.end());
    return result;
}
bool hall_set::letter(const hall_set::key_type &key) const noexcept
{
    return std::find(letters.begin(), letters.end(), key.index() + 1) != letters.end();
}
const hall_set::parent_type &hall_set::operator[](const hall_set::key_type &key) const noexcept
{
    return data[key.index() + size(key.degree()-1)];
}
const hall_set::key_type &hall_set::operator[](const hall_set::parent_type &parent) const
{
    auto found = reverse_map.find(parent);
    if (found != reverse_map.end()) {
        return found->second;
    }
    return root_element;
}

constexpr typename hall_set::key_type hall_set::root_element;
constexpr typename hall_set::parent_type hall_set::root_parent;

typename hall_set::find_result hall_basis::find(hall_basis::parent_type parents) const noexcept
{
    return p_hallset->find(parents);
}

template class basis_registry<hall_basis>

} // namespace lal
