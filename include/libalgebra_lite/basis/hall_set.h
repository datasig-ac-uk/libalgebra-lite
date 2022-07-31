//
// Created by user on 26/07/22.
//

#ifndef LIBALGEBRA_LITE_HALL_SET_H
#define LIBALGEBRA_LITE_HALL_SET_H

#include <libalgebra_lite/implementation_types.h>
#include <libalgebra_lite/basis/index_key.h>
#include "libalgebra_lite_export.h"
#include <utility>
#include <vector>

#include <boost/container/flat_map.hpp>
#include <mutex>
#include <memory>

namespace alg {

class hall_set {
public:
    using letter_type = let_t;
    using degree_type = deg_t;
    using key_type = index_key<>;
    using size_type = dimn_t;
    using parent_type = std::pair<key_type, key_type>;
private:
    using data_type = std::vector<parent_type>;
    using reverse_map_type = boost::container::flat_map<parent_type, key_type>;
    using l2k_map_type = std::vector<key_type>;
    using size_vector_type = std::vector<size_type>;
    using degree_range_map_type = std::vector<std::pair<size_type, size_type>>;

    degree_type width;
    degree_type current_degree;

    std::vector<letter_type> letters;
    data_type data;
    reverse_map_type reverse_map;
    l2k_map_type l2k;
    degree_range_map_type degree_ranges;
    size_vector_type sizes;

public:
    static constexpr key_type root_element {0, 0};
    static constexpr parent_type root_parent {root_element, root_element};

    explicit hall_set(degree_type width, degree_type depth=1);
    explicit hall_set(const hall_set& existing, degree_type deg);

    void grow_up(degree_type deg);


    key_type key_of_letter(let_t) const noexcept;
    size_type size(deg_t) const noexcept;
    bool letter(const key_type&) const noexcept;
    letter_type get_letter(dimn_t idx) const noexcept;

    const parent_type &operator[](const key_type&) const noexcept;
    const key_type& operator[](const parent_type&) const;
};



class LIBALGEBRA_LITE_EXPORT hall_basis
{
    std::shared_ptr<hall_set> p_hallset;
    deg_t m_width;
    deg_t m_depth;

public:

    using key_type = typename hall_set::key_type;
    using parent_type = typename hall_set::parent_type;


    deg_t width() const noexcept { return m_width; }
    deg_t depth() const noexcept { return m_depth; }

    static constexpr deg_t degree(const key_type& arg) noexcept
    { return deg_t(arg.degree()); }

    parent_type parents(const key_type& arg) const noexcept
    { return (*p_hallset)[arg]; }
    key_type lparent(const key_type& arg) const noexcept
    { return parents(arg).first; }
    key_type rparent(const key_type& arg) const noexcept
    { return parents(arg).second; }
    key_type key_of_letter(let_t letter) const noexcept
    { return p_hallset->key_of_letter(letter); }
    let_t first_letter(const key_type& key) const noexcept
    { return p_hallset->get_letter((*p_hallset)[key].first.index()); }
    dimn_t size(int deg) const noexcept
    {
        return p_hallset->size(deg < 0 ? static_cast<deg_t>(deg) : m_depth);
    }
    dimn_t start_of_degree(deg_t deg) const noexcept
    {
        return (deg == 0) ? 0 : p_hallset->size(deg-1);
    }


};




} // alg

#endif //LIBALGEBRA_LITE_HALL_SET_H
