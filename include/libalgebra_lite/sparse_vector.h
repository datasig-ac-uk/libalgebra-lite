//
// Created by user on 26/07/22.
//

#ifndef LIBALGEBRA_LITE_SPARSE_VECTOR_H
#define LIBALGEBRA_LITE_SPARSE_VECTOR_H

#include "implementation_types.h"

#include <iterator>
#include <map>
#include <utility>
#include <type_traits>

#include <boost/container/flat_map.hpp>

#include "basis_traits.h"
#include "coefficients.h"


namespace lal {
namespace dtl {

#define LAL_MUTABLE_REF_iOP(OP)                                             \
    template <typename Scalar>                                              \
    Self& operator OP(Scalar arg) noexcept(noexcept(m_tmp OP arg))          \
    {                                                                       \
        m_tmp OP arg;                                                       \
        return *this;                                                       \
    }

#define LAL_MUTABLE_REF_COMPARE(OP)                                         \
    template <typename Scalar>                                              \
    bool operator OP(Scalar arg) noexcept(noexcept(m_tmp OP arg))           \
    {                                                                       \
        return m_tmp OP arg;                                                \
    }


template <typename MapType, typename KeyType>
class sparse_mutable_reference
{
    using iterator_type = typename MapType::iterator;

    MapType& m_map;
    iterator_type m_it;
    typename MapType::mapped_type m_tmp;
    KeyType m_key;

    using Self = sparse_mutable_reference;

public:

    using key_type = KeyType;
    using scalar_type = typename MapType::mapped_type;

    sparse_mutable_reference(MapType& map, iterator_type it)
        : m_map(map), m_it(it), m_key(it->first), m_tmp(it->second)
    {
        assert(it != m_map.end());
    }

    sparse_mutable_reference(MapType& map, KeyType key)
        : m_map(map), m_key(key), m_it(map.find(key)), m_tmp(0)
    {
        if (m_it != m_map.end()) {
            m_tmp = m_it->second;
        }
    }

    ~sparse_mutable_reference()
    {
        if (m_tmp != scalar_type(0)) {
            if (m_it != m_map.end()) {
                m_it->second = m_tmp;
            } else {
                m_map[m_key] = m_tmp;
            }
        } else if (m_it != m_map.end()) {
            m_map.erase(m_it);
        }
    }

    operator const scalar_type& () const noexcept // NOLINT(google-explicit-constructor)
    {
        return m_tmp;
    }

    LAL_MUTABLE_REF_iOP(=)
    LAL_MUTABLE_REF_iOP(+=)
    LAL_MUTABLE_REF_iOP(-=)
    LAL_MUTABLE_REF_iOP(*=)
    LAL_MUTABLE_REF_iOP(/=)
    LAL_MUTABLE_REF_iOP(<<=)
    LAL_MUTABLE_REF_iOP(>>=)
    LAL_MUTABLE_REF_iOP(|=)
    LAL_MUTABLE_REF_iOP(&=)
    LAL_MUTABLE_REF_iOP(^=)
    LAL_MUTABLE_REF_iOP(%=)

    LAL_MUTABLE_REF_COMPARE(==)
    LAL_MUTABLE_REF_COMPARE(!=)
    LAL_MUTABLE_REF_COMPARE(<)
    LAL_MUTABLE_REF_COMPARE(<=)
    LAL_MUTABLE_REF_COMPARE(>)
    LAL_MUTABLE_REF_COMPARE(>=)

};

#undef LAL_MUTABLE_REF_COMPARE
#undef LAL_MUTABLE_REF_iOP

template <typename MapType, typename Iterator, typename Parent>
class sparse_iterator_base
{
protected:
    MapType* p_map = nullptr;
    Iterator m_it;

    using traits = std::iterator_traits<Iterator>;

    using key_type = typename traits::value_type::first_type;
    using scalar_type = typename traits::value_type::second_type;


public:

    sparse_iterator_base(MapType* map, Iterator it)
        : p_map(map), m_it(it)
    {
        assert(map != nullptr);
    }

    sparse_iterator_base(MapType& map, Iterator it)
        : p_map(&map), m_it(it)
    {}

    Parent& operator++() noexcept { ++m_it; return *this; }
    const Parent operator++(int) noexcept
    {
        auto result = Parent(p_map, m_it);
        ++m_it;
        return result;
    }

    Parent& operator*() const noexcept { return *this; }
    Parent* operator->() const noexcept { return this; }

    bool operator==(const sparse_iterator_base& other) const noexcept
    {
        return m_it == other.m_it;
    }
    bool operator!=(const sparse_iterator_base& other) const noexcept
    {
        return m_it != other.m_it;
    }
};



template <typename MapType, typename Iterator>
class sparse_iterator;

template <typename MapType>
class sparse_iterator<MapType, typename MapType::iterator>
        : public sparse_iterator_base<
                MapType,
                typename MapType::iterator,
                sparse_iterator<MapType, typename MapType::iterator>>
{
    using base = sparse_iterator_base<MapType, typename MapType::iterator, sparse_iterator>;
public:
    using reference = sparse_mutable_reference<MapType, const typename base::key_type&>;
    using base::base;

    const typename base::key_type& key() const noexcept
    {
        assert(base::p_map != nullptr);
        return base::m_it->first;
    }

    reference value() const noexcept
    {
        assert(base::p_map != nullptr);
        return reference(*base::p_map, base::m_it);
    }
};

template <typename MapType>
class sparse_iterator<MapType, typename MapType::const_iterator>
        : public sparse_iterator_base<MapType, typename MapType::const_iterator,
        sparse_iterator<MapType, typename MapType::const_iterator>>
{
    using base = sparse_iterator_base<MapType, typename MapType::const_iterator, sparse_iterator>;
public:

    using base::base;

    const typename base::key_type& key() const noexcept
    { return base::m_it->first; }
    const typename base::scalar_type& value() const noexcept
    { return base::m_it->second; }
};




} // namespace dtl


template <typename Basis, typename Coefficients>
class sparse_vector {
    using basis_traits = basis_trait<Basis>;
    using coeff_traits = coefficient_trait<Coefficients>;
public:

    using basis_type = Basis;
    using key_type = typename basis_traits::key_type;
    using coefficient_ring = typename coeff_traits::coefficient_ring;

    using scalar_type = typename coeff_traits::scalar_type;
    using rational_type = typename coeff_traits::rational_type;
private:

    using map_type = boost::container::flat_map<key_type, scalar_type>;
    map_type m_data;
    const basis_type* p_basis;
    deg_t m_degree = 0;

protected:
    sparse_vector(const basis_type* basis, map_type&& arg)
        : m_data(arg), p_basis(basis)
    {
        assert(p_basis != nullptr);
    }

public:

    using reference = dtl::sparse_mutable_reference<map_type, key_type>;
    using const_reference = const scalar_type&;

    using iterator = dtl::sparse_iterator<map_type, typename map_type::iterator>;
    using const_iterator = dtl::sparse_iterator<const map_type, typename map_type::const_iterator>;

    explicit sparse_vector(const basis_type* basis) : p_basis(basis)
    {
        assert(p_basis != nullptr);
    }


    constexpr dimn_t size() const noexcept { return m_data.size(); }
    constexpr bool empty() const noexcept { return m_data.empty(); }
    constexpr dimn_t dimension() const noexcept { return size(); }
    constexpr const basis_type& basis() const noexcept { return *p_basis; }
    dimn_t capacity() const noexcept { return basis_traits::max_dimension(*p_basis); }

    iterator begin() noexcept { return {m_data, m_data.begin()}; }
    iterator end() noexcept { return {m_data, m_data.end()}; }
    const_iterator begin() const noexcept { return {m_data, m_data.begin()}; }
    const_iterator end() const noexcept { return {m_data, m_data.end()}; }

    const_iterator cbegin() const noexcept { return begin(); }
    const_iterator cend() const noexcept { return end(); }

    const_reference operator[](const key_type& key) const noexcept
    {
        auto val = m_data.find(key);
        if (val != m_data.end()) {
            return val->second;
        }
        return coefficient_ring::zero();
    }

    reference operator[](const key_type& key) noexcept
    {
        return reference(m_data, key);
    }

    void clear() noexcept
    { m_data.clear(); }

    template <typename UnaryOp>
    sparse_vector unary_op(UnaryOp op) const
    {
        map_type data;
        data.reserve(m_data.size());
        for (const auto& item : m_data) {
            data.emplace(item.first, op(item.second));
        }
        return {p_basis, std::move(data)};
    }

    template <typename BinOp>
    sparse_vector binary_op(const sparse_vector& rhs, BinOp op) const
    {
        map_type data;
        data.reserve(m_data.size());

        auto lit = m_data.begin();
        auto lend = m_data.end();
        auto rit = rhs.m_data.begin();
        auto rend = rhs.m_data.end();

        scalar_type val;
        const auto& zero = coefficient_ring::zero();

        while (lit != lend && rit != rend) {
            if (lit->first == rit->first
                && (val = op(lit->second, rit->second) != zero)) {
                data[lit->first] = val;
                ++lit;
                ++rit;
            } else if (lit->first < rit->first) {
                val = op(lit->second, zero);
                if (val != zero) {
                    data[lit->first] = val;
                }
                ++lit;
            } else {
                val = op(zero, rit->second);
                if (val != zero) {
                    data[rit->first] = val;
                }
                ++rit;
            }
        }

        return {p_basis, std::move(data)};
    }

    template <typename BinOp>
    sparse_vector& inplace_binop(const sparse_vector& rhs, BinOp op)
    {
        auto new_op = [op](const scalar_type& lval, const scalar_type& rval) {
            auto tmp = lval;
            op(tmp, rval);
            return tmp;
        };
        auto tmp = binary_op(rhs, new_op);
        std::swap(m_data, tmp.m_data);
        return *this;
    }



};


} // namespace lal

#endif //LIBALGEBRA_LITE_SPARSE_VECTOR_H