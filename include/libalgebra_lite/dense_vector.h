//
// Created by user on 23/07/22.
//

#ifndef LIBALGEBRA_LITE_DENSE_VECTOR_H
#define LIBALGEBRA_LITE_DENSE_VECTOR_H

#include "implementation_types.h"
#include "basis_traits.h"
#include "coefficients.h"
#include "vector_traits.h"

#include <memory>
#include <type_traits>

namespace lal {


namespace dtl {

template <typename Basis, typename Scalar>
class dense_vector_const_iterator;

template <typename Basis, typename Scalar>
class dense_vector_iterator;

} // namespace dtl

template <typename Basis>
class key_range;


template <typename Basis, typename Coefficients,
          template <typename, typename...> class VectorType,
          typename... Args>
class dense_vector_base {
    using coeff_traits = coefficient_trait<Coefficients>;
    using basis_traits = basis_trait<Basis>;
public:

    using coefficient_ring = typename coeff_traits::coefficient_ring;
    using scalar_type = typename coeff_traits::scalar_type;
    using rational_type = typename coeff_traits::rational_type;

private:
    using storage_type = VectorType<scalar_type, Args...>;

    const Basis* p_basis;
    storage_type m_storage;

public:
    using basis_type        = Basis;
    using size_type         = typename storage_type::size_type;
    using difference_type   = typename storage_type::difference_type;
    using iterator          = dtl::dense_vector_iterator<Basis, scalar_type>;
    using const_iterator    = dtl::dense_vector_const_iterator<Basis, scalar_type>;
    using pointer           = typename storage_type::pointer;
    using const_pointer     = typename storage_type::const_pointer;
    using reference         = typename storage_type::reference;
    using const_reference   = typename storage_type::const_reference;

    dense_vector_base(const Basis* basis, std::initializer_list<scalar_type> args)
            : p_basis(basis), m_storage(args)
    { assert(basis != nullptr); }

    template<typename InputIt>
    dense_vector_base(const Basis* basis, InputIt begin, InputIt end)
            : m_storage(begin, end)
    { assert(basis != nullptr); }

    explicit dense_vector_base(const Basis* basis, size_type n)
        : p_basis(basis), m_storage(n)
    { assert(basis != nullptr); }

    template <typename S>
    explicit dense_vector_base(const Basis* basis, size_type n, const S& val)
        : p_basis(basis), m_storage(n, val)
    { assert(basis != nullptr); }

private:

    size_type adjust_size(size_type n) const noexcept
    {
        return std::min(basis_traits::max_dimension(*p_basis), n);
    }

public:

    void reserve(size_type n)
    {
        m_storage.reserve(adjust_size(n));
    }
    void resize(size_type n)
    {
        m_storage.resize(adjust_size(n));
    }

    template <typename S>
    void resize(size_type n, const S& val)
    {
        m_storage.resize(adjust_size(n), val);
    }

    constexpr size_type size() const noexcept { return m_storage.size(); }
    constexpr bool empty() const noexcept { return m_storage.empty(); }

    template <typename Index>
    reference operator[](Index idx) noexcept
    {
        return m_storage[idx];
    }

    template <typename Index>
    const_reference operator[](Index idx) const noexcept
    { return m_storage[idx]; }

    pointer as_mut_ptr() noexcept { return m_storage.data(); }
    const_pointer as_ptr() const noexcept { return m_storage.data(); }

    void clear() noexcept(noexcept(m_storage.clear())) { m_storage.clear(); }

    // these need to be implemented in terms of kernels.

    template <typename UnaryOp>
    dense_vector_base unary_op(UnaryOp op) const
    {
        dense_vector_base result;
        result.reserve(size());

        const auto begin = m_storage.begin();
        const auto end = m_storage.end();

        for (auto it = begin; it != end; ++it) {
            result.m_storage.emplace_back(op(*it));
        }

        return result;
    }

    template <typename UnaryOp>
    dense_vector_base& inplace_unary_op(UnaryOp op)
    {
        const auto begin = m_storage.begin();
        const auto end = m_storage.end();
        for (auto it=begin; it != end; ++it) {
            op(*it);
        }
        return *this;
    }

    template <typename BinaryOp>
    dense_vector_base binary_op(const dense_vector_base& arg, BinaryOp op)
    {
        dense_vector_base result;

        const difference_type lhs_size(size());
        const difference_type rhs_size(arg.size());

        result.reserve(std::max(lhs_size, rhs_size));

        const auto mid = std::min(lhs_size, rhs_size);
        const auto& zero = coefficient_ring::zero();

        for (difference_type i=0; i<mid; ++i) {
            result.m_storage.emplace_back(op(m_storage[i], arg.m_storage[i]));
        }

        for (auto i=mid; i<lhs_size; ++i) {
            result.emplace_back(op(m_storage[i], zero));
        }

        for (auto i=mid; i<rhs_size; ++i) {
            result.emplace_back(op(zero, arg.m_storage[i]));
        }

        return result;
    }

    template <typename InplaceBinaryOp>
    dense_vector_base& inplace_binary_op(const dense_vector_base& rhs, InplaceBinaryOp op)
    {
        const difference_type lhs_size(size());
        const difference_type rhs_size(rhs.size());

        if (rhs_size > lhs_size) {
            resize(rhs_size);
        }

        const auto& zero = coefficient_ring::zero();
        const auto mid = std::min(lhs_size, rhs_size);

        for (difference_type i=0; i<mid; ++i) {
            op(m_storage[i], rhs.m_storage[i]);
        }

        for (auto i=mid; i<lhs_size; ++i) {
            op(m_storage[i], zero);
        }

        difference_type rhs_max(std::min(size(), rhs_size()));
        for (auto i=mid; i<rhs_max; ++i) {
            op(m_storage, rhs.m_storage[i]);
        }

        return *this;
    }

};


template <typename Basis, typename Coefficients>
using dense_vector = dense_vector_base<Basis, Coefficients, std::vector>;

namespace dtl {



template <typename KeyRef, typename ScaRef>
class dense_iterator_item
{
    template <typename B, typename C>
    friend class dense_vector_iterator;

    template <typename B, typename C>
    friend class dense_vector_const_iterator;

    KeyRef m_key;
    ScaRef m_sca;

    dense_iterator_item(KeyRef key, ScaRef sca)
        : m_key(key), m_sca(sca)
    {}

public:
    KeyRef key() noexcept { return m_key; }
    ScaRef value() noexcept { return m_sca; }
};

template <typename Basis, typename Coefficients>
class dense_vector_iterator
{
    using basis_traits = basis_trait<Basis>;
    using key_type = typename basis_traits::key_type;
    using coeff_traits = coefficient_trait<Coefficients>;
    using scalar_type = typename coeff_traits::scalar_type;

    const Basis* p_basis = nullptr;
    scalar_type* p_data = nullptr;
    key_type m_key;


public:

    using value_type = dtl::dense_iterator_item<const key_type&, scalar_type&>;
    using reference = value_type;
    using pointer = value_type;
    using iterator_category = std::forward_iterator_tag;


    dense_vector_iterator() = default;

    dense_vector_iterator(const Basis* basis, scalar_type* data)
        : p_basis(basis), p_data(data), m_key()
    {}

    dense_vector_iterator& operator++()
    {
        p_data++;
        m_key++;
        return *this;
    }

    reference operator*() noexcept
    {
        return {m_key, *p_data};
    }

    pointer operator->() noexcept
    {
        return {m_key, *p_data};
    }

    bool operator==(const dense_vector_iterator& other) const noexcept
    {
        return p_data == other.p_data;
    }

    bool operator!=(const dense_vector_iterator& other) const noexcept
    {
        return p_data != other.p_data;
    }
};



template <typename Basis, typename Coefficients>
class dense_vector_const_iterator
{
    using basis_traits = basis_trait<Basis>;
    using key_type = typename basis_traits::key_type;
    using coeff_traits = coefficient_trait<Coefficients>;
    using scalar_type = typename coeff_traits::scalar_type;

    const Basis* p_basis = nullptr;
    const scalar_type* p_data = nullptr;
    key_type m_key;


public:

    using value_type = dtl::dense_iterator_item<const key_type&, const scalar_type&>;
    using reference = value_type;
    using pointer = value_type;
    using iterator_category = std::forward_iterator_tag;

    dense_vector_const_iterator(const Basis* basis, scalar_type* data)
        : p_basis(basis), p_data(data), m_key()
    {}

    dense_vector_const_iterator& operator++()
    {
        p_data++;
        m_key++;
        return *this;
    }

    reference operator*() noexcept
    {
        return {m_key, *p_data};
    }

    pointer operator->() noexcept
    {
        return {m_key, *p_data};
    }

    bool operator==(const dense_vector_const_iterator& other) const noexcept
    {
        return p_data == other.p_data;
    }

    bool operator!=(const dense_vector_const_iterator& other) const noexcept
    {
        return p_data != other.p_data;
    }
};


} // namespace dtl



} // namespace lal

#endif //LIBALGEBRA_LITE_DENSE_VECTOR_H
