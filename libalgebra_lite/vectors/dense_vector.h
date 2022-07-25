//
// Created by user on 23/07/22.
//

#ifndef LIBALGEBRA_LITE_DENSE_VECTOR_H
#define LIBALGEBRA_LITE_DENSE_VECTOR_H

#include <libalgebra_lite/implementation_types.h>

#include <memory>
#include <type_traits>

namespace alg {

template <typename Basis>
struct basis_trait;

template <typename Coefficients>
struct coefficient_trait;

namespace dtl {

template <typename Basis, typename Scalar>
class dense_vector_const_iterator;

template <typename Basis, typename Scalar>
class dense_vector_iterator;

} // namespace dtl

template <typename Keytype>
class key_range;

template <typename Basis, typename Coefficients, typename Allocator=typename coefficient_trait<Coefficients>::default_alloc>
class dense_vector;

template <typename Basis, typename Coefficients>
class dense_vector_slice
{
protected:
    using basis_traits = basis_trait<Basis>;
    using coeff_traits = coefficient_trait<Coefficients>;

public:

    using basis_type = Basis;
    using key_type = typename basis_traits::key_type;

    using coefficient_ring = typename coeff_traits::coeffificient_ring;
    using scalar_type = typename coeff_traits::scalar_type;
    using rational_type = typename coeff_traits::rational_type;

    using const_iterator = dtl::dense_vector_const_iterator<Basis, scalar_type>;
    using iterator = dtl::dense_vector_iterator<Basis, scalar_type>;

    using dense_vector = dense_vector<Basis, Coefficients>;

protected:
    scalar_type* p_data;
    const basis_type* p_basis;
    dimn_t m_dimension;
    deg_t m_degree = 0;

    constexpr dense_vector_slice(scalar_type* data, const basis_type* basis, dimn_t dimension) noexcept
        : p_data(data), p_basis(basis), m_dimension(dimension)
    {}

    constexpr dense_vector_slice(scalar_type* data, const basis_type* basis, dimn_t dimension, deg_t degree) noexcept
        : p_data(data), p_basis(basis), m_dimension(dimension), m_degree(degree)
    {}

public:

    constexpr dimn_t dimension() const noexcept { return m_dimension; }
    constexpr basis_type& basis() const noexcept { return *p_basis; }
    dimn_t capacity() const noexcept { return basis_traits::max_dimension(*p_basis); }

    const scalar_type& operator[](const key_type& k) const noexcept
    {
        return p_data[basis_traits::key_to_index(*p_basis, k)];
    }
    scalar_type& operator[](const key_type& k) noexcept
    { return p_data[basis_traits::key_to_index(*p_basis, k)]; }

    template <typename KeyType>
    std::enable_if_t<std::is_constructible<key_range<key_type>, const KeyType&>::value, const dense_vector_slice>
    operator[](const KeyType& k) const noexcept
    {
        if (!p_data) {
            return {nullptr, p_basis, 0};
        }

        const key_range<key_type> range(k);
        const auto idx_begin = basis_traits::key_to_index(*p_basis, range.begin());
        const auto dimension = std::min(basis_traits::key_to_index(*p_basis, range.end()), m_dimension);
        return {
            p_data + idx_begin,
            p_basis,
            dimension
        };
    }

    template <typename KeyType>
    std::enable_if_t<std::is_constructible<key_range<key_type>, const KeyType&>::value, dense_vector_slice>
    operator[](const KeyType& k) noexcept
    {
        if (!p_data) {
            return {nullptr, p_basis, 0};
        }

        const key_range<key_type> range(k);
        const auto idx_begin = basis_traits::key_to_index(*p_basis, range.begin());
        const auto dimension = std::min(basis_traits::key_to_index(*p_basis, range.end()), m_dimension);
        return {
                p_data+idx_begin,
                p_basis,
                dimension
        };
    }

    const scalar_type* as_ptr() const noexcept { return p_data; };
    scalar_type* as_mut_ptr() noexcept { return p_data; }

    constexpr const_iterator begin() const noexcept { return {p_basis, p_data}; }
    constexpr const_iterator end() const noexcept { return {p_basis, p_data + m_dimension}; }
    constexpr iterator begin() noexcept { return {p_basis, p_data}; }
    constexpr iterator end() noexcept { return {p_basis, p_data + m_dimension}; }


protected:

    void resize(dimn_t new_size, deg_t degree_hint=0) noexcept
    {
        const auto adjusted = basis_traits::get_dimension_degree(*p_basis, new_size, degree_hint);
        if (adjusted.first > m_dimension) {
            m_dimension = new_size;
            m_degree = adjusted.second;
        }
    }

    void destroy(dimn_t new_dim) noexcept(std::is_nothrow_destructible<scalar_type>::value)
    {
    }

    template <typename... Args>
    void construct_in_place(idimn_t index, Args&&... args) noexcept(noexcept(scalar_type(std::forward<Args>(args)...)))
    {
        ::new (p_data + index) scalar_type(std::forward<Args>(args)...);
    }

    template <typename UpdateFn>
    void update_in_place(dimn_t count, UpdateFn fn) noexcept(noexcept(fn(std::declval<scalar_type&>())))
    {
        for (dimn_t cnt=0; cnt<count; ++cnt) {
            fn(p_data[cnt]);
        }
    }

    void fill(dimn_t new_dim, const scalar_type& new_val) noexcept(std::is_nothrow_copy_constructible<scalar_type>::value)
    {
        std::uninitialized_fill(p_data + m_dimension, p_data + new_dim, new_val);
    }

public:

    void clear() noexcept(std::is_nothrow_destructible<scalar_type>::value)
    {
        destroy(0);
    }

    template <typename UnaryOp>
    dense_vector
    unary_op(UnaryOp op) const noexcept(noexcept(op(std::declval<const scalar_type&>())))
    {
        dense_vector result(p_basis);
        const scalar_type* data = p_data;
        for (idimn_t i=0; i<static_cast<idimn_t>(m_dimension); ++i) {
            construct_in_place(i, op(*(data++)));
        }
        result.resize(m_dimension);
    }

    template <typename UnaryOp>
    dense_vector_slice& inplace_unary_op(UnaryOp op) noexcept(noexcept(op(std::declval<scalar_type&>())))
    {
        update_in_place(m_dimension, op);
        return *this;
    }


    template <typename BinOp>
    dense_vector binary_op(const dense_vector_slice& rhs, BinOp op) const
        noexcept(noexcept(op(std::declval<const scalar_type&>(), std::declval<const scalar_type&>())))
    {
        dense_vector result(p_basis);

        if (basis_traits::compatible(*p_basis, *rhs.p_basis)) {
            const auto mid = std::min(m_dimension, rhs.m_dimension);
            const scalar_type* lptr = p_data, *rptr = rhs.p_data;
            for (idimn_t i=0; i<static_cast<idimn_t>(m_dimension); ++i) {
                construct_in_place(i, op(*(lptr++), *(rptr++)));
            }
            auto new_dim = mid;

            if (m_dimension > mid) {
                const scalar_type* dptr = p_data;
                const auto zero = coeff_traits::zero();
                for (idimn_t i=mid; i<static_cast<idimn_t>(m_dimension); ++i) {
                    construct_in_place(i, op(*(dptr++), zero));
                }
                new_dim = m_dimension;
            } else if (rhs.m_dimension > mid) {
                const scalar_type* dptr = rhs.p_data;
                const auto max = std::max(rhs.m_dimension, basis_traits::max_dimension(*p_basis));
                const auto zero = coeff_traits::zero();
                for (idimn_t i=m_dimension; i<static_cast<idimn_t>(max); ++i) {
                    construct_in_place(i, op(zero, *(dptr++)));
                }
            }

            resize(new_dim, m_degree);
        }
        return result;
    }

    /*
     * There is a small optimisation we can apply here if we assume that the binop is "addition-like"
     * in that if the rhs elements are zero then the corresponding lhs elements left untouched.
     * This is useful because it means we can ignore the case where the rhs dimension is smaller than
     * the lhs dimension.
     */

    template <typename BinaryOp>
    dense_vector_slice& inplace_binop(const dense_vector_slice& rhs, BinaryOp op)
        noexcept(noexcept(op(std::declval<scalar_type&>(), std::declval<const scalar_type&>())))
    {
        if (basis_traits::compatible(*p_basis, rhs.*p_basis)) {
            const auto mid = std::min(m_dimension, rhs.m_dimension);
            for (idimn_t i=0; i<static_cast<idimn_t>(mid); ++i) {
                op(p_data[i], rhs.p_data[i]);
            }
            auto new_dim = mid;

            if (m_dimension > mid) {
                const auto zero = coeff_traits::zero();
                for (idimn_t i=mid; i<m_dimension; ++i) {
                    op(p_data[i], zero);
                }
                new_dim = m_dimension;
            } else if (rhs.m_dimension > mid) {
                const auto max = std::min(rhs.m_dimension, basis_traits::max_dimension(*p_basis));
                /*
                 * This is an interesting problem because ideally I don't want to effectively
                 * iterate through the data twice to apply the operation. However, when the
                 * scalar_type is trivially constructible, filling is almost certainly what
                 * one wants to do. (A memset call will be vastly superior to looping.)
                 */
                fill(max, coeff_traits::zero());
                for (idimn_t i=mid; i<static_cast<idimn_t>(max); ++i) {
                    op(p_data[i], rhs.p_data[i]);
                }
                new_dim = max;
            }
            resize(new_dim, m_degree);
        }
        return *this;
    }

    template <typename TernaryOp>
    dense_vector_slice& ternary_op(const dense_vector_slice& lhs, const dense_vector_slice& rhs, TernaryOp op)
        noexcept(noexcept(op(
                std::declval<scalar_type*>(),
                std::declval<const scalar_type*>(),
                std::declval<const scalar_type*>())))
    {
        op(p_data, lhs.p_data, rhs.p_data);
        return *this;
    }


};





template <typename Basis, typename Coefficients, typename Allocator>
class dense_vector : public dense_vector_slice<Basis, Coefficients>
{
protected:
    using slice_type = dense_vector_slice<Basis, Coefficients>;
    using typename slice_type::basis_traits;
    using typename slice_type::coeff_traits;

public:
    using typename slice_type::basis_type;
    using typename slice_type::key_type;
    using typename slice_type::coefficient_ring;
    using typename slice_type::scalar_type;
    using typename slice_type::rational_type;

    using typename slice_type::iterator;
    using typename slice_type::const_iterator;

protected:

    using alloc_traits = std::allocator_traits<Allocator>;
    using allocator_type = Allocator;
    allocator_type m_alloc;

    void resize(dimn_t new_size);
    void alloc_and_copy(const scalar_type* existing, dimn_t count)
    {
        if (p_basis) {
            auto capacity = basis_traits::max_dimension(*p_basis);
            p_data = alloc_traits::allocate(m_alloc, capacity);

            if (existing) {
                std::uninitialized_copy_n(existing, count, p_data);
            }
        }
    }

    void dealloc();

    using slice_type::p_data;
    using slice_type::p_basis;
    using slice_type::m_dimension;
    using slice_type::m_degree;

public:

    dense_vector() : slice_type(nullptr, nullptr, 0)
    {}

    dense_vector(const slice_type& ref) : slice_type(nullptr, ref.p_basis, ref.m_dimension, ref.m_degree)
    {
        alloc_and_copy(ref.p_data, ref.m_dimension);
    }

    dense_vector(const basis_type* basis) : slice_type(nullptr, basis, 0)
    {
        alloc_and_copy(nullptr, 0);
    }

    ~dense_vector()
    {
        slice_type::destroy(0);
        dealloc();
    }




};





} // namespace alg

#endif //LIBALGEBRA_LITE_DENSE_VECTOR_H
