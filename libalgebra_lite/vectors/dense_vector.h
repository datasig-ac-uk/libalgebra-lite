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

protected:
    scalar_type* p_data;
    const basis_type* p_basis;
    const dimn_t m_dimension;
    const deg_t m_degree = 0;

    constexpr dense_vector_slice(scalar_type* data, const basis_type* basis, dimn_t dimension) noexcept
        : p_data(data), p_basis(basis), m_dimension(dimension)
    {}

    constexpr dense_vector_slice(scalar_type* data, const basis_type* basis, dimn_t dimension, deg_t degree) noexcept
        : p_data(data), p_basis(basis), m_dimension(dimension), m_degree(degree)
    {}

public:

    constexpr dimn_t dimension() const noexcept { return m_dimension; }
    constexpr basis_type& basis() const noexcept { return *p_basis; }

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
            return {nullptr, p_basis, 0}
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


};





template <typename Basis, typename Coefficients, typename Allocator=typename coefficient_traits<Coefficients>::allocator>
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
    void alloc(dimn_t new_size);


public:

    dense_vector() : slice_type(nullptr, nullptr, 0)
    {}

    dense_vector(const slice_type& ref) : slice_type(nullptr, ref.p_basis, 0, ref.m_degree)
    {
        if (ref.p_data && ref.m_dimension) {
            alloc(ref.m_dimension);
            std::copy(ref.p_data, ref.p_data+ref.m_dimension, p_data);
        }
    }

    dense_vector(const basis_type* basis) : slice_type(nullptr, basis, 0)
    {}






};





} // namespace alg

#endif //LIBALGEBRA_LITE_DENSE_VECTOR_H
