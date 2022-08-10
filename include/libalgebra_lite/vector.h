//
// Created by user on 08/08/22.
//

#ifndef LIBALGEBRA_LITE_VECTOR_H
#define LIBALGEBRA_LITE_VECTOR_H

#include <libalgebra_lite/implementation_types.h>
#include <libalgebra_lite/vectors/traits.h>
#include <libalgebra_lite/coefficients.h>
#include <libalgebra_lite/basis/traits.h>


#include <memory>


#define LAL_INLINE_ALWAYS __attribute__((always_inline))


namespace alg {

namespace dtl {

struct vector_base_access;

} // namespace dtl


template <typename Basis,
        typename Coefficients,
        template <typename, typename, typename...> class VectorType,
        typename... Args>
class vector
{
    using vector_type = VectorType<Basis, Coefficients, Args...>;
    std::shared_ptr<vector_type> p_impl;
    friend struct dtl::vector_base_access;

protected:
    using basis_traits = basis_trait<Basis>;
    using coeff_traits = coefficient_trait<Coefficients>;

public:
    using basis_type        = Basis;
    using key_type          = typename basis_traits::key_type;
    using coefficient_ring  = typename coeff_traits::coefficient_ring;
    using scalar_type       = typename coeff_traits::scalar_type;
    using rational_type     = typename coeff_traits::rational_type;

    using iterator          = typename vector_type::iterator;
    using const_iterator    = typename vector_type::const_iterator;
    using reference         = typename vector_type::reference;
    using const_reference   = typename vector_type::const_reference;


private:

    vector(vector_type&& arg)
        : p_impl(std::make_shared<vector_type>(std::move(arg)))
    {}

public:

    template <typename... CArgs>
    vector(CArgs&&... args)
        : p_impl(std::make_shared<vector_type>(std::forward<CArgs>(args)...))
    {}

    vector(const basis_type* basis, std::initializer_list<scalar_type> args)
        : p_impl(std::make_shared<vector_type>(basis, args))
    {}

    LAL_INLINE_ALWAYS
    dimn_t size() const noexcept
    {
        return p_impl->size();
    }
    LAL_INLINE_ALWAYS
    dimn_t dimension() const noexcept
    {
        return p_impl->dimension();
    }
    LAL_INLINE_ALWAYS
    dimn_t basis() const noexcept
    {
        return p_impl->basis();
    }
    LAL_INLINE_ALWAYS
    reference operator[](const key_type& key) noexcept
    {
        return (*p_impl)[key];
    }
    LAL_INLINE_ALWAYS
    const_reference operator[](const key_type& key) const noexcept
    {
        return (*p_impl)[key];
    }
    LAL_INLINE_ALWAYS
    void clear()
    { p_impl->clear(); }
    LAL_INLINE_ALWAYS
    iterator begin() noexcept
    { return p_impl->begin(); }
    LAL_INLINE_ALWAYS
    iterator end() noexcept
    { return p_impl->end(); }
    LAL_INLINE_ALWAYS
    const_iterator begin() const noexcept
    { return p_impl->begin(); }
    LAL_INLINE_ALWAYS
    const_iterator end() const noexcept { return p_impl->end(); }

    LAL_INLINE_ALWAYS
    explicit operator vector_type& ()
    { return *p_impl;}
    LAL_INLINE_ALWAYS
    explicit operator const vector_type& ()
    { return *p_impl; }
    LAL_INLINE_ALWAYS
    vector_type& base_vector() noexcept
    { return *p_impl; }
    LAL_INLINE_ALWAYS
    const vector_type& base_vector() const noexcept
    { return *p_impl; }

public:

    template <typename Key, typename Scal>
    vector& add_scal_prod(const Key& key, const Scal& scal);
    template <typename Key, typename Rat>
    vector& add_scal_div(const Key& key, const Rat& scal);
    template <typename Key, typename Scal>
    vector& sub_scal_prod(const Key& key, const Scal& scal);
    template <typename Key, typename Rat>
    vector& sub_scal_div(const Key& key, const Rat& scal);

    template <typename Scal>
    vector& add_scal_prod(const vector& rhs, const Scal& scal);
    template <typename Rat>
    vector& add_scal_div(const vector& rhs, const Rat& scal);
    template <typename Scal>
    vector& sub_scal_prod(const vector& rhs, const Scal& scal);
    template <typename Rat>
    vector& sub_scal_div(const vector& rhs, const Rat& scal);

    template <template <typename, typename, typename...> class AltVecType,
            typename... AltArgs,
            typename Scal>
    vector& add_scal_prod(
            const vector<Basis, Coefficients, AltVecType, AltArgs...>& rhs,
            const Scal& scal
    );
    template <template <typename, typename, typename...> class AltVecType,
            typename... AltArgs,
            typename Rat>
    vector& add_scal_div(
            const vector<Basis, Coefficients, AltVecType, AltArgs...>& rhs,
            const Rat& scal
    );
    template <template <typename, typename, typename...> class AltVecType,
            typename... AltArgs,
            typename Scal>
    vector& sub_scal_prod(
            const vector<Basis, Coefficients, AltVecType, AltArgs...>& rhs,
            const Scal& scal
    );
    template <template <typename, typename, typename...> class AltVecType,
            typename... AltArgs,
            typename Rat>
    vector& sub_scal_div(
            const vector<Basis, Coefficients, AltVecType, AltArgs...>& rhs,
            const Rat& scal
    );

    template <
            typename AltBasis,
            typename AltCoeffs,
            template <typename, typename, typename...> class AltVecType,
            typename... AltArgs,
            typename Scal>
    vector& add_scal_prod(
            const vector<AltBasis, AltCoeffs, AltVecType, AltArgs...>& rhs,
            const Scal& scal
    );
    template <
            typename AltBasis,
            typename AltCoeffs,
            template <typename, typename, typename...> class AltVecType,
            typename... AltArgs,
            typename Rat>
    vector& add_scal_div(
            const vector<AltBasis, AltCoeffs, AltVecType, AltArgs...>& rhs,
            const Rat& scal
    );
        template <
            typename AltBasis,
            typename AltCoeffs,
            template <typename, typename, typename...> class AltVecType,
            typename... AltArgs,
            typename Scal>
    vector& sub_scal_prod(
            const vector<AltBasis, AltCoeffs, AltVecType, AltArgs...>& rhs,
            const Scal& scal
    );
    template <
            typename AltBasis,
            typename AltCoeffs,
            template <typename, typename, typename...> class AltVecType,
            typename... AltArgs,
            typename Rat>
    vector& sub_scal_div(
            const vector<AltBasis, AltCoeffs, AltVecType, AltArgs...>& rhs,
            const Rat& scal
    );


    template <typename Vector>
    friend std::enable_if_t<std::is_base_of<vector, Vector>::value, Vector>
    operator-(const Vector& arg)
    {
        using coeffs = typename Vector::coefficient_ring;
        return arg>-unary_op([](const scalar_type& s) { return coeffs::uminus(s); });
    }

    template <typename Vector, typename Scal>
    friend std::enable_if_t<std::is_base_of<vector, Vector>::value, Vector>
    operator*(const Vector& arg, const Scal& scalar)
    {
        using coeffs = typename Vector::coefficient_ring;
        scalar_type multiplier(scalar);
        return arg->unary_op([multiplier](const scalar_type& s) { return coeffs::mul(s, multiplier); });
    };

    template <typename Vector, typename Scal>
    friend std::enable_if_t<std::is_base_of<vector, Vector>::value, Vector>
    operator*(const Scal& scalar, const Vector& arg)
    {
        using coeffs = typename Vector::coefficient_ring;
        scalar_type multiplier(scalar);
        return arg->unary_op([multiplier](const scalar_type& s) { return coeffs::mul(multiplier, s); });
    };

    template <typename Vector, typename Rat>
    friend std::enable_if_t<std::is_base_of<vector, Vector>::value, Vector>
    operator/(const Vector& arg, const Rat& scalar)
    {
        using coeffs = typename Vector::coefficient_ring;
        scalar_type multiplier(coefficient_ring::one()/scalar);
        return arg->unary_op([multiplier](const scalar_type& s) { return coeffs::mul(s, multiplier); });
    }

    template <typename LVector>
    friend std::enable_if_t<std::is_base_of<vector, LVector>::value, LVector>
    operator+(const LVector& lhs, const vector& rhs)
    {
        using coeffs = typename LVector::coefficient_ring;
        return lhs->binary_op(rhs, coeffs::add);
    }




    template <typename LVector,
            typename RCoefficients,
            template <typename, typename, typename...> class RVecType,
            typename... RArgs>
    friend std::enable_if_t<
        std::is_base_of<vector, LVector>::value,
        LVector
    >
    operator+(const LVector& lhs, const vector<Basis, RCoefficients, RVecType, RArgs...>& rhs)
    {
        using rscalar_type = typename coefficient_trait<RCoefficients>::scalar_type;
        return lhs->binary_op(rhs, [](const scalar_type& l, const rscalar_type& r) { return l + scalar_type(r); });
    }





};

namespace dtl {

template <typename B, typename C, template <typename, typename, typename...> class V, typename... A>
struct owned_vector_impl<vector<B, C, V, A...>>
{
    using type = owned_type_of<V<B, C, A...>>;
};

template <typename B, typename C, template <typename, typename, typename...> class V, typename... A>
struct view_vector_impl<vector<B, C, V, A...>>
{
    using type = view_type_of<V<B, C, A...>>;
};



} // namespace dtl





// Implementations of inplace fused methods


template<typename Basis,
        typename Coefficients,
        template <typename, typename, typename...> class VectorType,
        typename... Args>
template<typename Key, typename Scal>
vector<Basis, Coefficients, VectorType, Args...>&
vector<Basis, Coefficients, VectorType, Args...>::add_scal_prod(const Key& key, const Scal& scal)
{
    coefficient_ring::add_inplace((*p_impl)[key_type(key)], scalar_type(scal));
    return *this;
}
template<typename Basis,
        typename Coefficients,
        template <typename, typename, typename...> class VectorType,
        typename... Args>
template<typename Key, typename Rat>
vector<Basis, Coefficients, VectorType, Args...>&
vector<Basis, Coefficients, VectorType, Args...>::add_scal_div(const Key& key, const Rat& scal)
{
    coefficient_ring::add_inplace((*p_impl)[key_type(key)],
            coefficient_ring::div(coefficient_ring::one(), rational_type(scal)));
    return *this;
}
template<typename Basis,
        typename Coefficients,
        template<typename, typename, typename...> class VectorType,
        typename... Args>
template<typename Key, typename Scal>
vector<Basis, Coefficients, VectorType, Args...>&
vector<Basis, Coefficients, VectorType, Args...>::sub_scal_prod(const Key& key, const Scal& scal)
{
    coefficient_ring::sub_inplace((*p_impl)[key_type(key)], scalar_type(scal));
    return *this;
}
template<typename Basis,
        typename Coefficients,
        template<typename, typename, typename...> class VectorType,
        typename... Args>
template<typename Key, typename Rat>
vector<Basis, Coefficients, VectorType, Args...>&
vector<Basis, Coefficients, VectorType, Args...>::sub_scal_div(const Key& key, const Rat& scal)
{
    coefficient_ring::sub_inplace((*p_impl)[key_type(key)],
            coefficient_ring::div(coefficient_ring::one(), rational_type(scal)));
    return *this;
}
template<typename Basis,
        typename Coefficients,
        template<typename, typename, typename...> class VectorType,
        typename... Args>
template<typename Scal>
vector<Basis, Coefficients, VectorType, Args...>&
vector<Basis, Coefficients, VectorType, Args...>::add_scal_prod(const vector& rhs, const Scal& scal)
{
    return *this;
}
template<typename Basis,
        typename Coefficients,
        template<typename, typename, typename...> class VectorType,
        typename... Args>
template<typename Rat>
vector<Basis, Coefficients, VectorType, Args...>&
vector<Basis, Coefficients, VectorType, Args...>::add_scal_div(const vector& rhs, const Rat& scal)
{
    return *this;
}
template<typename Basis,
        typename Coefficients,
        template<typename, typename, typename...> class VectorType,
        typename... Args>
template<typename Scal>
vector<Basis, Coefficients, VectorType, Args...>&
vector<Basis, Coefficients, VectorType, Args...>::sub_scal_prod(const vector& rhs, const Scal& scal)
{
    return *this;
}
template<typename Basis,
        typename Coefficients,
        template<typename, typename, typename...> class VectorType,
        typename... Args>
template<typename Rat>
vector<Basis, Coefficients, VectorType, Args...>&
vector<Basis, Coefficients, VectorType, Args...>::sub_scal_div(const vector& rhs, const Rat& scal)
{
    return *this;
}
template<typename Basis,
        typename Coefficients,
        template<typename, typename, typename...> class VectorType,
        typename... Args>
template<template <typename, typename, typename...> class AltVecType,
        typename... AltArgs,
        typename Scal>
vector<Basis, Coefficients, VectorType, Args...>&
vector<Basis, Coefficients, VectorType, Args...>::add_scal_prod(
        const vector<Basis, Coefficients, AltVecType, AltArgs...>& rhs, const Scal& scal)
{
    return *this;
}
template<typename Basis,
        typename Coefficients,
        template<typename, typename, typename...> class VectorType,
        typename... Args>
template<template <typename, typename, typename...> class AltVecType,
        typename... AltArgs,
        typename Rat>
vector<Basis, Coefficients, VectorType, Args...>&
vector<Basis, Coefficients, VectorType, Args...>::add_scal_div(
        const vector<Basis, Coefficients, AltVecType, AltArgs...>& rhs, const Rat& scal)
{
    return *this;
}
template<typename Basis,
        typename Coefficients,
        template<typename, typename, typename...> class VectorType,
        typename... Args>
template<template <typename, typename, typename...> class AltVecType,
        typename... AltArgs,
        typename Scal>
vector<Basis, Coefficients, VectorType, Args...>&
vector<Basis, Coefficients, VectorType, Args...>::sub_scal_prod(
        const vector<Basis, Coefficients, AltVecType, AltArgs...>& rhs, const Scal& scal)
{
    return *this;
}
template<typename Basis,
        typename Coefficients,
        template<typename, typename, typename...> class VectorType,
        typename... Args>
template<template <typename, typename, typename...> class AltVecType,
        typename... AltArgs,
        typename Rat>
vector<Basis, Coefficients, VectorType, Args...>&
vector<Basis, Coefficients, VectorType, Args...>::sub_scal_div(
        const vector<Basis, Coefficients, AltVecType, AltArgs...>& rhs, const Rat& scal)
{
    return *this;
}
template<typename Basis,
        typename Coefficients,
        template<typename, typename, typename...> class VectorType,
        typename... Args>
template<typename AltBasis,
        typename AltCoeffs,
        template <typename, typename, typename...> class AltVecType,
        typename... AltArgs,
        typename Scal>
vector<Basis, Coefficients, VectorType, Args...>&
vector<Basis, Coefficients, VectorType, Args...>::add_scal_prod(
        const vector<AltBasis, AltCoeffs, AltVecType, AltArgs...>& rhs, const Scal& scal)
{
    return *this;
}
template<typename Basis,
        typename Coefficients,
        template<typename, typename, typename...> class VectorType,
        typename... Args>
template<typename AltBasis,
        typename AltCoeffs,
        template <typename, typename, typename...> class AltVecType,
        typename... AltArgs,
        typename Rat>
vector<Basis, Coefficients, VectorType, Args...>&
vector<Basis, Coefficients, VectorType, Args...>::add_scal_div(
        const vector<AltBasis, AltCoeffs, AltVecType, AltArgs...>& rhs, const Rat& scal)
{
    return *this;
}
template<typename Basis,
        typename Coefficients,
        template<typename, typename, typename...> class VectorType,
        typename... Args>
template<typename AltBasis,
        typename AltCoeffs,
        template <typename, typename, typename...> class AltVecType,
        typename... AltArgs,
        typename Scal>
vector<Basis, Coefficients, VectorType, Args...>&
vector<Basis, Coefficients, VectorType, Args...>::sub_scal_prod(
        const vector<AltBasis, AltCoeffs, AltVecType, AltArgs...>& rhs, const Scal& scal)
{
    return *this;
}
template<typename Basis,
        typename Coefficients,
        template<typename, typename, typename...> class VectorType,
        typename ... Args>
template<typename AltBasis,
        typename AltCoeffs,
        template <typename, typename, typename...> class AltVecType,
        typename... AltArgs,
        typename Rat>
vector<Basis, Coefficients, VectorType, Args...>&
vector<Basis, Coefficients, VectorType, Args...>::sub_scal_div(
        const vector<AltBasis, AltCoeffs, AltVecType, AltArgs...>& rhs, const Rat& scal)
{
    return *this;
}

} // namespace alg


#endif //LIBALGEBRA_LITE_VECTOR_H
