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
    std::shared_ptr<const Basis> p_basis;
    std::shared_ptr<vector_type> p_impl;
    friend struct dtl::vector_base_access;

protected:
    using basis_traits = basis_trait<Basis>;
    using coeff_traits = coefficient_trait<Coefficients>;

public:
    using basis_type        = Basis;
    using basis_pointer     = std::shared_ptr<const Basis>;
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

    vector(basis_pointer basis) : p_basis(std::move(basis)), p_impl(nullptr)
    {}

    vector(basis_pointer basis, std::initializer_list<scalar_type> args)
        : p_basis(std::move(basis)), p_impl(std::make_shared<vector_type>(p_basis.get(), args))
    {}

    vector(basis_pointer basis, const vector_type& arg)
        : p_basis(basis), p_impl(std::make_shared<vector_type>(arg))
    {
        assert(p_basis.get() == arg.p_basis);
    }


    LAL_INLINE_ALWAYS
    dimn_t size() const noexcept
    {
        if (p_impl) {
            return p_impl->size();
        }
        return 0;
    }
    LAL_INLINE_ALWAYS
    dimn_t dimension() const noexcept
    {
        if (p_impl) {
            return p_impl->dimension();
        }
        return 0;
    }
    LAL_INLINE_ALWAYS
    bool empty() const noexcept
    {
        if (p_impl) {
            return p_impl->empty();
        }
        return true;
    }


    LAL_INLINE_ALWAYS
    basis_pointer basis() const noexcept
    {
        return p_basis;
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
    {
        if (p_impl) {
            p_impl->clear();
        }
    }
    LAL_INLINE_ALWAYS
    iterator begin() noexcept
    {
        if (p_impl) {
            return p_impl->begin();
        }
        return iterator();
    }
    LAL_INLINE_ALWAYS
    iterator end() noexcept
    {
        if (p_impl) {
            return p_impl->end();
        }
        return iterator();
    }
    LAL_INLINE_ALWAYS
    const_iterator begin() const noexcept
    {
        if (p_impl) {
            return p_impl->begin();
        }
        return const_iterator();
    }
    LAL_INLINE_ALWAYS
    const_iterator end() const noexcept
    {
        if (p_impl) {
            return p_impl->end();
        }
        return const_iterator();
    }

    LAL_INLINE_ALWAYS
    vector_type& base_vector()
    {
        if (!p_impl) {
            p_impl = std::make_shared<vector_type>(p_basis.get());
        }
        return *p_impl;
    }
    LAL_INLINE_ALWAYS
    const vector_type& base_vector() const noexcept
    { return *p_impl; }

    vector clone() const
    {
        if (p_impl) {
            return vector(p_basis, *p_impl);
        }
        return vector(p_basis);
    }

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
        if (arg.p_impl) {
            return arg->unary_op([](const scalar_type& s) { return coeffs::uminus(s); });
        }
        return Vector(arg.p_basis);
    }

    template <typename Vector, typename Scal>
    friend std::enable_if_t<std::is_base_of<vector, Vector>::value, Vector>
    operator*(const Vector& arg, const Scal& scalar)
    {
        using coeffs = typename Vector::coefficient_ring;
        if (arg.p_impl) {
            scalar_type multiplier(scalar);
            return arg->unary_op([multiplier](const scalar_type& s) { return coeffs::mul(s, multiplier); });
        }
        return Vector(arg.p_basis);
    };

    template <typename Vector, typename Scal>
    friend std::enable_if_t<std::is_base_of<vector, Vector>::value, Vector>
    operator*(const Scal& scalar, const Vector& arg)
    {
        using coeffs = typename Vector::coefficient_ring;
        if (arg.p_impl) {
            scalar_type multiplier(scalar);
            return arg->unary_op([multiplier](const scalar_type& s) { return coeffs::mul(multiplier, s); });
        }
        return Vector(arg.p_basis);
    };

    template <typename Vector, typename Rat>
    friend std::enable_if_t<std::is_base_of<vector, Vector>::value, Vector>
    operator/(const Vector& arg, const Rat& scalar)
    {
        using coeffs = typename Vector::coefficient_ring;
        if (arg.p_impl) {
            scalar_type multiplier(coefficient_ring::one()/scalar);
            return arg->unary_op([multiplier](const scalar_type& s) { return coeffs::mul(s, multiplier); });
        }
        return Vector(arg.p_basis);
    }

    template <typename LVector>
    friend std::enable_if_t<std::is_base_of<vector, LVector>::value, LVector>
    operator+(const LVector& lhs, const vector& rhs)
    {
        using coeffs = typename LVector::coefficient_ring;
        if (!rhs.p_basis) {
            return lhs.clone();
        }
        if (!lhs.p_basis) {
            return LVector(rhs.p_basis);
        }
        if (lhs.p_impl) {
            return lhs.p_impl->binary_op(rhs, coeffs::add);
        }
        return LVector(rhs.clone());
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
        if (!rhs.p_basis) {
            return lhs.clone();
        }
        if (!lhs.p_basis) {
            return LVector(vector(rhs.p_basis));
        }
        if (lhs.p_impl) {
            return lhs.p_impl->binary_op(rhs, [](const scalar_type& l, const rscalar_type& r) { return l + scalar_type(r); });
        }
        return LVector(lhs.p_basis);

    }

    template <typename LVector, typename Scal>
    friend std::enable_if<std::is_base_of<vector, LVector>::value, LVector&>
    operator*=(LVector& lhs, const Scal& scal)
    {
        if (lhs.p_basis && lhs.p_impl) {
            scalar_type multiplier(scal);
            lhs.p_impl->unary_op_inplace([multiplier](auto& value) {
                return coefficient_ring::mul_inplace(value, multiplier);
            });
        }
        return lhs;
    }

    template <typename LVector, typename Rat>
    friend std::enable_if<std::is_base_of<vector, LVector>::value, LVector&>
    operator/=(LVector& lhs, const Rat& scal)
    {
        if (lhs.p_basis && lhs.p_impl) {
            auto multiplier = coefficient_ring::div(coefficient_ring::one(), rational_type(scal));
            lhs.p_impl->unary_op_inplace([multiplier](auto& value) {
                return coefficient_ring::mul_inplace(value, multiplier);
            });
        }
        return lhs;
    }



    template <typename LVector>
    friend std::enable_if<std::is_base_of<vector, LVector>::value, LVector&>
    operator+=(LVector& lhs, const vector& rhs)
    {
        if (!rhs.p_basis || !rhs.p_impl) {
            return lhs;
        }
        if (!lhs.p_basis || !lhs.p_impl) {
            lhs.vector::operator=(rhs);
            return lhs;
        }
        return lhs.p_impl->binary_op_inplace(rhs, LVector::coefficient_ring::add_inplace);
    }

    template <typename LVector>
    friend std::enable_if<std::is_base_of<vector, LVector>::value, LVector&>
    operator-=(LVector& lhs, const vector& rhs)
    {
        if (!rhs.p_basis || !rhs.p_impl) {
            return lhs;
        }
        if (!lhs.p_basis || !lhs.p_impl) {
            lhs.vector::operator=(-rhs);
            return lhs;
        }
        return lhs.p_impl->binary_op_inplace(rhs, LVector::coefficient_ring::sub_inplace);
    }





};






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
