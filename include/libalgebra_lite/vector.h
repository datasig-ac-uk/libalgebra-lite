//
// Created by user on 08/08/22.
//

#ifndef LIBALGEBRA_LITE_VECTOR_H
#define LIBALGEBRA_LITE_VECTOR_H

#include "implementation_types.h"

#include <memory>

#include "vector_traits.h"
#include "coefficients.h"
#include "basis_traits.h"
#include "registry.h"


namespace lal {

namespace dtl {


template <typename VectorType>
struct storage_base
{
    using vector_type = VectorType;
    using vect_traits = vector_traits<VectorType>;

    using basis_type = typename vect_traits::basis_type;
    using coefficient_ring = typename vect_traits::coefficient_ring;

    using basis_traits      = basis_trait<basis_type>;
    using coeff_traits      = coefficient_trait<coefficient_ring>;

    using basis_pointer     = std::shared_ptr<const basis_type>;
    using key_type          = typename basis_traits::key_type;
    using scalar_type       = typename coeff_traits::scalar_type;
    using rational_type     = typename coeff_traits::rational_type;

    using iterator          = typename vector_type::iterator;
    using const_iterator    = typename vector_type::const_iterator;
    using reference         = typename vector_type::reference;
    using const_reference   = typename vector_type::const_reference;


    basis_pointer p_basis;
    std::shared_ptr<vector_type> p_impl;

    storage_base(storage_base&& other) noexcept
        : p_basis(std::move(other.p_basis)), p_impl(std::move(other.p_impl))
    {}

    explicit storage_base(basis_pointer basis, VectorType&& arg)
        : p_basis(std::move(basis)), p_impl(std::make_shared<VectorType>(std::move(arg)))
    {}

    template <typename... Args>
    explicit storage_base(key_type k, scalar_type s, Args&&... args)
        : p_basis(basis_registry<basis_type>::get(std::forward<Args>(args)...)),
          p_impl(p_basis.get(), k, s)
    {}

    template <typename... BasisArgs>
    explicit storage_base(VectorType&& arg, BasisArgs&&... basis_args)
        : p_basis(basis_traits::basis_factory(std::forward<BasisArgs>(basis_args)...)),
          p_impl(std::make_shared<VectorType>(std::move(arg)))
    {}

    template <typename... VectorArgs>
    explicit storage_base(basis_pointer basis, VectorArgs&&... args)
        : p_basis(std::move(basis)),
          p_impl(std::make_shared<VectorType>(p_basis.get(), std::forward<VectorArgs>(args)...))
    {}


    dimn_t size() const noexcept
    {
        if (p_impl) {
            return p_impl->size();
        }
        return 0;
    }

    dimn_t dimension() const noexcept
    {
        if (p_impl) {
            return p_impl->dimension();
        }
        return 0;
    }

    bool empty() const noexcept
    {
        if (p_impl) {
            return p_impl->empty();
        }
        return true;
    }

    template <typename KeyType>
    const_reference operator[](const KeyType& key) const
    {
        if (p_impl) {
            return (*p_impl)[key];
        }
        return coefficient_ring::zero();
    }

    const basis_type& basis() const noexcept
    {
        assert(static_cast<bool>(p_basis));
        return *p_basis;
    }

    basis_pointer get_basis() const noexcept
    {
        return p_basis;
    }

    const_iterator begin() const noexcept
    {
        if (p_impl) {
            return p_impl->cbegin();
        }
        return const_iterator();
    }

    const_iterator end() const noexcept
    {
        if (p_impl) {
            return p_impl->cend();
        }
        return const_iterator();
    }

    const vector_type& base_vector() const noexcept
    { return *p_impl; }
};


template <typename VectorType>
struct standard_storage : public storage_base<VectorType>
{
    using base = storage_base<VectorType>;

    using typename base::vector_type;

    using typename base::vect_traits;
    using typename base::basis_traits;
    using typename base::coeff_traits;

    using typename base::basis_type;
    using typename base::key_type;
    using typename base::coefficient_ring;
    using typename base::scalar_type;
    using typename base::rational_type;
    using typename base::iterator;
    using typename base::const_iterator;
    using typename base::reference;
    using typename base::const_reference;

    using typename base::basis_pointer;

    using base::p_impl;
    using base::p_basis;

    // Inherit the constructors from the storage_base
    using base::base;


    using base::operator[];
    using base::begin;
    using base::end;
    using base::base_vector;

    template <typename KeyType>
    reference operator[](const KeyType& key)
    {
        ensure_created();
        return (*p_impl)[key];
    }

    void clear()
    {
        if (p_impl) {
            p_impl->clear();
        }
    }

    iterator begin() noexcept
    {
        if (p_impl) {
            return p_impl->begin();
        }
        return iterator();
    }

    iterator end() noexcept
    {
        if (p_impl) {
            return p_impl->end();
        }
        return iterator();
    }

    vector_type& base_vector()
    {
        ensure_created();
        return *p_impl;
    }
    
    void ensure_created()
    {
        if (!p_impl) {
            p_impl = std::make_shared<vector_type>(p_basis.get());
        }
    }


};

template <typename VectorType>
class copy_on_write_storage : public storage_base<VectorType>
{
    using base = storage_base<VectorType>;
    enum {
        BORROWED,
        OWNED
    } m_state;

    void convert_to_copy()
    {
        if (m_state == BORROWED) {
            base::p_impl = std::make_shared<VectorType>(*base::p_impl);
        }
    }

public:

    using typename base::vector_type;

    using typename base::vect_traits;
    using typename base::basis_traits;
    using typename base::coeff_traits;

    using typename base::basis_type;
    using typename base::key_type;
    using typename base::coefficient_ring;
    using typename base::scalar_type;
    using typename base::rational_type;
    using typename base::iterator;
    using typename base::const_iterator;
    using typename base::reference;
    using typename base::const_reference;

    using typename base::basis_pointer;

    using base::p_impl;
    using base::p_basis;

    using base::base;





};



} // namespace dtl


template <typename Basis,
        typename Coefficients,
        template <typename, typename> class VectorType,
        template <typename> class StorageModel=dtl::standard_storage>
class vector : protected StorageModel<VectorType<Basis, Coefficients>>
{
    using base_type = StorageModel<VectorType<Basis, Coefficients>>;
    using registry = basis_registry<Basis>;

public:
    using typename base_type::vector_type;
    using typename base_type::basis_type;
    using typename base_type::basis_pointer;
    using typename base_type::key_type;
    using typename base_type::coefficient_ring;
    using typename base_type::scalar_type;
    using typename base_type::rational_type;

    using typename base_type::iterator;
    using typename base_type::const_iterator;
    using typename base_type::reference;
    using typename base_type::const_reference;

protected:

    using base_type::p_basis;
    using base_type::p_impl;

    vector(basis_pointer basis, vector_type&& arg)
        : base_type(basis, std::make_shared<vector_type>(std::move(arg)))
    {}

public:

    vector() : base_type(registry::get())
    {}

    vector(const vector& other) : base_type(other.p_basis)
    {
        if (other.p_impl) {
            p_impl = std::make_shared<vector_type>(*other.p_impl);
        }
    }

    vector(vector&& other) noexcept : base_type(std::move(other.p_basis))
    {
        p_impl = std::move(other.p_impl);
    }

    template <typename Key, typename Scalar>
    explicit vector(Key k, Scalar s)
        : base_type(registry::get(), key_type(k), scalar_type(s))
    {}

    template <typename Key, typename Scalar>
    explicit vector(basis_pointer basis, Key key, Scalar s)
        : base_type(std::move(basis), key_type(key), scalar_type(s))
    {}

    explicit vector(basis_pointer basis) : base_type(std::move(basis))
    {}

    vector(basis_pointer basis, std::initializer_list<scalar_type> args)
        : base_type(basis, vector_type(basis.get(), args))
    {}

    vector(basis_pointer basis, const vector_type& arg)
        : base_type(basis, arg)
    {
        assert(p_basis.get() == arg.p_basis);
    }
    vector clone() const
    {
        if (p_impl) {
            return vector(p_basis, *p_impl);
        }
        return vector(p_basis);
    }

    using base_type::basis;
    using base_type::size;
    using base_type::dimension;
    using base_type::empty;
    using base_type::clear;
    using base_type::operator[];
    using base_type::begin;
    using base_type::end;
    using base_type::base_vector;
    using base_type::get_basis;


    vector& operator=(const vector& other)
    {
        if (p_basis.get() != other.p_basis.get()) {
            p_basis = other.p_basis;
        }
        if (other.p_impl){
            p_impl = std::make_shared<vector_type>(*other.p_impl);
        } else {
            p_impl = nullptr;
        }
        return *this;
    }

    vector& operator=(vector&& other) noexcept
    {
        p_basis = std::move(other.p_basis);
        p_impl = std::move(other.p_impl);
        return *this;
    }

private:

    template <typename, typename, typename, typename=void>
    struct has_iterator_inplace_binop : std::false_type {};

    template <typename V, typename F, typename I>
    struct has_iterator_inplace_binop<V, F, I, std::void_t<
            decltype(V::template inplace_binop(
                    std::declval<I>(),
                    std::declval<I>(),
                    std::declval<F>()
                    ))>>
        : std::true_type
    {};


public:

    template <typename Iter, typename C>
    std::enable_if_t<
            has_iterator_inplace_binop<
                    vector_type,
                    decltype(coefficient_ring::template add_inplace<>),
                    Iter
                    >::value,
            vector&>
    add_inplace(Iter begin, Iter end)
    {
        base_type::ensure_created();
        return p_impl->inplace_binop(begin, end, coefficient_ring::add_inplace);
    }

    template <typename Iter>
    std::enable_if_t<
            has_iterator_inplace_binop<
                    vector_type,
                    decltype(coefficient_ring::template add_inplace<>),
                    Iter>::value,
            vector&>
    sub_inplace(Iter begin, Iter end)
    {
        base_type::ensure_created();
        return p_impl->inplace_binop(begin, end, coefficient_ring::sub_inplace);
    }


    template <typename Iter>
    std::enable_if_t<
            !has_iterator_inplace_binop<
                    vector_type,
                    decltype(coefficient_ring::template add_inplace<>),
                    Iter>::value,
            vector&>
    add_inplace(Iter begin, Iter end)
    {
        base_type::ensure_created();
        for (auto it = begin; it !=end; ++it) {
            coefficient_ring::add_inplace((*p_impl)[it->first], it->second);
        }
        return *this;
    }



    template <typename Iter>
    std::enable_if_t<
            !has_iterator_inplace_binop<
                    vector_type,
                    decltype(coefficient_ring::template sub_inplace<>),
                    Iter>::value,
            vector&>
    sub_inplace(Iter begin, Iter end)
    {
        base_type::ensure_created();
        for (auto it = begin; it != end; ++it) {
            coefficient_ring::sub_inplace(*p_impl, it->second);
        }
    }


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

    template <template <typename, typename> class AltVecType,
            template <typename> class AltStorageModel,
            typename Scal>
    vector& add_scal_prod(
            const vector<Basis, Coefficients, AltVecType, AltStorageModel>& rhs,
            const Scal& scal
    );
    template <template <typename, typename> class AltVecType,
            template <typename> class AltStorageModel,
            typename Rat>
    vector& add_scal_div(
            const vector<Basis, Coefficients, AltVecType, AltStorageModel>& rhs,
            const Rat& scal
    );
    template <template <typename, typename> class AltVecType,
            template <typename> class AltStorageModel,
            typename Scal>
    vector& sub_scal_prod(
            const vector<Basis, Coefficients, AltVecType, AltStorageModel>& rhs,
            const Scal& scal
    );
    template <template <typename, typename> class AltVecType,
            template <typename> class AltStorageModel,
            typename Rat>
    vector& sub_scal_div(
            const vector<Basis, Coefficients, AltVecType, AltStorageModel>& rhs,
            const Rat& scal
    );

    template <
            typename AltBasis,
            typename AltCoeffs,
            template <typename, typename> class AltVecType,
            template <typename> class AltStorageModel,
            typename Scal>
    vector& add_scal_prod(
            const vector<AltBasis, AltCoeffs, AltVecType, AltStorageModel>& rhs,
            const Scal& scal
    );
    template <
            typename AltBasis,
            typename AltCoeffs,
            template <typename, typename> class AltVecType,
            template <typename> class AltStorageModel,
            typename Rat>
    vector& add_scal_div(
            const vector<AltBasis, AltCoeffs, AltVecType, AltStorageModel>& rhs,
            const Rat& scal
    );
        template <
            typename AltBasis,
            typename AltCoeffs,
            template <typename, typename> class AltVecType,
                template <typename> class AltStorageModel,
            typename Scal>
    vector& sub_scal_prod(
            const vector<AltBasis, AltCoeffs, AltVecType, AltStorageModel>& rhs,
            const Scal& scal
    );
    template <
            typename AltBasis,
            typename AltCoeffs,
            template <typename, typename> class AltVecType,
            template <typename> class AltStorageModel,
            typename Rat>
    vector& sub_scal_div(
            const vector<AltBasis, AltCoeffs, AltVecType, AltStorageModel>& rhs,
            const Rat& scal
    );


    template <typename Vector>
    friend std::enable_if_t<std::is_base_of<vector, Vector>::value, Vector>
    operator-(const Vector& arg)
    {
        using coeffs = typename Vector::coefficient_ring;
        vector result_inner(arg.p_basis);
        if (arg.p_impl) {
            result_inner.p_impl = std::make_shared<vector_type>(
                    arg.p_impl->unary_op([](const scalar_type& s) { return coeffs::uminus(s); }));
        }
        return Vector(std::move(result_inner));
    }

    template <typename Vector, typename Scal>
    friend std::enable_if_t<std::is_base_of<vector, Vector>::value, Vector>
    operator*(const Vector& arg, const Scal& scalar)
    {
        using coeffs = typename Vector::coefficient_ring;
        if (arg.p_impl) {
            scalar_type multiplier(scalar);
            return Vector(arg.p_basis,
                    arg.p_impl->unary_op([multiplier](const scalar_type& s) { return coeffs::mul(s, multiplier); }));
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
            return Vector(arg.p_basis,
                    arg.p_impl->unary_op([multiplier](const scalar_type& s) { return coeffs::mul(multiplier, s); }));
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
            return Vector(arg.p_basis,
                    arg.p_impl->unary_op([multiplier](const scalar_type& s) { return coeffs::mul(s, multiplier); }));
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
            template <typename, typename> class RVecType,
            template <typename> class RStorageModel>
    friend std::enable_if_t<
        std::is_base_of<vector, LVector>::value,
        LVector
    >
    operator+(const LVector& lhs, const vector<Basis, RCoefficients, RVecType, RStorageModel>& rhs)
    {
        using rscalar_type = typename coefficient_trait<RCoefficients>::scalar_type;
        if (!rhs.p_basis) {
            return lhs.clone();
        }
        if (!lhs.p_basis) {
            return LVector(vector(rhs.p_basis));
        }
        if (lhs.p_impl) {
            return LVector(lhs.p_basis,
                    lhs.p_impl->binary_op(rhs, [](const scalar_type& l, const rscalar_type& r) { return l + scalar_type(r); }));
        }
        return LVector(lhs.p_basis);

    }

    template <typename LVector, typename Scal>
    friend std::enable_if_t<std::is_base_of<vector, LVector>::value, LVector&>
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
    friend std::enable_if_t<std::is_base_of<vector, LVector>::value, LVector&>
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
    friend std::enable_if_t<std::is_base_of<vector, LVector>::value, LVector&>
    operator+=(LVector& lhs, const vector& rhs)
    {
        if (!rhs.p_basis || !rhs.p_impl) {
            return lhs;
        }
        if (!lhs.p_basis) {
            lhs.p_basis = rhs.p_basis;
        }
        lhs.ensure_created();
        lhs.p_impl->inplace_binary_op(rhs.base_vector(), [](auto& l, const auto& r) { l += r; });
        return lhs;
    }

    template <typename LVector>
    friend std::enable_if_t<std::is_base_of<vector, LVector>::value, LVector&>
    operator-=(LVector& lhs, const vector& rhs)
    {
        if (!rhs.p_basis || !rhs.p_impl) {
            return lhs;
        }
        lhs.ensure_created();
        lhs.p_impl->inplace_binary_op(rhs.base_vector(), LVector::coefficient_ring::template sub_inplace<>);
        return lhs;
    }



    template <typename Vector>
    static Vector new_like(const Vector& arg)
    {
        return Vector(arg.p_basis);
    }


};

template <typename B, typename C, template <typename, typename> class VT, template <typename> class SM>
bool operator==(const vector<B, C, VT, SM>& lhs, const vector<B, C, VT, SM>& rhs) noexcept
{
    return false;
}
template <typename B, typename C, template <typename, typename> class VT, template <typename> class SM>
bool operator!=(const vector<B, C, VT, SM>& lhs, const vector<B, C, VT, SM>& rhs) noexcept
{
    return !(lhs == rhs);
}




template <typename Basis, typename Coeff, template <typename, typename> class VectorType, template <typename> class StorageModel>
std::ostream& operator<<(std::ostream& os, const vector<Basis, Coeff, VectorType, StorageModel>& vect)
{
    const auto& basis = vect.basis();
    const auto& zero = coefficient_trait<Coeff>::coefficient_ring::zero();
    os << "{ ";
    for (auto item : vect) {
        auto val = item.value();
        if (item.value() != zero) {
            os << item.value() << '(';
            basis.print_key(os, item.key());
            os << ") ";
        }
    }
    os <<'}';
    return os;
}


// Implementations of inplace fused methods


template<typename Basis,
        typename Coefficients,
        template <typename, typename> class VectorType,
        template <typename> class StorageModel>
template<typename Key, typename Scal>
vector<Basis, Coefficients, VectorType, StorageModel>&
vector<Basis, Coefficients, VectorType, StorageModel>::add_scal_prod(const Key& key, const Scal& scal)
{
    base_type::ensure_created();
    coefficient_ring::add_inplace((*p_impl)[key_type(key)], scalar_type(scal));
    return *this;
}
template<typename Basis,
        typename Coefficients,
        template <typename, typename> class VectorType,
        template <typename> class StorageModel>
template<typename Key, typename Rat>
vector<Basis, Coefficients, VectorType, StorageModel>&
vector<Basis, Coefficients, VectorType, StorageModel>::add_scal_div(const Key& key, const Rat& scal)
{
    base_type::ensure_created();
    coefficient_ring::add_inplace((*p_impl)[key_type(key)],
            coefficient_ring::div(coefficient_ring::one(), rational_type(scal)));
    return *this;
}
template<typename Basis,
        typename Coefficients,
        template<typename, typename> class VectorType,
        template <typename> class StorageModel>
template<typename Key, typename Scal>
vector<Basis, Coefficients, VectorType, StorageModel>&
vector<Basis, Coefficients, VectorType, StorageModel>::sub_scal_prod(const Key& key, const Scal& scal)
{
    base_type::ensure_created();
    coefficient_ring::sub_inplace((*p_impl)[key_type(key)], scalar_type(scal));
    return *this;
}
template<typename Basis,
        typename Coefficients,
        template<typename, typename> class VectorType,
        template <typename> class StorageModel>
template<typename Key, typename Rat>
vector<Basis, Coefficients, VectorType, StorageModel>&
vector<Basis, Coefficients, VectorType, StorageModel>::sub_scal_div(const Key& key, const Rat& scal)
{
    base_type::ensure_created();
    coefficient_ring::sub_inplace((*p_impl)[key_type(key)],
            coefficient_ring::div(coefficient_ring::one(), rational_type(scal)));
    return *this;
}
template<typename Basis,
        typename Coefficients,
        template<typename, typename> class VectorType,
        template <typename> class StorageModel>
template<typename Scal>
vector<Basis, Coefficients, VectorType, StorageModel>&
vector<Basis, Coefficients, VectorType, StorageModel>::add_scal_prod(const vector& rhs, const Scal& scal)
{
    base_type::ensure_created();
    return *this;
}
template<typename Basis,
        typename Coefficients,
        template<typename, typename> class VectorType,
        template <typename> class StorageModel>
template<typename Rat>
vector<Basis, Coefficients, VectorType, StorageModel>&
vector<Basis, Coefficients, VectorType, StorageModel>::add_scal_div(const vector& rhs, const Rat& scal)
{
    base_type::ensure_created();
    return *this;
}
template<typename Basis,
        typename Coefficients,
        template<typename, typename> class VectorType,
        template <typename> class StorageModel>
template<typename Scal>
vector<Basis, Coefficients, VectorType, StorageModel>&
vector<Basis, Coefficients, VectorType, StorageModel>::sub_scal_prod(const vector& rhs, const Scal& scal)
{
    base_type::ensure_created();
    return *this;
}
template<typename Basis,
        typename Coefficients,
        template<typename, typename> class VectorType,
        template <typename> class StorageModel>
template<typename Rat>
vector<Basis, Coefficients, VectorType, StorageModel>&
vector<Basis, Coefficients, VectorType, StorageModel>::sub_scal_div(const vector& rhs, const Rat& scal)
{
    base_type::ensure_created();
    return *this;
}
template<typename Basis,
        typename Coefficients,
        template<typename, typename> class VectorType,
        template <typename> class StorageModel>
template<template <typename, typename> class AltVecType,
        template <typename> class AltStorageModel,
        typename Scal>
vector<Basis, Coefficients, VectorType, StorageModel>&
vector<Basis, Coefficients, VectorType, StorageModel>::add_scal_prod(
        const vector<Basis, Coefficients, AltVecType,  AltStorageModel>& rhs, const Scal& scal)
{
    base_type::ensure_created();
    return *this;
}
template<typename Basis,
        typename Coefficients,
        template<typename, typename> class VectorType,
        template <typename> class StorageModel>
template<template <typename, typename> class AltVecType,
        template <typename> class AltStorageModel,
        typename Rat>
vector<Basis, Coefficients, VectorType, StorageModel>&
vector<Basis, Coefficients, VectorType, StorageModel>::add_scal_div(
        const vector<Basis, Coefficients, AltVecType, AltStorageModel>& rhs, const Rat& scal)
{
    base_type::ensure_created();
    return *this;
}
template<typename Basis,
        typename Coefficients,
        template<typename, typename> class VectorType,
        template <typename> class StorageModel>
template<template <typename, typename> class AltVecType,
        template <typename> class AltStorageModel,
        typename Scal>
vector<Basis, Coefficients, VectorType, StorageModel>&
vector<Basis, Coefficients, VectorType, StorageModel>::sub_scal_prod(
        const vector<Basis, Coefficients, AltVecType, AltStorageModel>& rhs, const Scal& scal)
{
    base_type::ensure_created();
    return *this;
}
template<typename Basis,
        typename Coefficients,
        template<typename, typename> class VectorType,
        template <typename> class StorageModel>
template<template <typename, typename> class AltVecType,
        template <typename> class AltStorageModel,
        typename Rat>
vector<Basis, Coefficients, VectorType, StorageModel>&
vector<Basis, Coefficients, VectorType, StorageModel>::sub_scal_div(
        const vector<Basis, Coefficients, AltVecType, AltStorageModel>& rhs, const Rat& scal)
{
    base_type::ensure_created();
    return *this;
}
template<typename Basis,
        typename Coefficients,
        template<typename, typename> class VectorType,
        template <typename> class StorageModel>
template<typename AltBasis,
        typename AltCoeffs,
        template <typename, typename> class AltVecType,
        template <typename> class AltStorageModel,
        typename Scal>
vector<Basis, Coefficients, VectorType, StorageModel>&
vector<Basis, Coefficients, VectorType, StorageModel>::add_scal_prod(
        const vector<AltBasis, AltCoeffs, AltVecType, AltStorageModel>& rhs, const Scal& scal)
{
    base_type::ensure_created();
    return *this;
}
template<typename Basis,
        typename Coefficients,
        template<typename, typename> class VectorType,
        template <typename> class StorageModel>
template<typename AltBasis,
        typename AltCoeffs,
        template <typename, typename> class AltVecType,
        template <typename> class AltStorageModel,
        typename Rat>
vector<Basis, Coefficients, VectorType, StorageModel>&
vector<Basis, Coefficients, VectorType, StorageModel>::add_scal_div(
        const vector<AltBasis, AltCoeffs, AltVecType, AltStorageModel>& rhs, const Rat& scal)
{
    base_type::ensure_created();
    return *this;
}
template<typename Basis,
        typename Coefficients,
        template<typename, typename> class VectorType,
        template <typename> class StorageModel>
template<typename AltBasis,
        typename AltCoeffs,
        template <typename, typename> class AltVecType,
        template <typename> class AltStorageModel,
        typename Scal>
vector<Basis, Coefficients, VectorType, StorageModel>&
vector<Basis, Coefficients, VectorType, StorageModel>::sub_scal_prod(
        const vector<AltBasis, AltCoeffs, AltVecType, AltStorageModel>& rhs, const Scal& scal)
{
    base_type::ensure_created();
    return *this;
}
template<typename Basis,
        typename Coefficients,
        template<typename, typename> class VectorType,
        template <typename> class StorageModel>
template<typename AltBasis,
        typename AltCoeffs,
        template <typename, typename> class AltVecType,
        template <typename> class AltStorageModel,
        typename Rat>
vector<Basis, Coefficients, VectorType, StorageModel>&
vector<Basis, Coefficients, VectorType, StorageModel>::sub_scal_div(
        const vector<AltBasis, AltCoeffs, AltVecType, AltStorageModel>& rhs, const Rat& scal)
{
    base_type::ensure_created();
    return *this;
}

} // namespace alg


#endif //LIBALGEBRA_LITE_VECTOR_H
