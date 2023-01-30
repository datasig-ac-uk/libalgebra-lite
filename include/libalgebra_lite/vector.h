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

    using basis_type        = typename vect_traits::basis_type;
    using registry          = basis_registry<basis_type>;
    using coefficient_ring  = typename vect_traits::coefficient_ring;

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

    template <typename... Args>
    storage_base(Args... args) : p_basis(registry::get(std::forward<Args>(args)...))
    {}

    storage_base(basis_pointer&& basis) noexcept : p_basis(std::move(basis))
    {}

    storage_base(storage_base&& other) noexcept : p_basis(std::move(other.p_basis))
    {}

    const basis_type& basis() const noexcept
    {
        assert(static_cast<bool>(p_basis));
        return *p_basis;
    }

    basis_pointer get_basis() const noexcept
    {
        return p_basis;
    }

};


template <typename VectorType>
class standard_storage : public storage_base<VectorType>
{
public:
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

    using base::p_basis;

private:
    vector_type m_instance;

protected:
    const vector_type& instance() const noexcept { return m_instance; }
    vector_type& instance() noexcept { return m_instance; }

public:

    explicit standard_storage(basis_pointer&& basis) : base(std::move(basis)),
        m_instance(&base::basis())
    {}

    explicit standard_storage(basis_pointer basis, vector_type&& data)
        : base(std::move(basis)), m_instance(std::move(data))
    {
        assert(m_instance.basis() == &*p_basis);
    }

    template <typename... Args>
    explicit standard_storage(Args... args) : base(std::forward<Args>(args)...),
        m_instance(&base::basis())
    {}

    template <typename... VArgs>
    explicit standard_storage(basis_pointer basis) : base(std::move(basis)),
        , m_instance(&base::basis(), std::forward<VArgs>(args)...)
    {}

    vector_type& base_vector() noexcept { return instance(); }
    const vector_type& base_vector() const noexcept { return instance(); }



};



} // namespace dtl


template <typename Basis,
        typename Coefficients,
        template <typename, typename> class VectorType,
        template <typename> class StorageModel=dtl::standard_storage>
class vector : protected StorageModel<VectorType<Basis, Coefficients>>
{
    using base_type = StorageModel<VectorType<Basis, Coefficients>>;

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


    template <typename Key, typename Scalar>
    explicit vector(Key k, Scalar s) : base_type(p_basis, key_type(k), scalar_type(s))
    {
    }

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
    }

    using base_type::basis;
    using base_type::get_basis;

    dimn_t size() const noexcept
    {
        return base_type::instance().size();
    }

    dimn_t dimension() const noexcept
    {
        return base_type::instance().dimension();
    }

    bool empty() const noexcept
    {
        return base_type::instance().empty();
    }

    template<typename KeyType>
    const_reference operator[](const KeyType& key) const
    {
        return base_type::instance()[key_type(key)];
    }

    template<typename KeyType>
    reference operator[](const KeyType& key)
    {
        return base_type::instance()[key];
    }

    void clear()
    {
        base_type::instance().clear();
    }

    iterator begin() noexcept
    {
        return base_type::instance().begin();
    }

    iterator end() noexcept
    {
        return base_type::instance().end();
    }
    const_iterator begin() const noexcept
    {
        return base_type::instance().cbegin();
    }
    const_iterator end() const noexcept
    {
        return base_type::instance().cend();
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
        return base_type::instance().inplace_binop(begin, end, coefficient_ring::add_inplace);
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
        return base_type::instance().inplace_binop(begin, end, coefficient_ring::sub_inplace);
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
        const auto& self = base_type::instance();
        for (auto it = begin; it !=end; ++it) {
            self[it->first] += it->second;
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
        const auto& self = base_type::instance();
        for (auto it = begin; it != end; ++it) {
            self[it->first] -= it->second;
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
        // TODO: implement
        return Vector(std::move(result_inner));
    }

    template <typename Vector, typename Scal>
    friend std::enable_if_t<std::is_base_of<vector, Vector>::value, Vector>
    operator*(const Vector& arg, const Scal& scalar)
    {
        using coeffs = typename Vector::coefficient_ring;
        // TODO: implement

        return Vector(arg.p_basis);
    };

    template <typename Vector, typename Scal>
    friend std::enable_if_t<std::is_base_of<vector, Vector>::value, Vector>
    operator*(const Scal& scalar, const Vector& arg)
    {
        using coeffs = typename Vector::coefficient_ring;
        scalar_type multiplier(scalar);
            // TODO: implement
        return Vector(arg.p_basis);
    };

    template <typename Vector, typename Rat>
    friend std::enable_if_t<std::is_base_of<vector, Vector>::value, Vector>
    operator/(const Vector& arg, const Rat& scalar)
    {
        using coeffs = typename Vector::coefficient_ring;
        // TODO: implement
        return Vector(arg.p_basis);
    }

    template <typename LVector>
    friend std::enable_if_t<std::is_base_of<vector, LVector>::value, LVector>
    operator+(const LVector& lhs, const vector& rhs)
    {
        using coeffs = typename LVector::coefficient_ring;
        // TODO: implement
        return LVector(rhs.clone());
    }

    template <typename LVector>
    friend std::enable_if_t<std::is_base_of<vector, LVector>::value, LVector>
    operator-(const LVector& lhs, const vector& rhs)
    {
        using coeffs = typename LVector::coefficient_ring;
        // TODO: implmenet
        return LVector(-rhs);
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
        // TODO: implement
        return LVector(lhs.p_basis);

    }

    template <typename LVector, typename Scal>
    friend std::enable_if_t<std::is_base_of<vector, LVector>::value, LVector&>
    operator*=(LVector& lhs, const Scal& scal)
    {
        // TODO: implement
        return lhs;
    }

    template <typename LVector, typename Rat>
    friend std::enable_if_t<std::is_base_of<vector, LVector>::value, LVector&>
    operator/=(LVector& lhs, const Rat& scal)
    {
        // TODO: implement
        return lhs;
    }



    template <typename LVector>
    friend std::enable_if_t<std::is_base_of<vector, LVector>::value, LVector&>
    operator+=(LVector& lhs, const vector& rhs)
    {
        // TODO: implement
        return lhs;
    }

    template <typename LVector>
    friend std::enable_if_t<std::is_base_of<vector, LVector>::value, LVector&>
    operator-=(LVector& lhs, const vector& rhs)
    {
        // TODO: implement
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
    if (lhs.dimension() != rhs.dimension()) {
        return false;
    }

    const auto& lhs_bv = lhs.base_vector();
    const auto& rhs_bv = rhs.base_vector();

    for (auto right : rhs_bv) {
        if (lhs[right.key()] != right.value()) {
            return false;
        }
    }
    return true;

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
    base_type::instance()[key_type(key)] += scalar_type(scal);
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
    base_type::instance()[key_type(key)] += coefficient_ring::one() / rational_type(scal);
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
    base_type::instance()[key_type(key)] -= scalar_type(scal);
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
    base_type::instance()[key_type(key)] -= coefficient_ring::one() / rational_type(scal);
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
    // TODO: implement
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
    // TODO: implement
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
    // TODO: implement
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
    // TODO: implement
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
    // TODO: Implement
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
    // TODO: implement
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
    // TODO: implement
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
    // TODO: implement
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
    // TODO: implement
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
    // TODO: implement
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
    // TODO: implement
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
    // TODO: implement
    return *this;
}

} // namespace alg


#endif //LIBALGEBRA_LITE_VECTOR_H
