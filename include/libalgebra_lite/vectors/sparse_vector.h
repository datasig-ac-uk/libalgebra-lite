//
// Created by user on 26/07/22.
//

#ifndef LIBALGEBRA_LITE_SPARSE_VECTOR_H
#define LIBALGEBRA_LITE_SPARSE_VECTOR_H

#include <libalgebra_lite/implementation_types.h>
#include <libalgebra_lite/basis/traits.h>
#include <libalgebra_lite/coefficients.h>

#include <boost/container/flat_map.hpp>

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
    MapType& m_map;
    typename MapType::iterator m_it;
    typename MapType::mapped_type m_tmp;
    KeyType m_key;

    using Self = sparse_mutable_reference;

public:

    using key_type = KeyType;
    using scalar_type = typename MapType::mapped_type;

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

} // namespace dtl


template <typename Basis, typename Coefficients>
class sparse_vector
{
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


public:

    using reference = dtl::sparse_mutable_reference<map_type, key_type>;



};


} // namespace lal

#endif //LIBALGEBRA_LITE_SPARSE_VECTOR_H
