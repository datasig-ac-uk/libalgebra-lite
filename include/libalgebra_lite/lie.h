//
// Created by user on 12/08/22.
//

#ifndef LIBALGEBRA_LITE_LIE_H
#define LIBALGEBRA_LITE_LIE_H

#include "implementation_types.h"
#include "libalgebra_lite_export.h"

#include <mutex>
#include <unordered_map>

#include <boost/functional/hash.hpp>
#include <boost/container/small_vector.hpp>

#include "algebra.h"
#include "hall_set.h"


namespace lal {



class LIBALGEBRA_LITE_EXPORT lie_multiplier : public base_multiplier<lie_multiplier, hall_basis, 2>
{
    using base_type = base_multiplier<lie_multiplier, hall_basis, 2>;

    std::shared_ptr<const hall_basis> p_basis;

    using typename base_type::key_type;
    using typename base_type::product_type;
    using typename base_type::reference;

    using parent_type = std::pair<key_type, key_type>;

    mutable std::unordered_map<parent_type, product_type, boost::hash<parent_type>> m_cache;
    mutable std::recursive_mutex m_lock;

    product_type key_prod_impl(key_type lhs, key_type rhs) const;


public:


    explicit lie_multiplier(std::shared_ptr<const hall_basis> basis)
        : p_basis(std::move(basis)), m_lock(), m_cache()
    {}

    reference operator()(key_type lhs, key_type rhs) const;


};

extern template class LIBALGEBRA_LITE_EXPORT base_multiplier<lie_multiplier, hall_basis>;


struct LIBALGEBRA_LITE_EXPORT lie_multiplication : public base_multiplication<lie_multiplier>
{
    using base = base_multiplication<lie_multiplier>;

    lie_multiplication(std::shared_ptr<const hall_basis> basis)
        : base(std::move(basis))
    {}
};



template<typename Coefficients,
        template<typename, typename> class VectorType,
        template<typename> class StorageModel>
using lie = algebra<hall_basis,
                    Coefficients,
                    lie_multiplication,
                    VectorType,
                    StorageModel>;


} // namespace lal

#endif //LIBALGEBRA_LITE_LIE_H