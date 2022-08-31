//
// Created by user on 12/08/22.
//

#ifndef LIBALGEBRA_LITE_LIE_H
#define LIBALGEBRA_LITE_LIE_H

#include <libalgebra_lite/implementation_types.h>
#include <libalgebra_lite/algebra.h>
#include <libalgebra_lite/basis/hall_set.h>

#include <mutex>
#include <unordered_map>

#include <boost/functional/hash.hpp>
#include <boost/container/small_vector.hpp>

#ifndef LIBALGEBRA_LITE_EXPORT
#define LIBALGEBRA_LITE_EXPORT
#endif


namespace lal {


class LIBALGEBRA_LITE_EXPORT lie_multiplier
{
    std::shared_ptr<const hall_basis> p_basis;


    using key_type = typename hall_basis::key_type;
    using parent_type = std::pair<key_type, key_type>;

    using product_type = boost::container::small_vector<std::pair<key_type, int>, 1>;
    using product_ref_type = const boost::container::small_vector_base<std::pair<key_type, int>>&;

    mutable std::unordered_map<parent_type, product_type, boost::hash<parent_type>> m_cache;
    mutable std::recursive_mutex m_lock;

    product_type key_prod_impl(key_type lhs, key_type rhs) const;
    static product_type generic_minus(product_ref_type arg);
    static product_type generic_sub(product_ref_type lhs, product_ref_type rhs);

    product_type generic_by_key(product_ref_type lhs, key_type rhs) const;


public:

    explicit lie_multiplier(std::shared_ptr<const hall_basis> basis)
        : p_basis(std::move(basis)), m_lock(), m_cache()
    {}

    product_ref_type operator()(key_type lhs, key_type rhs) const;


};

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
