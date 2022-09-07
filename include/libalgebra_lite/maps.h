//
// Created by user on 12/08/22.
//

#ifndef LIBALGEBRA_LITE_MAPS_H
#define LIBALGEBRA_LITE_MAPS_H

#include "implementation_types.h"
#include "libalgebra_lite_export.h"

#include <boost/container/small_vector.hpp>

#include "tensor_basis.h"
#include "free_tensor.h"
#include "lie.h"

namespace lal {


namespace dtl {

class generic_commutator
{
    const tensor_basis& m_basis;
    const free_tensor_multiplication& m_mul;

public:
    using key_type = typename tensor_basis::key_type;
    using pair_type = std::pair<key_type, int>;
    using tensor_type = boost::container::small_vector<pair_type, 1>;
    using ref_type = const boost::container::small_vector_base<pair_type>&;

    generic_commutator(const tensor_basis& basis, const free_tensor_multiplication& mul)
        : m_basis(basis), m_mul(mul)
    {}

    tensor_type operator()(ref_type lhs, ref_type rhs);


};


} // namespace dtl


class LIBALGEBRA_LITE_EXPORT maps
{
    std::shared_ptr<const tensor_basis> p_tensor_basis;
    std::shared_ptr<const hall_basis> p_lie_basis;
    std::shared_ptr<const lie_multiplication> p_lie_mul;
    std::shared_ptr<const free_tensor_multiplication> p_ftensor_mul;

    dtl::generic_commutator m_commutator;

public:
    using lkey_type = typename hall_basis::key_type;
    using tkey_type = typename tensor_basis::key_type;
    using generic_scalar_type = int;

    using lie_pair = std::pair<lkey_type, generic_scalar_type>;
    using tensor_pair = std::pair<tkey_type, generic_scalar_type>;
    using generic_lie = boost::container::small_vector<lie_pair, 1>;
    using generic_tensor = boost::container::small_vector<tensor_pair, 1>;

    using gtensor_ref = const boost::container::small_vector_base<tensor_pair>&;
    using glie_ref = const boost::container::small_vector_base<lie_pair>&;

private:

    static generic_tensor expand_letter(let_t letter);


    hall_extension<decltype(&maps::expand_letter),
            dtl::generic_commutator,
            gtensor_ref> m_expand;


    generic_lie rbracketing_impl(tkey_type arg) const;

public:

    maps(deg_t width, deg_t depth)
        : p_tensor_basis(basis_registry<tensor_basis>::get(width, depth)),
          p_lie_basis(basis_registry<hall_basis>::get(width, depth)),
          p_ftensor_mul(multiplication_registry<free_tensor_multiplication>::get(width)),
          p_lie_mul(multiplication_registry<lie_multiplication>::get(width)),
          m_commutator(*p_tensor_basis, *p_ftensor_mul)
    {}


    glie_ref rbracketing(tkey_type tkey) const;
    gtensor_ref expand(lkey_type lkey) const;


    template <typename Coefficients,
              template <typename, typename> class VectorType,
              template <typename> class StorageModel>
    free_tensor<Coefficients, VectorType, StorageModel>
    lie_to_tensor(const lie<Coefficients, VectorType, StorageModel>& arg) const
    {
        if (arg.basis().width() != p_lie_basis->width()) {
            throw std::invalid_argument("mismatched width");
        }

        auto max_deg = p_lie_basis->depth();
        free_tensor<Coefficients, VectorType, StorageModel> result(p_tensor_basis);
        if (arg.basis().depth() <= max_deg) {
            for (auto outer : arg) {
                auto val = outer.value();
                for (auto inner : expand(outer.key())) {
                    result.add_scal_prod(inner.first, Coefficients::mul(inner.second, val));
                }
            }
        } else {
            for (auto outer : arg) {
                auto key = outer.key();
                auto val = outer.value();
                if (p_lie_basis->degree(key) <= max_deg) {
                    for (auto inner : expand(key)) {
                        result.add_scal_prod(inner.first, Coefficients::mul(inner.second, val));
                    }
                }
            }
        }
        return result;
    }

    template <typename Coefficients,
            template <typename, typename> class VectorType,
            template <typename> class StorageModel>
    lie<Coefficients, VectorType, StorageModel>
    tensor_to_lie(const free_tensor<Coefficients, VectorType, StorageModel>& arg) const
    {
        using scalar_type = typename Coefficients::scalar_type;

        if (arg.basis().width() != p_tensor_basis->width()) {
            throw std::invalid_argument("mismatched width");
        }
        auto max_deg = p_tensor_basis->depth();

        lie<Coefficients, VectorType, StorageModel> result(p_lie_basis);
        if (arg.basis.depth() <= max_deg) {
            for (auto outer : arg) {
                auto val = outer.value();
                for (auto inner : rbracketing(outer.key())) {
                    result.add_scal_prod(inner.key(), Coefficients::mul(scalar_type(inner.second), val));
                }
            }
        } else {
            for (auto outer : arg) {
                auto key = outer.key();
                auto val = outer.value();
                if (p_tensor_basis->degree(key) <= max_deg) {
                    for (auto inner : rbracketing(key)) {
                        result.add_scal_prod(inner.key(), Coefficients::mul(scalar_type(inner.second), val));
                    }
                }
            }
        }
        return result;
    }
};


} // namespace alg



#endif //LIBALGEBRA_LITE_MAPS_H
