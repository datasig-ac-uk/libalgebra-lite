//
// Created by user on 12/08/22.
//

#ifndef LIBALGEBRA_LITE_MAPS_H
#define LIBALGEBRA_LITE_MAPS_H

#include "implementation_types.h"
#include "tensor_basis.h"
#include "free_tensor.h"
#include "lie.h"

namespace lal {

class LIBALGEBRA_LITE_EXPORT maps
{
    std::shared_ptr<const tensor_basis> p_tensor_basis;
    std::shared_ptr<const hall_basis> p_lie_basis;

public:
    using lkey_type = typename hall_basis::key_type;
    using tkey_type = typename tensor_basis::key_type;
    using generic_scalar_type = int;
    using lie_pair = std::pair<lkey_type, generic_scalar_type>;
    using tensor_pair = std::pair<tkey_type, generic_scalar_type>;
    using generic_lie = std::vector<lie_pair>;
    using generic_tensor = std::vector<tensor_pair>;

private:

    static generic_lie expand_letter(let_t letter);
    static generic_lie mulitply_generic_lie(const generic_lie&, const generic_lie&);

public:

    const generic_lie& rbracketing(tkey_type tkey) const;
    const generic_tensor& expand(lkey_type lkey) const;


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
