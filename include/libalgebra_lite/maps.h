//
// Created by user on 12/08/22.
//

#ifndef LIBALGEBRA_LITE_MAPS_H
#define LIBALGEBRA_LITE_MAPS_H

#include "implementation_types.h"
#include <libalgebra_lite/basis/tensor_basis.h>
#include "free_tensor.h"
#include "lie.h"

namespace alg {

class maps
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
        free_tensor<Coefficients, VectorType, StorageModel> result(p_tensor_basis);
        for (auto outer : arg) {
            for (auto inner : expand(outer.key())) {
                result.add_scal_prod(inner.first, inner.second*outer.value());
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
        lie<Coefficients, VectorType, StorageModel> result(p_lie_basis);
        for (auto outer : arg) {
            for (auto inner : rbracketing(outer.key())) {
                result.add_scal_prod()
            }
        }
    }
};


} // namespace alg



#endif //LIBALGEBRA_LITE_MAPS_H
