//
// Created by user on 12/08/22.
//

#ifndef LIBALGEBRA_LITE_LIE_H
#define LIBALGEBRA_LITE_LIE_H

#include <libalgebra_lite/implementation_types.h>
#include <libalgebra_lite/algebra.h>
#include <libalgebra_lite/basis/hall_set.h>


#ifndef LIBALGEBRA_LITE_EXPORT
#define LIBALGEBRA_LITE_EXPORT
#endif


namespace alg {

class LIBALGEBRA_LITE_EXPORT lie_multiplication
{
    deg_t m_width, m_depth;

};


template<typename Coefficients,
        template<typename, typename> class VectorType,
        template<typename> class StorageModel>
using lie = algebra<void,
                    Coefficients,
                    lie_multiplication,
                    VectorType,
                    StorageModel>;


} // namespace alg

#endif //LIBALGEBRA_LITE_LIE_H
