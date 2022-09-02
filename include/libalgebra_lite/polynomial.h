//
// Created by user on 30/08/22.
//

#ifndef LIBALGEBRA_LITE_POLYNOMIAL_H
#define LIBALGEBRA_LITE_POLYNOMIAL_H

#include "implementation_types.h"
#include "polynomial_basis.h"
#include "sparse_vector.h"
#include "algebra.h"
#include "libalgebra_lite_export.h"
#include "coefficients.h"

#include <boost/container/small_vector.hpp>

#ifndef LIBALGEBRA_LITE_EXPORT
#define LIBALGEBRA_LITE_EXPORT
#endif

namespace lal {

class LIBALGEBRA_LITE_EXPORT polynomial_multiplier
{
    using letter_type = typename polynomial_basis::letter_type;
    using key_type = typename polynomial_basis::key_type;

    using product_type = boost::container::small_vector<std::pair<key_type, int>, 1>;

public:

    product_type operator()(const key_type& lhs, const key_type& rhs) const;

};



template <typename Coefficients>
class polynomial : public algebra<polynomial_basis,
        Coefficients,
        base_multiplication<polynomial_multiplier>,
        sparse_vector,
        dtl::standard_storage
        >
{
    using base = algebra<polynomial_basis,
                         Coefficients,
                         base_multiplication<polynomial_multiplier>,
                         sparse_vector,
                         dtl::standard_storage>;

public:

    using base::base;

    template <typename IndeterminateMap>
    typename polynomial::scalar_type operator()(const IndeterminateMap& arg) const noexcept
    {
        using ring = typename polynomial::coefficient_ring;
        auto ans = ring::zero();
        for (const auto& item : *this) {
            auto key_result = item.first.template eval<ring>(arg);
            ring::add_inplace(ans, ring::mul(item.second, key_result));
        }
        return ans;
    }


};

extern template class LIBALGEBRA_LITE_EXPORT polynomial<double_field>;
extern template class LIBALGEBRA_LITE_EXPORT polynomial<float_field>;

extern template class LIBALGEBRA_LITE_EXPORT coefficient_ring<polynomial<double_field>, double>;
extern template class LIBALGEBRA_LITE_EXPORT coefficient_ring<polynomial<float_field>, float>;
extern template class LIBALGEBRA_LITE_EXPORT coefficient_ring<polynomial<rational_field>,
        typename rational_field::scalar_type>;

}



#endif //LIBALGEBRA_LITE_POLYNOMIAL_H
