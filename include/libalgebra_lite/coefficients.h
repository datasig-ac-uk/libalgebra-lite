//
// Created by user on 26/07/22.
//

#ifndef LIBALGEBRA_LITE_COEFFICIENTS_H
#define LIBALGEBRA_LITE_COEFFICIENTS_H

#include "implementation_types.h"
#include "libalgebra_lite_export.h"

#include <memory>

#include "rationals.h"

namespace lal {

template <typename Coeff>
struct coefficient_trait {
    using coefficient_ring = Coeff;
    using scalar_type = typename Coeff::scalar_type;
    using rational_type = typename Coeff::rational_type;
    using default_alloc = std::allocator<scalar_type>;
};

#define LAL_RING_GENERATE_BINOP(NAME, OP, RET_T, LHS_T, RHS_T)          \
    template <typename Lhs = LHS_T, typename Rhs = RHS_T>               \
    static constexpr RET_T NAME(const Lhs &lhs, const Rhs &rhs) {       \
        return lhs OP rhs;                                              \
    }                                                                   \
                                                                        \
    template <typename Rhs = RHS_T>                                     \
    static constexpr RET_T NAME##_inplace(RET_T &lhs, const Rhs &rhs) { \
        return (lhs OP## = rhs);                                        \
    }

template <typename Scalar, typename Rational>
struct coefficient_ring {
    using scalar_type = Scalar;
    using rational_type = Rational;

    static const scalar_type &zero() noexcept {
        static const scalar_type zero{};
        return zero;
    }
    static const scalar_type &one() noexcept {
        static const scalar_type one(1);
        return one;
    }
    static const scalar_type &mone() noexcept {
        static const scalar_type mone(-1);
        return mone;
    }

    static constexpr scalar_type uminus(const scalar_type &arg) {
        return -arg;
    }

    LAL_RING_GENERATE_BINOP(add, +, scalar_type, scalar_type, scalar_type)
    LAL_RING_GENERATE_BINOP(sub, -, scalar_type, scalar_type, scalar_type)
    LAL_RING_GENERATE_BINOP(mul, *, scalar_type, scalar_type, scalar_type)
    LAL_RING_GENERATE_BINOP(div, /, scalar_type, scalar_type, rational_type)
};

#undef LAL_RING_GENERATE_BINOP

template <typename Scalar>
struct coefficient_field : public coefficient_ring<Scalar, Scalar> {
};

LAL_EXPORT_TEMPLATE(coefficient_field, double)
LAL_EXPORT_TEMPLATE(coefficient_field, float)

using double_field = coefficient_field<double>;
using float_field = coefficient_field<float>;

LAL_EXPORT_TEMPLATE(coefficient_field, dtl::rational_scalar_type)
using rational_field = coefficient_field<dtl::rational_scalar_type>;

template <>
struct coefficient_trait<float> {
    using coefficient_ring = float_field;
    using scalar_type = float;
    using rational_type = float;
};

template <>
struct coefficient_trait<double> {
    using coefficient_ring = double_field;
    using scalar_type = double;
    using rational_type = double;
};

template <>
struct coefficient_trait<dtl::rational_scalar_type> {
    using coefficient_ring = rational_field;
    using scalar_type = dtl::rational_scalar_type;
    using rational_type = dtl::rational_scalar_type;
};

}// namespace lal

#endif//LIBALGEBRA_LITE_COEFFICIENTS_H
