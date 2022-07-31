//
// Created by user on 26/07/22.
//

#ifndef LIBALGEBRA_LITE_COEFFICIENTS_H
#define LIBALGEBRA_LITE_COEFFICIENTS_H

#include <libalgebra_lite/implementation_types.h>

#include <memory>

namespace alg {


template <typename Coeff>
struct coefficient_trait
{
    using coefficient_ring = Coeff;
    using scalar_type = typename Coeff::scalar_type;
    using rational_type = typename Coeff::rational_type;
    using default_alloc = std::allocator<scalar_type>;

    static const scalar_type& zero() noexcept
    {
        static const scalar_type zero;
        return zero;
    }
    static const scalar_type& one() noexcept
    {
        static const scalar_type one(1);
        return one;
    }
    static const scalar_type& mone() noexcept
    {
        static const scalar_type mone(-1);
        return mone;
    }


};


template <typename Scalar, typename Rational>
struct coefficient_ring
{
    using scalar_type = Scalar;
    using rational_type = Rational;
};

template <typename Scalar>
struct coefficient_field : public coefficient_ring<Scalar, Scalar>
{
};


using double_field = coefficient_field<double>;
using float_field = coefficient_field<float>;


template <>
struct coefficient_trait<float>
{
    using coefficient_ring = float_field;
    using scalar_type = float;
    using rational_type = float;
};

template <>
struct coefficient_trait<double>
{
    using coefficient_ring = double_field;
    using scalar_type = double;
    using rational_type = double;
};


}

#endif //LIBALGEBRA_LITE_COEFFICIENTS_H
