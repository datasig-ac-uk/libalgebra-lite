//
// Created by user on 06/02/23.
//

#ifndef LIBALGEBRA_LITE_INCLUDE_LIBALGEBRA_LITE_OPERATORS_H
#define LIBALGEBRA_LITE_INCLUDE_LIBALGEBRA_LITE_OPERATORS_H

#include "implementation_types.h"

#include <utility>
#include <type_traits>

#include "vector.h"

namespace lal {

namespace dtl {
template <typename T> struct is_vector;
}

namespace operators {

#define LAL_LO_DECLRESULT(IMPL) \
    decltype(std::declval<IMPL>()(std::declval<const ArgumentType&>()))

namespace dtl {

template <typename Impl1, typename Impl2>
class sum_operator;

template <typename Impl, typename ArgumentType>
class left_scalar_multiply_operator;

template <typename Impl, typename ArgumentType>
class right_scalar_multiply_operator;

template <typename Operator>
struct operator_traits;


} // namespace dtl


template <typename Impl>
class linear_operator : protected Impl
{
    using traits = dtl::operator_traits<Impl>;
public:
    using result_type       = typename traits::result_type;
    using argument_type     = typename traits::argument_type;


protected:
    using implementation_type = Impl;

public:

    template <typename... Args>
    explicit linear_operator(Args&&... args) : Impl(std::forward<Args>(args)...)
    {}

    using implementation_type::operator();


    template <typename Impl2>
    friend dtl::sum_operator<Impl, Impl2>
    operator+(const linear_operator& left,
              const linear_operator<Impl2>& right)
    {
        // Sum passes on the implementation to the sum operator
        return { static_cast<const Impl&>(left), static_cast<const Impl2&>(right) };
    }

//    template <typename Scalar>
//    friend dtl::left_scalar_multiply_operator<Impl, ArgumentType>
//    operator*(const linear_operator& op, Scalar scal)
//    {
//        return { scal, static_cast<const Impl&>(op) };
//    }




};




namespace dtl {

template <typename Impl1, typename Impl2>
class sum_operator
{
    Impl1 m_left;
    Impl2 m_right;

    template <typename LeftType, typename RightType>
    static LeftType add_results(LeftType&& left, RightType&& right) {
        for (auto rit : right) {
            left.add_scal_prod(rit.key(), rit.value());
        }
        return left;
    }

    template <typename Type>
    static Type add_results(Type&& left, Type&& right) {
        left += right;
        return left;
    }

public:

    sum_operator(const Impl1& left, const Impl2& right)
        : m_left(left), m_right(right)
    {}

    sum_operator(Impl1&& left, Impl2&& right)
        : m_left(std::move(left)), m_right(std::move(right))
    {}

    template <typename Argument>
    auto operator()(const Argument &arg) const -> decltype(add_results(m_left(arg), m_right(arg)))
    {
        return add_results(m_left(arg), m_right(arg));
    }

};


template <typename Impl, typename ArgumentType>
class left_scalar_multiply_operator
{
    static_assert(lal::dtl::is_vector<ArgumentType>::value,
        "argument type must be a vector");

    using scalar_type = typename ArgumentType::scalar_type;

    Impl m_operator;
    scalar_type m_scalar;

public:
    using argument_type = ArgumentType;
    using result_type = LAL_LO_DECLRESULT(Impl);

    template <typename Scalar>
    left_scalar_multiply_operator(Scalar&& s, Impl&& op)
        : m_operator(std::move(op)), m_scalar(std::forward<Scalar>(s))
    {}

    result_type operator()(const argument_type& arg) const
    {
        return m_scalar*m_operator(arg);
    }

};

template <typename Impl, typename ArgumentType>
class right_scalar_multiply_operator
{
    static_assert(lal::dtl::is_vector<ArgumentType>::value,
        "argument type must be a vector");

    using scalar_type = typename ArgumentType::scalar_type;

    Impl m_operator;
    scalar_type m_scalar;

public:
    using argument_type = ArgumentType;
    using result_type = LAL_LO_DECLRESULT(Impl);

    template <typename Scalar>
    right_scalar_multiply_operator(Scalar &&s, const Impl &op)
        : m_operator(op), m_scalar(std::forward<Scalar>(s)) {}

    template <typename Scalar>
    right_scalar_multiply_operator(Scalar&& s, Impl&& op)
        : m_operator(std::move(op)), m_scalar(std::forward<Scalar>(s))
    {}

    result_type operator()(const argument_type& arg) const
    {
        return m_operator(arg) * m_scalar;
    }

};



} // namespace dtl


#undef LAL_LO_DECLRESULT





} // namespace operators
} // namespace lal


#endif //LIBALGEBRA_LITE_INCLUDE_LIBALGEBRA_LITE_OPERATORS_H
