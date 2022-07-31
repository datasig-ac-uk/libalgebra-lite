//
// Created by user on 26/07/22.
//

#ifndef LIBALGEBRA_LITE_ALGEBRA_H
#define LIBALGEBRA_LITE_ALGEBRA_H


#include <libalgebra_lite/implementation_types.h>
#include <libalgebra_lite/basis/traits.h>
#include <libalgebra_lite/coefficients.h>
#include <libalgebra_lite/vectors/traits.h>

#include <memory>
#include <type_traits>

namespace alg {


template <typename Vector, typename Multiplication>
class algebra : public Vector
{
    using base_vector = vector_base<Vector>;

public:
    using vector_type = Vector;
    using multiplication_type = Multiplication;

    using typename vector_type::basis_type;
    using typename vector_type::key_type;
    using typename vector_type::coefficient_ring;
    using typename vector_type::scalar_type;
    using typename vector_type::rational_type;

private:
    std::shared_ptr<const multiplication_type> p_mult;

public:
    algebra(const base_vector& base, std::shared_ptr<Multiplication> mult)
        : Vector(base), p_mult(std::move(mult))
    {}

    template <typename OtherVector,
            typename=typename std::enable_if_t<std::is_base_of<base_vector, OtherVector>::value>>
    explicit algebra(const algebra<OtherVector, Multiplication>& other)
        : Vector(static_cast<const OtherVector&>(other)), p_mult(other.p_mult)
    {}

    template <typename Base=base_vector,
            typename=typename std::enable_if_t<!std::is_same<Vector, Base>::value>>
    operator algebra<Base, Multiplication>()
    {
        return algebra<Base, Multiplication>(static_cast<base_vector&>(*this), p_mult);
    }

    template <typename Base=base_vector,
            typename=typename std::enable_if_t<!std::is_same<Vector, Base>::value>>
    operator const algebra<Base, Multiplication> () const
    {
        return algebra<Base, Multiplication>(static_cast<const base_vector&>(*this), p_mult);
    }

    template <typename Owned=typename base_vector::owned_vector_type,
            typename=typename std::enable_if_t<!std::is_same<Vector, Owned>::value>>
    explicit operator algebra<Owned, Multiplication>() const
    {
        using new_alg = algebra<Owned, Multiplication>;
        return new_alg(static_cast<const base_vector&>(*this), p_mult);
    }

    const multiplication_type& multiplication() const noexcept { return *p_mult; }

};


namespace dtl {

template <typename Algebra>
class is_algebra
{
    template <typename Vector, typename Multiplication>
    static std::true_type test(algebra<Vector, Multiplication>&);

    static std::false_type test(...);

public:
    static constexpr bool value = decltype(test(std::declval<Algebra&>()))::value;
};

} // namespace dtl

template <typename Algebra>
std::enable_if_t<dtl::is_algebra<Algebra>::value,
    typename Algebra::vector_type::owned_vector_type>
operator*(const Algebra& lhs, const Algebra& rhs)
{
    using base_t = vector_base<typename Algebra::vector_type>;
    using owned_t = typename base_t::owned_vector_type;

    owned_t result;
    const auto& mult = lhs.multiplication();
    mult.multiply_and_add(result, lhs, rhs);
    return result;
}

template <typename Algebra>
std::enable_if_t<dtl::is_algebra<Algebra>::value, Algebra&>
operator*=(Algebra& lhs, const Algebra& rhs)
{
    const auto& mult = lhs.multiplication();
    mult.multiply_inplace(lhs, rhs);
}




} // namespace alg

#endif //LIBALGEBRA_LITE_ALGEBRA_H
