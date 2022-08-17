//
// Created by user on 26/07/22.
//

#ifndef LIBALGEBRA_LITE_ALGEBRA_H
#define LIBALGEBRA_LITE_ALGEBRA_H


#include <libalgebra_lite/implementation_types.h>
#include <libalgebra_lite/basis/traits.h>
#include <libalgebra_lite/coefficients.h>
#include <libalgebra_lite/vectors/traits.h>
#include <libalgebra_lite/vector.h>

#include <algorithm>
#include <memory>
#include <type_traits>
#include <utility>
#include <vector>

#include <boost/tti/has_template.hpp>

namespace alg {


template <typename Multiplication>
struct multiplication_traits
{

    using mult_ptr = std::shared_ptr<const Multiplication>;


private:

    template <typename, typename, typename, typename, typename = void>
    struct has_fma_inplace : std::false_type {};

    template <typename L, typename R, typename F, typename M>
    struct has_fma_inplace<L, R, F, M,
            std::void_t<decltype(M::template fma_inplace(std::declval<L&>(), std::declval<const R&>(), std::declval<F>()))>>
        : std::true_type
    {};

    template <typename, typename, typename, typename, typename=void>
    struct has_fma_inplace_with_deg : std::false_type {};

    template <typename L, typename R, typename F, typename M>
    struct has_fma_inplace_with_deg<L, R, F, M,
            std::void_t<decltype(M::template fma_inplace(std::declval<L&>(), std::declval<const R&>(), std::declval<F>()), std::declval<deg_t>())>>
        : std::true_type
    {};

    template <typename V>
    using has_degree_t = std::is_same<typename basis_trait<typename V::basis_type>::degree_tag, with_degree_tag>;

public:
    template <typename Result,
              typename Vector1,
              typename Vector2,
              typename Fn>
    static std::enable_if_t<!has_degree_t<Result>::value>
    multiply_and_add(
            mult_ptr mult,
            Result& result,
            const Vector1& lhs,
            const Vector2& rhs,
            Fn fn)
    {
        if (lhs.empty() || rhs.empty()) {
            return;
        }

        mult->fma(result.base_vector(),
                lhs.base_vector(),
                rhs.base_vector(),
                fn);
    }

    template <typename Result,
            typename Vector1,
            typename Vector2>
    static std::enable_if_t<!has_degree_t<Result>::value>
    multiply_and_add(mult_ptr mult,
            Result& result,
            const Vector1& lhs,
            const Vector2& rhs)
    {
        using scalar_type = typename Result::scalar_type;
        mult->fma(result.base_vector(),
                lhs.base_vector(),
                rhs.base_vector(),
                [](scalar_type s) { return s; });
    }

    template <typename Result,
              typename Vector1,
              typename Vector2,
              typename Fn>
    static std::enable_if_t<has_degree_t<Result>::value>
    multiply_and_add(
            mult_ptr mult,
            Result& result,
            const Vector1& lhs,
            const Vector2& rhs,
            Fn fn)
    {
        using traits = basis_trait<typename Result::basis_type>;
        mult->fma(result.base_vector(),
                lhs.base_vector(),
                rhs.base_vector(),
                fn,
                traits::max_degree(result.basis()));
    }

    template <typename Result,
            typename Vector1,
            typename Vector2>
    static std::enable_if_t<has_degree_t<Result>::value>
    multiply_and_add(mult_ptr mult,
            Result& result,
            const Vector1& lhs,
            const Vector2& rhs)
    {
        using scalar_type = typename Result::scalar_type;
        using traits = basis_trait<typename Result::basis_type>;
        mult->fma(result.base_vector(),
                lhs.base_vector(),
                rhs.base_vector(),
                [](scalar_type s) { return s; },
                traits::max_degree(result.basis()));
    }

    template <typename Left, typename Right, typename Fn, typename Mult=Multiplication>
    static std::enable_if_t<!has_degree_t<Left>::value && has_fma_inplace<Left, Right, Fn, Mult>::value>
    multiply_and_add_inplace(mult_ptr mult, Left& lhs, const Right& rhs, Fn fn)
    {
        mult->fma_inplace(lhs.vector_type(), rhs.vector_type(), fn);
    }

    template <typename Left, typename Right, typename Fn, typename Mult=Multiplication>
    static std::enable_if_t<has_degree_t<Left>::value && has_fma_inplace<Left, Right, Fn, Mult>::value>
    multiply_and_add_inplace(mult_ptr mult, Left& lhs, const Right& rhs, Fn fn)
    {
        using traits = basis_trait<typename Left::basis_type>;
        mult->fma_inplace(lhs.base_vector(), rhs.base_vector(), fn,
                traits::max_degree(lhs.basis()));
    }

    template <typename Left, typename Right, typename Fn, typename Mult=Multiplication>
    static std::enable_if_t<has_degree_t<Left>::value && has_fma_inplace<Left, Right, Fn, Mult>::value>
    multiply_and_add_inplace(mult_ptr mult, Left& lhs, const Right& rhs, Fn fn, deg_t max_deg)
    {
        mult->fma_inplace(lhs.base_vector(), rhs.base_vector(), fn, max_deg);
    }

    template <typename Left, typename Right, typename Fn>
    static std::enable_if_t<!has_degree_t<Left>::value>
    multiply_inplace(mult_ptr mult, Left& lhs, const Right& rhs, Fn fn)
    {
        Left tmp(lhs.basis(), mult);
        mult->fma_inplace(lhs.base_vector(), rhs.base_vector(), fn);
        lhs.swap(tmp);
    }


    template <typename Left, typename Right, typename Fn>
    static std::enable_if_t<has_degree_t<Left>::value>
    multiply_inplace(mult_ptr mult, Left& lhs, const Right& rhs, Fn fn, deg_t max_deg)
    {
        Left tmp(lhs.basis(), mult);
        mult->fma_inplace(lhs.base_vector(), rhs.base_vector(), fn, max_deg);
        lhs.swap(tmp);
    }

    template <typename Left, typename Right, typename Fn>
    static std::enable_if_t<has_degree_t<Left>::value>
    multiply_inplace(mult_ptr mult, Left& lhs, const Right& rhs, Fn fn)
    {
        auto max_deg = basis_trait<typename Left::basis_type>::max_degree(*lhs.mult());
        multiply_inplace(std::move(mult), lhs, rhs, fn, max_deg);
    }




};



template <typename Basis>
class basic_multiplier;

namespace dtl {

template <typename Basis, typename Coefficients>
class general_multiplication_helper
{
protected:
    using basis_traits = basis_trait<Basis>;
    using key_type = typename basis_traits::key_type;
    using coeff_traits = coefficient_trait<Coefficients>;
    using scalar_type = typename coeff_traits::scalar_type;

    using key_value = std::pair<key_type, scalar_type>;
    std::vector<key_value> right_buffer;

public:

    using const_iterator = typename std::vector<key_value>::const_iterator;

    template<template<typename, typename, typename...> class VTR, typename... RArgs>
    explicit general_multiplication_helper(const VTR<Basis, Coefficients, RArgs...>& rhs)
        : right_buffer()
    {
        right_buffer.reserve(rhs.size());
        for (auto item : rhs) {
            right_buffer.emplace_back(item.key(), item.value());
        }
    }

    const_iterator begin() const noexcept { return right_buffer.begin(); }
    const_iterator end() const noexcept { return right_buffer.end(); }


};


template <typename It>
class degree_range_iterator_helper
{
    It m_begin, m_end;

public:
    using iterator = It;

    degree_range_iterator_helper(It begin, It end)
        : m_begin(begin), m_end(end)
    {}

    iterator begin() noexcept { return m_begin; }
    iterator end() noexcept { return m_end; }
};

template <typename Basis, typename Coefficients>
class graded_multiplication_helper : protected general_multiplication_helper<Basis, Coefficients>
{
   using base_type = general_multiplication_helper<Basis, Coefficients>;
   using ordering = typename base_type::basis_traits::kv_ordering;

   using typename base_type::basis_traits;
   using typename base_type::key_type;
   using typename base_type::scalar_type;
   using typename base_type::key_value;
   using base_type::right_buffer;

   using base_iter = typename std::vector<key_value>::const_iterator;

   std::vector<base_iter> degree_ranges;
   deg_t max_degree;

public:

    using iterable = degree_range_iterator_helper<base_iter>;

    template<template<typename, typename, typename...> class VTR, typename... RArgs>
    explicit graded_multiplication_helper(const VTR<Basis, Coefficients, RArgs...>& rhs)
        : base_type(rhs)
    {
        ordering order;
        std::sort(right_buffer.begin(), right_buffer.end(), order);

        const auto& basis = rhs.basis();
        max_degree = basis_traits::max_degree(basis);
        degree_ranges.reserve(max_degree + 1);
        auto it = right_buffer.cbegin();
        auto end = right_buffer.cend();
        degree_ranges.push_back(it);
        auto degree = 0;

        for (; it!=end; ++it) {
            auto current = basis_traits::degree(basis, it->first);
            for (;degree < current; ++degree) {
                degree_ranges.push_back(it);
            }
        }
        for (; degree<=max_degree; ++degree) {
            degree_ranges.push_back(end);
        }
    }

    iterable degree_range(deg_t degree) noexcept
    {
        if (degree < 0) {
            return iterable(right_buffer.begin(), right_buffer.begin());
        }
        if (degree > max_degree) {
            return iterable(right_buffer.end(), right_buffer.end());
        }
        return iterable(degree_ranges[degree], degree_ranges[degree+1]);
    }

};

} // namespace dtl


template <typename Multiplier>
class base_multiplication
{

    Multiplier m_mult;

    template <typename V>
    using basis_t = typename V::basis_type;

    template <typename V>
    using key_tp = typename basis_trait<basis_t<V>>::key_type;

    template <typename V>
    using scal_t = typename V::scalar_type;

    template <typename V, typename Sca>
    using key_vect = std::vector<
            std::pair<typename basis_trait<typename V::basis_type>::key_type,
                      Sca>>;

    template <typename V>
    using helper_type = dtl::general_multiplication_helper<
            typename V::basis_type, typename V::coefficient_ring>;
    template <typename V>
    using graded_helper_type = dtl::graded_multiplication_helper<
            typename V::basis_type, typename V::coefficient_ring>;

    template <typename OutVector, typename KeyProd, typename Sca>
    void asp_helper(OutVector& out, KeyProd&& key_prod, Sca&& scalar) const
    {
        using scalar_type = scal_t<OutVector>;
        out.add_scal_prod(std::forward<KeyProd>(key_prod),
                scalar_type(std::forward<Sca>(scalar)));
    }

    template <typename OutVector, typename KSca, typename PSca>
    void asp_helper(OutVector& out, key_vect<OutVector, KSca>&& key_prod, PSca&& scalar) const
    {
        using scalar_type = scal_t<OutVector>;
        for (const auto& kv_pair : key_prod) {
            out.add_scal_prod(kv_pair.first,
                    scalar_type(kv_pair.second)*scalar_type(std::forward<PSca>(scalar)));

        }
    }


public:

    template <typename OutVector,
            typename LeftVector,
            typename RightVector,
            typename Fn>
    void fma(OutVector& out, const LeftVector& lhs, const RightVector& rhs, Fn fn) const
    {
        // The helper makes a contiguous copy of rhs key-value pairs
        // so the inner loop is always contiguous, rather than potentially
        // a linked list or other cache unfriendly data structure.
        helper_type<RightVector> helper(rhs);

        for (auto litem : lhs) {
            for (auto ritem : helper) {
                asp_helper(out, m_mult(litem.key(), ritem.first),
                                fn(litem.value()*ritem.second));
            }
        }
    }

    template <typename OutVector,
            typename LeftVector,
            typename RightVector,
            typename Fn>
    void fma(OutVector& out, const LeftVector& lhs, const RightVector& rhs, Fn fn, deg_t max_deg) const
    {
        using out_basis_traits = basis_trait<typename OutVector::basis_type>;
        using lhs_basis_traits = basis_trait<typename LeftVector::basis_type>;

        // The helper makes a contiguous copy of rhs key-value pairs
        // so the inner loop is always contiguous, rather than potentially
        // a linked list or other cache unfriendly data structure.
        // The graded helper also sorts the keys and constructs a buffer
        // of degree ranges that we can use to truncate products that
        // would overflow the max degree.
        graded_helper_type<RightVector> helper(rhs);

        auto out_deg = std::min(out_basis_traits::max_degree(out.basis()), lhs.degree() + rhs.degree());

        const auto& lhs_basis = lhs.basis();
        for (auto litem : lhs) {
            auto lhs_degree = lhs_basis_traits::degree(lhs_basis, lhs_basis);
            auto rhs_degree = out_deg - lhs_degree;
            for (auto ritem : helper.degree_range(rhs_degree)) {
                asp_helper(out, m_mult(litem.key(), ritem.first),
                                fn(litem.value(), ritem.second));
            }
        }
    }






};



template <typename Basis,
        typename Coefficients,
        typename Multiplication,
        template <typename, typename> class VectorType,
        template <typename> class StorageModel
        >
class algebra : public vector<Basis, Coefficients, VectorType, StorageModel>
{
public:

    using vector_type = vector<Basis, Coefficients, VectorType, StorageModel>;
    using multiplication_type = Multiplication;

    using typename vector_type::basis_type;
    using typename vector_type::key_type;
    using typename vector_type::coefficient_ring;
    using typename vector_type::scalar_type;
    using typename vector_type::rational_type;

private:
    std::shared_ptr<const multiplication_type> p_mult;

public:
    algebra(const vector_type& base, std::shared_ptr<Multiplication> mult)
        : vector_type(base), p_mult(std::move(mult))
    {}

    template <template <typename, typename> class OtherVectorType,
    template <typename> class OtherStorageModel>
    explicit algebra(const algebra<Basis, Coefficients, Multiplication, OtherVectorType, OtherStorageModel>& other)
        : vector_type((other)), p_mult(other.p_mult)
    {}


    std::shared_ptr<multiplication_type> multiplication() const noexcept { return p_mult; }


    algebra& add_mul(const algebra& lhs, const algebra& rhs);
    algebra& sub_mul(const algebra& lhs, const algebra& rhs);

    algebra& mul_scal_prod(const algebra& lhs, const scalar_type& scal);
    algebra& mul_scal_div(const algebra& lhs, const rational_type& scal);



};


namespace dtl {

template <typename Algebra>
class is_algebra
{
    template <typename Basis,
            typename Coefficients,
            typename Multiplication,
            template <typename, typename> class VType,
            template <typename> class SModel>
    static std::true_type test(algebra<Basis, Coefficients, Multiplication, VType, SModel>&);

    static std::false_type test(...);

public:
    static constexpr bool value = decltype(test(std::declval<Algebra&>()))::value;
};

} // namespace dtl

template <typename Algebra>
std::enable_if_t<dtl::is_algebra<Algebra>::value, Algebra>
operator*(const Algebra& lhs, const Algebra& rhs)
{
    using traits = multiplication_traits<typename Algebra::multiplication_type>;

    auto multiplication = lhs.multiplication();
    if (!multiplication) {
        multiplication = rhs.multiplication();
    }
    Algebra result(lhs.basis, multiplication);
    if (multiplication && !lhs.empty() && !rhs.empty()) {
        traits::multiply_and_add(std::move(multiplication), lhs, rhs);
    }
    return result;
}

template <typename Algebra>
std::enable_if_t<dtl::is_algebra<Algebra>::value, Algebra&>
operator*=(Algebra& lhs, const Algebra& rhs)
{
    using traits = multiplication_traits<typename Algebra::multiplication_type>;
    if (rhs.empty()) {
        lhs.clear();
    }
    auto multiplication = lhs.multiplication();
    if (!multiplication) {
        multiplication = rhs.multiplication();
    }

    if (multiplication && !lhs.empty()) {
        traits::multiply_inplace(lhs.multiplication, lhs, rhs);
    }
    return lhs;
}

template <typename Algebra>
std::enable_if_t<dtl::is_algebra<Algebra>::value, Algebra>
commutator(const Algebra& lhs, const Algebra& rhs)
{
    using traits = multiplication_traits<typename Algebra::multiplication_type>;
    auto multiplication = lhs.multiplication();
    if (!multiplication) {
        multiplication = rhs.multiplication();
    }

    Algebra result(lhs.basis(), multiplication);
    if (multiplication && !lhs.empty() && !rhs.empty()) {
        traits::mulitply_and_add(*multiplication, result, lhs, rhs);
        traits::multiply_and_add(*multiplication, result, rhs, lhs, Algebra::coefficient_ring::uminus);
    }
    return result;
}



} // namespace alg

#endif //LIBALGEBRA_LITE_ALGEBRA_H
