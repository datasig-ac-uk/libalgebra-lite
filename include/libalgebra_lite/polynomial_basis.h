//
// Created by user on 28/08/22.
//

#ifndef LIBALGEBRA_LITE_POLYNOMIAL_BASIS_H
#define LIBALGEBRA_LITE_POLYNOMIAL_BASIS_H

#include "implementation_types.h"

#include <functional>
#include <numeric>

#include <boost/container/small_vector.hpp>
#include <boost/container/flat_map.hpp>

#include "libalgebra_lite_export.h"
#include "packed_integer.h"

namespace lal {


class LIBALGEBRA_LITE_EXPORT monomial
{
    using letter_type = dtl::packed_integer<dimn_t, char>;

    using small_vec = boost::container::small_vector<std::pair<letter_type, deg_t>, 1>;
    using map_type = boost::container::flat_map<letter_type, deg_t, std::less<>, small_vec>;

    map_type m_data;

    template <typename Coeff>
    static typename Coeff::scalar_type
    power(typename Coeff::scalar_type arg, deg_t exponent) noexcept
    {
        if (exponent == 0) {
            return Coeff::one();
        }
        if (exponent == 1) {
            return arg;
        }
        auto result1 = power<Coeff>(arg, exponent/2);
        auto result2 = Coeff::mul(result1);
        return (exponent % 2==0) ? result2 : Coeff::mul(arg, result2);
    }

public:

    using iterator = typename map_type::iterator;
    using const_iterator = typename map_type::const_iterator;

    explicit monomial(letter_type let, deg_t power=1)
    {
        assert(power > 0);
        m_data[let] = power;
    }

    template <typename MapType>
    explicit monomial(const MapType& arg) : m_data(arg.begin(), arg.end())
    {}

    template <typename InputIt>
    explicit monomial(InputIt begin, InputIt end) : m_data(begin, end)
    {}

    deg_t degree() const noexcept;

    iterator begin() noexcept { return m_data.begin(); }
    iterator end() noexcept { return m_data.end(); }
    const_iterator begin() const noexcept { return m_data.begin(); }
    const_iterator end() const noexcept { return m_data.end(); }

    template <typename Coefficients, typename MapType>
    typename Coefficients::scalar_type eval(const MapType& arg) const noexcept
    {
        auto result = Coefficients::zero();

        for (const auto& item : m_data) {
            Coefficients::add_inplace(power<Coefficients>(arg[item.first], item.second));
        }
        return result;
    }

};



struct LIBALGEBRA_LITE_EXPORT polynomial_basis
{
    using letter_type = dtl::packed_integer<dimn_t, char>;
    using key_type = monomial;

    static key_type key_of_letter(letter_type letter)
    {
        return key_type(letter);
    }

    static deg_t degree(const key_type& key)
    {
        return std::accumulate(key.begin(), key.end(), 0,
                [](const auto& curr, const auto& key) { return curr + key.second; });
    }

    struct key_order
    {
        bool operator()(const key_type& lhs, const key_type& rhs) const;
    };

};


} // namespace lal

#endif //LIBALGEBRA_LITE_POLYNOMIAL_BASIS_H