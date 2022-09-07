//
// Created by user on 05/09/22.
//

#ifndef LIBALGEBRA_LITE_REGISTRY_H
#define LIBALGEBRA_LITE_REGISTRY_H

#include "implementation_types.h"
#include "libalgebra_lite_export.h"

#include <memory>
#include <mutex>
#include <unordered_map>
#include <utility>

#include <boost/functional/hash.hpp>

namespace lal {

template <typename Basis>
class basis_registry
{
    static std::mutex m_lock;
    static std::unordered_map<std::pair<deg_t, deg_t>,
        std::shared_ptr<const Basis>, boost::hash<std::pair<deg_t, deg_t>>> m_cache;

public:

    /*
     * Most of the basis types require a width and depth pair
     */

    static std::shared_ptr<const Basis> get(deg_t width, deg_t depth)
    {
        std::lock_guard<std::mutex> access(m_lock);

        auto& found = m_cache[{width, depth}];
        if (found) {
            return found;
        }
        return found = std::make_shared<const Basis>(width, depth);
    }

};

template <typename Basis>
std::mutex basis_registry<Basis>::m_lock;

template <typename Basis>
std::unordered_map<std::pair<deg_t, deg_t>, std::shared_ptr<const Basis>,
        boost::hash<std::pair<deg_t, deg_t>>>
    basis_registry<Basis>::m_cache;


template <typename Multiplication>
class multiplication_registry
{
    static std::mutex m_lock;
    static std::unordered_map<deg_t, std::shared_ptr<const Multiplication>> m_cache;

public:

    /*
     * Most of the multiplication classes take width as their only
     * parameter, which is really only there to help validate the
     * bases of the vector arguments.
     */

    static std::shared_ptr<const Multiplication> get(deg_t width)
    {
        std::lock_guard<std::mutex> access(m_lock);

        auto& found = m_cache[width];
        if (found) {
            return found;
        }

        found = std::make_shared<const Multiplication>(width);
    }

};


template <typename Multiplication>
std::mutex multiplication_registry<Multiplication>::m_lock;

template <typename Multiplication>
std::unordered_map<deg_t, std::shared_ptr<const Multiplication>>
    multiplication_registry<Multiplication>::m_cache;



}






#endif //LIBALGEBRA_LITE_REGISTRY_H
