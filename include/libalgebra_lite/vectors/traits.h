//
// Created by user on 27/07/22.
//

#ifndef LIBALGEBRA_LITE_VECTOR_TRAITS_H
#define LIBALGEBRA_LITE_VECTOR_TRAITS_H

namespace alg {
namespace dtl {

template <typename Vector>
struct get_view_type_impl;

template <typename VectorView>
struct get_vector_type_impl;
}

template <template <typename, typename, typename...> class Vector>
struct associated_view_template;

template <template <typename, typename, typename...> class VectorView>
struct associated_vector_template;



template <typename Vector>
using associated_view_type = dtl::get_view_type_impl<Vector>;

template <typename VectorView>
using associated_vec_type = dtl::get_vector_type_impl<VectorView>;

namespace dtl {
template<typename Vector>
struct vector_base_trait;

template <typename Vector>
class vector_type_of_impl
{
    using owned_t = typename Vector::owned_vector_type;
    using view_t = typename Vector::view_vector_type;

    static_assert(std::is_base_of<view_t, owned_t>::value, "mismatch in view and owned vector types");

    template <typename Arg>
    static owned_t test(Arg&);

    template <typename>
    static view_t test(...);


public:
    using type = decltype(test<owned_t>(std::declval<Vector&>()));
};




} // namespace dtl

template <typename Vector>
using vector_base = typename dtl::vector_base_trait<Vector>::type;

template <typename Vector>
using vector_type_of = typename dtl::vector_type_of_impl<Vector>;


}

#endif //LIBALGEBRA_LITE_VECTOR_TRAITS_H
