//
// Created by user on 23/07/22.
//

#ifndef LIBALGEBRA_LITE_IMPLEMENTATION_TYPES_H
#define LIBALGEBRA_LITE_IMPLEMENTATION_TYPES_H

#include <cstddef>
#include <cstdint>

#include <boost/predef.h>

#include "config.h"

#ifdef LAL_USE_LIBAGEBRA
#include <libalgebra/libalgebra.h>
#endif


namespace lal {

using dimn_t = std::size_t;
using idimn_t = std::ptrdiff_t;
using deg_t = std::int32_t;

using let_t = std::size_t;



} // namespace lal


#if BOOST_COMP_MSVC
#define LAL_INLINE_ALWAYS
#define LAL_INLINE_NEVER
#define LAL_RESTRICT
#define LAL_UNUSED
#elif BOOST_COMP_GNUC
#define LAL_INLINE_ALWAYS __attribute__((always_inline))
#define LAL_INLINE_NEVER __attribute__((noinline))
#define LAL_RESTRICT __restrict
#define LAL_UNUSED __attribute__((unused))
#else
#define LAL_INLINE_ALWAYS
#define LAL_INLINE_NEVER
#define LAL_RESTRICT
#define LAL_UNUSED
#endif


#if WIN32
#ifdef Libalgebra_Lite_EXPORTS
#define LAL_EXPORT_TEMPLATE(TMPL, ...) \
    extern template class TMPL<__VA_ARGS__>;
#else
#define LAL_EXPORT_TEMPLATE(TMPL, ...) \
    template class LIBALGEBRA_LITE_EXPORT TMPL<__VA_ARGS__>;
#endif
#else
#define LAL_EXPORT_TEMPLATE(TMPL, ...) \
    extern template class LIBALGEBRA_LITE_EXPORT TMPL<__VA_ARGS__>;
#endif

#endif //LIBALGEBRA_LITE_IMPLEMENTATION_TYPES_H
