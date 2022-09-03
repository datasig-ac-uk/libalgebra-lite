//
// Created by user on 23/07/22.
//

#ifndef LIBALGEBRA_LITE_IMPLEMENTATION_TYPES_H
#define LIBALGEBRA_LITE_IMPLEMENTATION_TYPES_H

#include <cstddef>
#include <cstdint>

#include <boost/predef.h>

namespace lal {

using dimn_t = std::size_t;
using idimn_t = std::ptrdiff_t;
using deg_t = std::int32_t;

using let_t = std::size_t;

} // namespace lal


#if BOOST_COMP_MSVC
#define LAL_INLINE_ALWAYS
#define LAL_INLINE_NEVER
#elif BOOST_COMP_GNUC
#define LAL_INLINE_ALWAYS __attribute__((always_inline))
#define LAL_INLINE_NEVER __attribute__((noinline))
#endif





#endif //LIBALGEBRA_LITE_IMPLEMENTATION_TYPES_H
