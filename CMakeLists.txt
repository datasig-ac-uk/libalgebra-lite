cmake_minimum_required(VERSION 3.22)
project(Libalgebra_lite)

#set(CMAKE_CXX_STANDARD 14)
#set(CMAKE_CXX_STANDARD_REQUIRED ON)
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

include(GenerateExportHeader)
include(GNUInstallDirs)

option(LIBALGEBRA_LITE_BUILD_TESTS OFF "Build tests for libalgebra lite")

find_package(Boost REQUIRED)

add_library(Libalgebra_lite SHARED
        src/algebra/lie_multiplier.cpp
        src/algebra/polynomial_multiplier.cpp
        src/algebra/polynomial.cpp
        src/algebra/half_shuffle_multiplier.cpp
        src/algebra/free_tensor_multiplier.cpp
        src/algebra/shuffle_multiplier.cpp
        src/algebra/polynomial_ring.cpp
        src/basis/hall_set.cpp
        src/basis/tensor_basis.cpp
        src/basis/monomial.cpp
        src/basis/polynomial_basis.cpp
        src/basis/unpacked_tensor_word.cpp
        src/coefficients/floating_fields.cpp
        src/coefficients/rational_field.cpp
        src/maps.cpp
        include/libalgebra_lite/detail/integer_maths.h
        include/libalgebra_lite/detail/macros.h
        include/libalgebra_lite/detail/notnull.h
        include/libalgebra_lite/detail/traits.h
        include/libalgebra_lite/algebra.h
        include/libalgebra_lite/basis.h
        include/libalgebra_lite/basis_traits.h
        include/libalgebra_lite/coefficients.h
        include/libalgebra_lite/dense_vector.h
        include/libalgebra_lite/free_tensor.h
        include/libalgebra_lite/hall_set.h
        include/libalgebra_lite/implementation_types.h
        include/libalgebra_lite/index_key.h
        include/libalgebra_lite/key_range.h
        include/libalgebra_lite/lie.h
        include/libalgebra_lite/maps.h
        include/libalgebra_lite/operators.h
        include/libalgebra_lite/packed_integer.h
        include/libalgebra_lite/polynomial.h
        include/libalgebra_lite/polynomial_basis.h
        include/libalgebra_lite/registry.h
        include/libalgebra_lite/shuffle_tensor.h
        include/libalgebra_lite/sparse_vector.h
        include/libalgebra_lite/tensor_basis.h
        include/libalgebra_lite/unpacked_tensor_word.h
        include/libalgebra_lite/vector.h
        include/libalgebra_lite/vector_base.h
        include/libalgebra_lite/vector_bundle.h
        include/libalgebra_lite/vector_traits.h
        )

add_library(Libalgebra_lite::Libalgebra_lite ALIAS Libalgebra_lite)
generate_export_header(Libalgebra_lite)

find_package(Bignum QUIET)
if (NOT TARGET Bignum::Bignum)
    set(LAL_NO_USE_GMP ON)
else ()
    target_link_libraries(Libalgebra_lite PUBLIC Bignum::Bignum)
    set(LAL_NO_USE_GMP OFF)
endif ()

target_include_directories(Libalgebra_lite PRIVATE
        ${CMAKE_CURRENT_LIST_DIR}
        ${CMAKE_CURRENT_BINARY_DIR})
target_include_directories(Libalgebra_lite PUBLIC
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/include/>"
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>"
        "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>")

target_link_libraries(Libalgebra_lite PUBLIC Boost::boost)

#if (TARGET Libalgebra::Libalgebra)
#    set(LAL_USE_LIBALGEBRA ON)
#    target_link_libraries(Libalgebra_lite PUBLIC Libalgebra::Libalgebra)
#endif ()

set_target_properties(Libalgebra_lite PROPERTIES
        INSTALL_RPATH $ORIGIN)


if (NOT DEFINED LAL_GMP_ALTERNATIVE_HEADER)
    set(LAL_GMP_ALTERNATIVE_HEADER cpp_int.hpp)
    set(LAL_GMP_ALTERNATIVE_TYPE cpp_rational)
endif ()

configure_file(cmake/config.h.in config.h NO_SOURCE_PERMISSIONS @ONLY)
configure_file(cmake/rationals.h.in rationals.h NO_SOURCE_PERMISSIONS @ONLY)


install(TARGETS Libalgebra_lite
        EXPORT Libalgebra_liteTargets
        PUBLIC_HEADER DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/libalgebra_lite"
        RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
        LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
        ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
        )
install(DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/include/libalgebra_lite
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
        COMPONENT Development
        EXCLUDE_FROM_ALL
        )
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/config.h
        ${CMAKE_CURRENT_BINARY_DIR}/rationals.h
        ${CMAKE_CURRENT_BINARY_DIR}/libalgebra_lite_export.h
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/libalgebra_lite
        COMPONENT Development
        EXCLUDE_FROM_ALL
        )
install(EXPORT Libalgebra_liteTargets
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/Libalgebra_lite
        NAMESPACE Libalgebra_lite::
        COMPONENT Development
        EXCLUDE_FROM_ALL
        )

include(CMakePackageConfigHelpers)
write_basic_package_version_file(Libalgebra_liteVersion.cmake
        VERSION 1.0.0
        COMPATIBILITY SameMajorVersion
        )
set(cmakeModulesDir cmake)
configure_package_config_file(${CMAKE_CURRENT_LIST_DIR}/cmake/Libalgebra_liteConfig.cmake.in
        Libalgebra_liteConfig.cmake
        INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/Libalgebra_lite
        PATH_VARS cmakeModulesDir
        )

install(FILES
        ${CMAKE_CURRENT_BINARY_DIR}/Libalgebra_liteConfig.cmake
        ${CMAKE_CURRENT_BINARY_DIR}/Libalgebra_liteVersion.cmake
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/Libalgebra_lite
        COMPONENT Development
        EXCLUDE_FROM_ALL
        )

if(CMAKE_IS_TOP_LEVEL OR LIBALGEBRA_LITE_BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()
