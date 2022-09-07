cmake_minimum_required(VERSION 3.22)
project(Libalgebra_lite)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

include(GenerateExportHeader)
include(GNUInstallDirs)



find_package(Boost CONFIG REQUIRED COMPONENTS headers)

add_library(Libalgebra_lite SHARED
        src/basis/hall_set.cpp
        src/basis/tensor_basis.cpp
        include/libalgebra_lite/free_tensor.h
        include/libalgebra_lite/lie.h
        include/libalgebra_lite/maps.h
        include/libalgebra_lite/vector_bundle.h
        src/algebra/lie_multiplier.cpp
        include/libalgebra_lite/polynomial_basis.h
        include/libalgebra_lite/polynomial.h
        src/algebra/polynomial_multiplier.cpp
        src/algebra/polynomial.cpp src/algebra/half_shuffle_multiplier.cpp src/algebra/shuffle_multiplier.cpp src/algebra/free_tensor_multiplier.cpp include/libalgebra_lite/registry.h src/basis/polynomial_basis.cpp src/maps.cpp)
generate_export_header(Libalgebra_lite)


target_include_directories(Libalgebra_lite PRIVATE
        ${CMAKE_CURRENT_LIST_DIR}
        ${CMAKE_CURRENT_BINARY_DIR})
target_include_directories(Libalgebra_lite PUBLIC
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/include/>"
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>"
        "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>")

target_link_libraries(Libalgebra_lite PRIVATE Boost::boost)





enable_testing()
add_subdirectory(tests)
