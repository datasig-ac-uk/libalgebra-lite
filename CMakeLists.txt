cmake_minimum_required(VERSION 3.22)
project(Libalgebra_lite)

set(CMAKE_CXX_STANDARD 14)

add_library(Libalgebra_lite STATIC
        library.cpp
        libalgebra_lite/vectors/dense_vector.h
        libalgebra_lite/implementation_types.h)


target_include_directories(Libalgebra_lite PRIVATE
        ${CMAKE_CURRENT_LIST_DIR})
