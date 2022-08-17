cmake_minimum_required(VERSION 3.22)
project(Libalgebra_lite)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

include(GenerateExportHeader)
include(GNUInstallDirs)



find_package(Boost CONFIG REQUIRED COMPONENTS headers)

add_library(Libalgebra_lite SHARED
        library.cpp
        src/basis/hall_set.cpp
        src/basis/tensor_basis.cpp
        include/libalgebra_lite/free_tensor.h
        include/libalgebra_lite/lie.h
        include/libalgebra_lite/maps.h
        )
generate_export_header(Libalgebra_lite)


target_include_directories(Libalgebra_lite PRIVATE
        ${CMAKE_CURRENT_LIST_DIR}
        ${CMAKE_CURRENT_BINARY_DIR})
target_include_directories(Libalgebra_lite PUBLIC
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/include/>"
        "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>")

target_link_libraries(Libalgebra_lite PRIVATE Boost::boost)




add_executable(playground playground.cpp)
target_include_directories(playground PRIVATE ${CMAKE_CURRENT_LIST_DIR}/include)
target_link_libraries(playground PRIVATE Boost::boost)



enable_testing()
add_subdirectory(tests)
