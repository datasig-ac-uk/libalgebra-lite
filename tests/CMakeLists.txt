cmake_minimum_required(VERSION 3.18)


include(GoogleTest)
find_package(GTest CONFIG REQUIRED)



add_library(lal_testing_common STATIC)

target_sources(lal_testing_common PRIVATE
        common/main.cpp
        )
target_link_libraries(lal_testing_common PRIVATE
        GTest::gtest_main
        )
