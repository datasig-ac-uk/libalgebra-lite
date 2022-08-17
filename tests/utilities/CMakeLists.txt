


add_executable(utilities_tests packed_integer.cpp)

target_link_libraries(utilities_tests PRIVATE
        GTest::gtest_main
        Libalgebra_lite
        )

include(GoogleTest)

gtest_discover_tests(utilities_tests)
