

add_executable(test_polynomials
        polynomial_fixture.h
        polynomial_basis.cpp)

target_link_libraries(test_polynomials PRIVATE
        GTest::gtest_main
        Libalgebra_lite)


gtest_discover_tests(test_polynomials)
