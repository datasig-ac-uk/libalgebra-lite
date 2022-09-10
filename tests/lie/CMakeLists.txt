

add_executable(test_lie hall_basis.cpp lie_fixture.h)


target_link_libraries(test_lie PRIVATE
        GTest::gtest_main
        Libalgebra_lite
        )

target_sources(test_lie PRIVATE
        lie_multiplier.cpp
        )

gtest_discover_tests(test_lie)
