

add_executable(test_lie)


target_link_libraries(test_lie PRIVATE
        GTest::gtest_main
        Libalgebra_lite
        )

target_sources(test_lie PRIVATE
        lie_multiplier.cpp
        )
