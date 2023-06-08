

add_executable(test_polynomials
        polynomial_fixture.h
        polynomial_basis.cpp)

target_link_libraries(test_polynomials PRIVATE
        GTest::gtest_main
        Libalgebra_lite)


if (WIN32)
    add_custom_command(TARGET test_polynomials POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy_if_different $<TARGET_RUNTIME_DLLS:test_polynomials> $<TARGET_FILE_DIR:test_polynomials>
            COMMAND_EXPAND_LISTS
            )
endif ()

gtest_discover_tests(test_polynomials)
