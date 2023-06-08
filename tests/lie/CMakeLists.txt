

add_executable(test_lie hall_basis.cpp lie_fixture.h)


target_link_libraries(test_lie PRIVATE
        GTest::gtest_main
        Libalgebra_lite
        )

target_sources(test_lie PRIVATE
        lie_multiplier.cpp
        )

if (WIN32)
    add_custom_command(TARGET test_lie POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy_if_different $<TARGET_RUNTIME_DLLS:test_lie> $<TARGET_FILE_DIR:test_lie>
            COMMAND_EXPAND_LISTS
            )
endif()

gtest_discover_tests(test_lie)
