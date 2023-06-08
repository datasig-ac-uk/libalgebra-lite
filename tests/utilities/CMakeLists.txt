


add_executable(utilities_tests packed_integer.cpp)

target_link_libraries(utilities_tests PRIVATE
        GTest::gtest_main
        Libalgebra_lite
        )

include(GoogleTest)


if (WIN32)
    add_custom_command(TARGET utilities_tests POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy_if_different $<TARGET_RUNTIME_DLLS:utilities_tests> $<TARGET_FILE_DIR:utilities_tests>
            COMMAND_EXPAND_LISTS
            )
endif ()

gtest_discover_tests(utilities_tests)
