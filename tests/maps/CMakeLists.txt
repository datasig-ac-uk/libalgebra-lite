

add_executable(test_maps maps_fixture.h maps.cpp)


target_link_libraries(test_maps PRIVATE Libalgebra_lite GTest::gtest_main)

if (WIN32)
    add_custom_command(TARGET test_maps POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy_if_different $<TARGET_RUNTIME_DLLS:test_maps> $<TARGET_FILE_DIR:test_maps>
            COMMAND_EXPAND_LISTS
            )
endif ()
gtest_discover_tests(test_maps)
