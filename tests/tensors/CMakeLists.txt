

add_executable(test_tensors
        tensor_fixture.h
        shuffle_tensor_multiplier.cpp
        half_shuffle_multiplier.cpp
        tensor_multiplication.cpp
        tensor_basis.cpp
        dense_tensor.cpp
        sparse_tensor.cpp)


target_link_libraries(test_tensors PRIVATE GTest::gtest_main Libalgebra_lite)


target_sources(test_tensors PRIVATE
        free_tensor_multiplier.cpp
        )


if (WIN32)
    add_custom_command(TARGET test_tensors POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy_if_different $<TARGET_RUNTIME_DLLS:test_tensors> $<TARGET_FILE_DIR:test_tensors>
            COMMAND_EXPAND_LISTS
            )
endif ()

gtest_discover_tests(test_tensors)
