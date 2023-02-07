

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




gtest_discover_tests(test_tensors)
