

add_executable(test_maps maps_fixture.h maps.cpp)


target_link_libraries(test_maps PRIVATE Libalgebra_lite GTest::gtest_main)

gtest_discover_tests(test_maps)
