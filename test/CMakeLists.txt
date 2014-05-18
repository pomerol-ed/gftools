include_directories (. ../include ${gtest_ROOT}/include)

enable_testing()

# Here all the tests are set. The source is file is assumed to be ${test}.cpp
#add_library(gftools_extern STATIC externs.hpp externs.cpp)

set (tests 
num_io_test
tuple_tools_test
grid_test
kmesh_test
matsubara_grid_test
real_grid_test
enum_grid_test
grid_tools_test
eval_expression_test
container_test
grid_object_test
#SaveLoadTest
#ShiftTest
#InterpolateTest
)

foreach (test ${tests})
    set(test_src ${test}.cpp)
    add_executable(${test} ${test_src} $<TARGET_OBJECTS:libgftools>)
    target_link_libraries(${test} gtest_main)
    add_test(${test} ${test})
endforeach(test)
