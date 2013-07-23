include_directories (. ../include)

enable_testing()

# Here all the tests are set. The source is file is assumed to be ${test}.cpp
#add_library(gftools_extern STATIC externs.hpp externs.cpp)

set (tests 
GridTest
ContainerTest
GridObjectTest
SaveLoadTest
ShiftTest
InterpolateTest
)

foreach (test ${tests})
    set(test_src ${test}.cpp)
    add_executable(${test} ${test_src} $<TARGET_OBJECTS:libgftools>)
#    target_link_libraries(${test} gftools_extern)
    add_test(${test} ${test})
endforeach(test)
