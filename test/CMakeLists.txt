include_directories (. ../include)

enable_testing()

# Here all the tests are set. The source is file is assumed to be ${test}.cpp
set (tests 
GridTest
ContainerTest
GridObjectTest
ShiftTest
InterpolateTest
)

foreach (test ${tests})
    set(test_src ${test}.cpp)
    add_executable(${test} ${test_src} $<TARGET_OBJECTS:libgftools>)
    add_test(${test} ${test})
endforeach(test)
