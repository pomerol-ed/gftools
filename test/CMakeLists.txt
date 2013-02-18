include_directories (. ../include)

enable_testing()

set (loglib_src
    Logger.h Logger.cpp
    )
add_library(loglib STATIC ${loglib_src})


# Here all the tests are set. The source is file is assumed to be ${test}.cpp
set (tests 
GridTest
ContainerTest
GridObjectTest
ShiftTest
)

foreach (test ${tests})
    set(test_src ${test}.cpp)
    add_executable(${test} ${test_src} $<TARGET_OBJECTS:libgftools>)
    target_link_libraries(${test} loglib)
    add_test(${test} ${test})
endforeach(test)
