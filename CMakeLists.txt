cmake_minimum_required (VERSION 2.8)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)
project (GFTools)

# Select default build type
IF(DEFINED CMAKE_BUILD_TYPE)
   SET(CMAKE_BUILD_TYPE ${CMAKE_BUILD_TYPE} CACHE STRING "Choose the type of build, options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug Release RelWithDebInfo MinSizeRel.")
ELSE()
   SET(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build, options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug Release RelWithDebInfo MinSizeRel.")
ENDIF()

# begin :: C++11 check
set(CMAKE_CXX_FLAGS "-Wall -std=c++11")
if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU")
    execute_process(
        COMMAND ${CMAKE_CXX_COMPILER} -dumpversion OUTPUT_VARIABLE GCC_VERSION)
    if (NOT (GCC_VERSION VERSION_GREATER 4.6 OR GCC_VERSION VERSION_EQUAL 4.6))
        message(FATAL_ERROR "${PROJECT_NAME} requires g++ 4.6 or greater.")
    endif ()
elseif ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
    if ("${CMAKE_SYSTEM_NAME}" MATCHES "Linux")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libstdc++")
    elseif ("${CMAKE_SYSTEM_NAME}" MATCHES "Darwin")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
    else()
        message(FATAL_ERROR "Platform undefined")
    endif()
elseif ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Intel")
#    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
else ()
    message(FATAL_ERROR "Your C++ compiler does not support C++11.")
endif ()
#end :: C++11 check

# begin :: doxygen
set(DOXYFILE_SOURCE_DIR "./include")
include(UseDoxygen)
# end :: doxygen

find_package (Eigen3 3.1)

message(STATUS "Eigen3 includes: " ${EIGEN3_INCLUDE_DIR} )

include_directories(
    ${EIGEN3_INCLUDE_DIR}
    ./include
)

add_subdirectory(include)
option(Testing "Enable testing" OFF)
if (Testing)
    message(STATUS "Building tests")
    add_subdirectory(test)
    enable_testing()
endif (Testing)

option(Examples "Build examples" OFF)
if (Examples)
    message(STATUS "Building examples")
    add_subdirectory(example)
endif (Examples)



