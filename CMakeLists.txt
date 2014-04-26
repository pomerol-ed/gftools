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
        set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-local-typedefs")
    if (NOT (GCC_VERSION VERSION_GREATER 4.6 OR GCC_VERSION VERSION_EQUAL 4.6))
        message(FATAL_ERROR "${PROJECT_NAME} requires g++ 4.6 or greater.")
    endif ()
elseif ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang" OR "${CMAKE_CXX_COMPILER_ID}" MATCHES "Intel")
    if ("${CMAKE_SYSTEM_NAME}" MATCHES "Linux")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libstdc++")
    elseif ("${CMAKE_SYSTEM_NAME}" MATCHES "Darwin")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
    else()
        message(FATAL_ERROR "Platform undefined")
    endif()
else ()
    message(FATAL_ERROR "Your C++ compiler does not support C++11.")
endif ()
#end :: C++11 check

# begin :: doxygen
set(DOXYFILE_SOURCE_DIR "./gftools")
include(UseDoxygen)
# end :: doxygen

# Eigen3 dependencies
find_package (Eigen3 3.1)
message(STATUS "Eigen3 includes: " ${EIGEN3_INCLUDE_DIR} )

# Boost
#set(Boost_USE_STATIC_LIBS        ON)
#set(Boost_USE_STATIC_RUNTIME    OFF)
find_package (Boost) 

include_directories(
    ${EIGEN3_INCLUDE_DIR}
    #${EIGEN3_INCLUDE_DIR}/unsupported
    ${Boost_INCLUDE_DIRS}
    ./gftools
)

# Build and enable tests
add_subdirectory(gftools)
option(Testing "Enable testing" OFF)
if (Testing)
    message(STATUS "Building tests")

    # Find gtest or otherwise fetch it into the build_dir/gtest
    if (DEFINED gtest_ROOT)
        message(STATUS "gtest source specified at ${gtest_ROOT}")
        find_path(gtest_ROOT NAMES "include/gtest/gtest.h" HINTS ${gtest_ROOT})
        if (NOT IS_DIRECTORY ${gtest_ROOT})
            message(WARNING "Provided wrong gtest_ROOT. Please unset gtest_ROOT - gtest will be fetched")
            unset(gtest_ROOT CACHE)
        endif()
    endif()

    if (NOT DEFINED gtest_ROOT)
        find_path(gtest_ROOT1 NAMES "include/gtest/gtest.h" HINTS ${CMAKE_SOURCE_DIR}/gtest-1.6.0  ${CMAKE_SOURCE_DIR}/gtest-1.7.0  ${CMAKE_SOURCE_DIR}/gtest   ${CMAKE_BINARY_DIR}/gtest)
        if (IS_DIRECTORY ${gtest_ROOT1})
            set (gtest_ROOT ${gtest_ROOT1})
        else()
            message(STATUS "Trying to fetch gtest via subversion")
            find_package(Subversion)
            execute_process(COMMAND "${Subversion_SVN_EXECUTABLE}" "checkout" "http://googletest.googlecode.com/svn/trunk/" "gtest" WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
            set (gtest_ROOT "${CMAKE_BINARY_DIR}/gtest")
        endif()
    endif()
    message(STATUS "googletest is in ${gtest_ROOT}")
    add_subdirectory(${gtest_ROOT})

    add_subdirectory(test)
    enable_testing()
endif (Testing)

# Build example
option(Examples "Build examples" OFF)
if (Examples)
    message(STATUS "Building examples")
    add_subdirectory(example)
endif (Examples)

# Generate pkg-config file
configure_file("${CMAKE_SOURCE_DIR}/gftools.pc.in" "${CMAKE_BINARY_DIR}/gftools.pc")
install(FILES "${CMAKE_BINARY_DIR}/gftools.pc" DESTINATION "lib/pkgconfig")
install(FILES "${CMAKE_SOURCE_DIR}/gftools.hpp" DESTINATION include)

