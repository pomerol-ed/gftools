### gftools cmake configuration

cmake_minimum_required (VERSION 2.8)

# gftools use only CXX compiler 
project (GFTools CXX)

# Disable build in source
if (${CMAKE_BINARY_DIR} STREQUAL ${CMAKE_SOURCE_DIR})
    message(FATAL_ERROR "In source builds are disabled. Please use a separate build directory.")
endif()
set(CMAKE_DISABLE_SOURCE_CHANGES ON)
set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)
 # find additional cmake module in cmake dir
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

# Print build type
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

# set c++-11 flags
option(AutoSetCXX11 "Auto set c++11 flags" ON)
if (AutoSetCXX11)
    include(SetCXX11) 
endif (AutoSetCXX11)

# doxygen
set(DOXYFILE_SOURCE_DIR "./gftools")
include(UseDoxygen)

### Dependencies
# Eigen3 
find_package (Eigen3 3.1)
message(STATUS "Eigen3 includes: " ${EIGEN3_INCLUDE_DIR} )
# Boost
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
    include(EnableGtests)
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
configure_file("${CMAKE_SOURCE_DIR}/gftools.hpp" "${CMAKE_BINARY_DIR}/gftools.hpp")
install(FILES "${CMAKE_BINARY_DIR}/gftools.pc" DESTINATION "lib/pkgconfig")
install(FILES "${CMAKE_BINARY_DIR}/gftools.hpp" DESTINATION include)

