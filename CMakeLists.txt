### gftools cmake configuration

cmake_minimum_required (VERSION 3.1)

# gftools use only CXX compiler
project (GFTools CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

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

add_library(gftools INTERFACE)
target_include_directories(gftools INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:include>
)

install(TARGETS gftools EXPORT gftools-config)
install(DIRECTORY gftools DESTINATION ${CMAKE_INSTALL_PREFIX}/include)
install(EXPORT gftools-config DESTINATION share/gftools/cmake)

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
