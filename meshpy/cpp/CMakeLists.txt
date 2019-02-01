cmake_minimum_required(VERSION 3.0)

# Project name.
project(meshpy_cpp)

# Find python library and set the include directories.
find_package(PythonLibs REQUIRED)
include_directories( ${PYTHON_INCLUDE_DIRS} )

# Path to pybind11.
add_subdirectory(lib/pybind11)
include_directories( lib/pybind11/include )

# Compiler flags.
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -Wall -Wextra -O3")

# Add libraries.
add_library(meshpy_cpp MODULE src/meshpy_cpp.cpp)

# Define target.
set_target_properties(
  meshpy_cpp
  PROPERTIES
  PREFIX "${PYTHON_MODULE_PREFIX}"
  SUFFIX "${PYTHON_MODULE_EXTENSION}"
)
