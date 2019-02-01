#include <iostream>

#include <stdexcept>

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <tuple>

#include "find_close_nodes.H"

PYBIND11_MODULE(meshpy_cpp, m)
{
  m.doc() = "Add two vectors using pybind11";  // optional module docstring

  m.def("find_close_nodes", &find_close_nodes, "Add two NumPy arrays");
}
