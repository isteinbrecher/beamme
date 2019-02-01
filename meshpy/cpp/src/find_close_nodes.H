

#ifndef HEADER__
#define HEADER__

#include <iostream>

#include <stdexcept>

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <tuple>

namespace py = pybind11;


std::tuple<py::array_t<int>, int> find_close_nodes(
    py::array_t<double> coordinates_numpy, double eps);

#endif
