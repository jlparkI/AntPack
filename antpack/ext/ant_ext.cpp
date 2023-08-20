/* Contains the wrapper code for the C++ extension for alignment
 * calculations.
 */

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>    // Enables automatic type conversion for C++, python containers
#include <string>
#include "aligners.h"

namespace py = pybind11;
using namespace std;

PYBIND11_MODULE(ant_ext, m){
    py::class_<IMGTAligner>(m, "IMGTAligner")
        .def(py::init<py::array_t<double>>())
        .def("align", &IMGTAligner::align);
}