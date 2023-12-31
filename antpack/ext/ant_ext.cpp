/* Contains the wrapper code for the C++ extension for alignment
 * calculations.
 */

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>    // Enables automatic type conversion for C++, python containers
#include <string>
#include "aligner.h"
#include "utilities.h"

namespace py = pybind11;
using namespace std;

PYBIND11_MODULE(ant_ext, m){
    m.def("validate_sequence", &validate_sequence);

    py::class_<BasicAligner>(m, "BasicAligner")
        .def(py::init<py::array_t<double>,
                std::vector<std::vector<std::string>>,
                std::string, std::string>())
        .def("align", &BasicAligner::align);
}