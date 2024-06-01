/* Contains the wrapper code for the C++ extension for alignment
 * calculations.
 */

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>    // Enables automatic type conversion for C++, python containers
#include <string>
#include "ig_aligner.h"
#include "vj_match_counter.h"
#include "responsibility_calcs.h"
#include "utilities.h"

namespace py = pybind11;
using namespace std;

PYBIND11_MODULE(antpack_cpp_ext, m){
    m.def("validate_sequence", &validate_sequence);
    m.def("sort_position_codes_cpp", &sort_position_codes_cpp);

    py::class_<IGAligner>(m, "IGAligner")
        .def(py::init<py::array_t<double>,
                std::vector<std::vector<std::string>>,
                std::string, std::string,
                double, double, bool>())
        .def("align", &IGAligner::align);

    py::class_<VJMatchCounter>(m, "VJMatchCounter")
        .def(py::init<std::vector<std::string>,
             std::vector<std::string>>() )
        .def("vjMatch", &VJMatchCounter::vjMatch)
        .def("findVJSequenceByName", &VJMatchCounter::findVJSequenceByName)
        .def("getSeqLists", &VJMatchCounter::getSeqLists);

    m.def("getProbsCExt", &getProbsCExt, py::call_guard<py::gil_scoped_release>());
    m.def("mask_terminal_deletions", &mask_terminal_deletions, py::call_guard<py::gil_scoped_release>());
    m.def("getProbsCExt_masked", &getProbsCExt_masked, py::call_guard<py::gil_scoped_release>());
}
