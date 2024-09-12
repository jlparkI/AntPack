/* Contains the wrapper code for the C++ extension for alignment
 * calculations.
 */

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>    // Enables automatic type conversion for C++, python containers
#include <string>
#include "annotator_classes/single_chain_annotator.h"
#include "annotator_classes/paired_chain_annotator.h"
#include "annotator_classes/annotator_base_class.h"
#include "annotator_classes/ig_aligner.h"
#include "annotator_classes/cterm_finder.h"
#include "vj_assignment/vj_match_counter.h"
#include "humanness_calcs/responsibility_calcs.h"
#include "utilities/utilities.h"

namespace py = pybind11;
using namespace std;

PYBIND11_MODULE(antpack_cpp_ext, m){
    m.def("validate_sequence", &validate_sequence);

    py::class_<AnnotatorBaseClassCpp>(m, "AnnotatorBaseClassCpp")
        .def(py::init<std::string, std::string>())
        .def("sort_position_codes", &AnnotatorBaseClassCpp::sort_position_codes)
        .def("build_msa", &AnnotatorBaseClassCpp::build_msa)
        .def("assign_cdr_labels", &AnnotatorBaseClassCpp::assign_cdr_labels)
        .def("trim_alignment", &AnnotatorBaseClassCpp::trim_alignment);

    py::class_<SingleChainAnnotatorCpp, AnnotatorBaseClassCpp>(m, "SingleChainAnnotatorCpp")
        .def(py::init<std::vector<std::string>,
                std::string, bool, std::string>())
        .def("analyze_seq", &SingleChainAnnotatorCpp::analyze_seq)
        .def("analyze_seqs", &SingleChainAnnotatorCpp::analyze_seqs);
        //.def("_test_needle_scoring", &SingleChainAnnotatorCpp::_test_needle_scoring);

    py::class_<PairedChainAnnotatorCpp, AnnotatorBaseClassCpp>(m, "PairedChainAnnotatorCpp")
        .def(py::init<std::string, std::string>())
        .def("analyze_seq", &PairedChainAnnotatorCpp::analyze_seq);

    py::class_<IGAligner>(m, "IGAligner")
        .def(py::init<py::array_t<double>,
                std::vector<std::vector<std::string>>,
                std::string, std::string,
                double, double, bool>());

    py::class_<VJMatchCounter>(m, "VJMatchCounter")
        .def(py::init<std::map<std::string, std::vector<std::string>>,
                std::map<std::string, std::vector<std::string>>,
                py::array_t<double, py::array::c_style>,
                std::string>() )
        .def("assign_vj_genes", &VJMatchCounter::assign_vj_genes)
        .def("get_vj_gene_sequence", &VJMatchCounter::get_vj_gene_sequence)
        .def("get_seq_lists", &VJMatchCounter::get_seq_lists);

    m.def("getProbsCExt", &getProbsCExt, py::call_guard<py::gil_scoped_release>());
    m.def("mask_terminal_deletions", &mask_terminal_deletions, py::call_guard<py::gil_scoped_release>());
    m.def("getProbsCExt_masked", &getProbsCExt_masked, py::call_guard<py::gil_scoped_release>());
}
