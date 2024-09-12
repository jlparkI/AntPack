/* Contains the wrapper code for the C++ extension for alignment
 * calculations.
 */

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/vector.h>
#include <nanobind/stl/map.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/tuple.h>
#include <nanobind/stl/pair.h>

#include <string>
#include "annotator_classes/single_chain_annotator.h"
#include "annotator_classes/paired_chain_annotator.h"
#include "annotator_classes/annotator_base_class.h"
#include "annotator_classes/ig_aligner.h"
#include "annotator_classes/cterm_finder.h"
#include "vj_assignment/vj_match_counter.h"
#include "humanness_calcs/responsibility_calcs.h"
#include "utilities/utilities.h"

namespace nb = nanobind;
using namespace std;

NB_MODULE(antpack_cpp_ext, m){
    m.def("validate_sequence", &validate_sequence);

    nb::class_<AnnotatorBaseClassCpp>(m, "AnnotatorBaseClassCpp")
        .def(nb::init<std::string, std::string>())
        .def("sort_position_codes", &AnnotatorBaseClassCpp::sort_position_codes)
        .def("build_msa", &AnnotatorBaseClassCpp::build_msa)
        .def("assign_cdr_labels", &AnnotatorBaseClassCpp::assign_cdr_labels)
        .def("trim_alignment", &AnnotatorBaseClassCpp::trim_alignment);

    nb::class_<SingleChainAnnotatorCpp, AnnotatorBaseClassCpp>(m, "SingleChainAnnotatorCpp")
        .def(nb::init<std::vector<std::string>,
                std::string, bool, std::string>())
        .def("analyze_seq", &SingleChainAnnotatorCpp::analyze_seq)
        .def("analyze_seqs", &SingleChainAnnotatorCpp::analyze_seqs);
        //.def("_test_needle_scoring", &SingleChainAnnotatorCpp::_test_needle_scoring);

    nb::class_<PairedChainAnnotatorCpp, AnnotatorBaseClassCpp>(m, "PairedChainAnnotatorCpp")
        .def(nb::init<std::string, std::string>())
        .def("analyze_seq", &PairedChainAnnotatorCpp::analyze_seq);

    nb::class_<IGAligner>(m, "IGAligner")
        .def(nb::init<std::string,
                std::string, std::string,
                double, double, bool>());

    nb::class_<VJMatchCounter>(m, "VJMatchCounter")
        .def(nb::init<std::map<std::string, std::vector<std::string>>,
                std::map<std::string, std::vector<std::string>>,
                nb::ndarray<double, nb::shape<22,22>, nb::device::cpu, nb::c_contig>,
                std::string>() )
        .def("assign_vj_genes", &VJMatchCounter::assign_vj_genes)
        .def("get_vj_gene_sequence", &VJMatchCounter::get_vj_gene_sequence)
        .def("get_seq_lists", &VJMatchCounter::get_seq_lists);

    m.def("getProbsCExt", &getProbsCExt, nb::call_guard<nb::gil_scoped_release>());
    m.def("mask_terminal_deletions", &mask_terminal_deletions, nb::call_guard<nb::gil_scoped_release>());
    m.def("getProbsCExt_masked", &getProbsCExt_masked, nb::call_guard<nb::gil_scoped_release>());
}
