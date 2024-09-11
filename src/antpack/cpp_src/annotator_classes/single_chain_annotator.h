#ifndef SINGLE_CHAIN_ANNOTATOR_HEADER_H
#define SINGLE_CHAIN_ANNOTATOR_HEADER_H

#include <pybind11/numpy.h>
#include <vector>
#include <string>
#include <tuple>
#include <filesystem>
#include <unordered_map>
#include "annotator_base_class.h"
#include "cterm_finder.h"
#include "../utilities/utilities.h"
#include "../utilities/consensus_file_utilities.h"
#include "ig_aligner.h"
#include "../numbering_constants.h"



namespace py = pybind11;




class SingleChainAnnotatorCpp : public AnnotatorBaseClassCpp {
    public:
        SingleChainAnnotatorCpp(std::vector<std::string> chains = {"H", "K", "L"},
                std::string scheme = "imgt", bool compress_init_gaps = false,
                std::string consensus_filepath = "");

        std::tuple<std::vector<std::string>, double, std::string,
            std::string> analyze_seq(std::string);
        std::vector<std::tuple<std::vector<std::string>, double, std::string,
            std::string>> analyze_seqs(std::vector<std::string> sequences);
        std::unordered_map<std::string, std::vector<int>> get_cdr_breakpoints();

    protected:
        std::vector<std::string> chains;
        std::string scheme;
        bool compress_init_gaps;

        std::vector<std::unique_ptr<IGAligner>> scoring_tools;
        std::unique_ptr<CTermFinder> boundary_finder;

        int align_input_subregion(std::tuple<std::vector<std::string>, double,
                std::string, std::string> &best_result, double &best_identity,
                std::string &query_sequence);

        std::tuple<std::vector<std::string>, double, std::string,
                std::string, std::vector<std::string>> analyze_test_only(std::string query_sequence,
                bool retrieve_cdr_labeling, py::array_t<double> scoreMatrix,
                py::array_t<uint8_t> pathTrace);

};

#endif
