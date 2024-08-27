#ifndef SINGLE_CHAIN_ANNOTATOR_HEADER_H
#define SINGLE_CHAIN_ANNOTATOR_HEADER_H

#include <pybind11/numpy.h>
#include <vector>
#include <string>
#include <tuple>
#include <thread>
#include <future>
#include <functional>
#include <filesystem>
#include <iostream>
#include <unordered_map>
#include "annotator_base_class.h"
#include "utilities.h"
#include "consensus_file_utilities.h"
#include "ig_aligner.h"
#include "numbering_constants.h"

#include <iostream>


namespace py = pybind11;

// We have 22 AAs -- the 20 standard aminos then a gap in query penalty
// vs gap in template -- thus, 22 expected values.
#define NUM_AAS 22

// Codes for the pathways that can link a score
// to the best-scoring parent.
#define LEFT_TRANSFER 0
#define DIAGONAL_TRANSFER 1
#define UP_TRANSFER 2

// The columns of the score matrix that are accessed for gap penalties.
#define QUERY_GAP_COLUMN 20
#define TEMPLATE_GAP_COLUMN 21





class SingleChainAnnotatorCpp : public AnnotatorBaseClassCpp {
    public:
        SingleChainAnnotatorCpp(std::vector<std::string> chains = {"H", "K", "L"},
                std::string scheme = "imgt", bool compress_init_gaps = false,
                bool multithread = false,
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
        bool multithread;

        std::vector<std::unique_ptr<IGAligner>> scoring_tools;


        std::tuple<std::vector<std::string>, double, std::string,
                std::string, std::vector<std::string>> analyze_test_only(std::string query_sequence,
                bool retrieve_cdr_labeling, py::array_t<double> scoreMatrix,
                py::array_t<uint8_t> pathTrace);

};

#endif
