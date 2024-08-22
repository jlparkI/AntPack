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





class SingleChainAnnotatorCpp {
    public:
        SingleChainAnnotatorCpp(std::vector<std::string> chains = {"H", "K", "L"},
                std::string scheme = "imgt", bool compress_init_gaps = false,
                bool multithread = true, std::string project_filepath = "");
        std::tuple<std::vector<std::string>, double, std::string,
            std::string> analyze_seq(std::string);
        std::vector<std::tuple<std::vector<std::string>, double, std::string,
            std::string>> analyze_seqs(std::vector<std::string> sequences);

        std::vector<std::string> sort_position_codes(std::vector<std::string> position_code_list);
        std::tuple<std::vector<std::string>, std::vector<std::string>> build_msa(std::vector<std::string> sequences,
            std::vector<std::tuple<std::vector<std::string>, double, std::string, std::string>> annotations);
        std::tuple<std::string, std::vector<std::string>, int, int> trim_alignment(std::string sequence,
            std::tuple<std::vector<std::string>, double, std::string, std::string> alignment);
        std::vector<std::string> assign_cdr_labels(std::tuple<std::vector<std::string>, 
               double, std::string, std::string> alignment, std::string cdr_scheme = "");

    protected:
        std::vector<std::string> chains;
        std::string scheme;
        bool compress_init_gaps;
        bool multithread;
        std::unordered_map<std::string, std::vector<int>> cdr_breakpoints;

        std::vector<std::unique_ptr<IGAligner>> scoring_tools;

        std::array<std::string, 7> cdr_region_labels {{"fmwk1", "cdr1", "fmwk2", "cdr2",
                                "fmwk3", "cdr3", "fmwk4"}};

        std::tuple<std::vector<std::string>, double, std::string,
                std::string, std::vector<std::string>> analyze_test_only(std::string query_sequence,
                bool retrieve_cdr_labeling, py::array_t<double> scoreMatrix,
                py::array_t<uint8_t> pathTrace);

};

#endif
