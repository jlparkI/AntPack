#ifndef SINGLE_CHAIN_ANNOTATOR_HEADER_H
#define SINGLE_CHAIN_ANNOTATOR_HEADER_H

#include <nanobind/ndarray.h>
#include <vector>
#include <string>
#include <tuple>
#include <unordered_map>
#include "annotator_base_class.h"
#include "../utilities/utilities.h"
#include "ig_aligner.h"
#include "../numbering_constants.h"



namespace nb = nanobind;




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

        void _test_needle_scoring(std::string query_sequence,
                    nb::ndarray<double, nb::shape<-1,-1>, nb::device::cpu, nb::c_contig> scoreMatrix,
                    nb::ndarray<uint8_t, nb::shape<-1,-1>, nb::device::cpu, nb::c_contig> pathTraceMat,
                    std::string chain);
        

    protected:
        std::vector<std::string> chains;
        std::string scheme;
        bool compress_init_gaps;

        std::vector<IGAligner> scoring_tools;

        int align_input_subregion(std::tuple<std::vector<std::string>, double,
                std::string, std::string> &best_result, double &best_identity,
                std::string &query_sequence);

};

#endif
