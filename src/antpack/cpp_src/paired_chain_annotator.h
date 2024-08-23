#ifndef PAIRED_CHAIN_ANNOTATOR_HEADER_H
#define PAIRED_CHAIN_ANNOTATOR_HEADER_H

#include <pybind11/numpy.h>
#include <vector>
#include <string>
#include <algorithm>
#include <tuple>
#include <set>
#include <array>
#include "cterm_finder.h"
#include "single_chain_annotator.h"
#include "numbering_constants.h"



namespace py = pybind11;



class PairedChainAnnotatorCpp {
    public:
        PairedChainAnnotatorCpp(std::string scheme = "imgt",
                bool multithread = false,
                std::string project_filepath = "");

        std::tuple<std::vector<std::string>, double, std::string,
            std::string> analyze_seq(std::string);

        std::vector<std::string> sort_position_codes(std::vector<std::string> position_code_list);
        std::tuple<std::vector<std::string>, std::vector<std::string>> build_msa(std::vector<std::string> sequences,
            std::vector<std::tuple<std::vector<std::string>, double, std::string, std::string>> annotations);
        std::tuple<std::string, std::vector<std::string>, int, int> trim_alignment(std::string sequence,
            std::tuple<std::vector<std::string>, double, std::string, std::string> alignment);
        std::vector<std::string> assign_cdr_labels(std::tuple<std::vector<std::string>, 
               double, std::string, std::string> alignment, std::string cdr_scheme = "");


    protected:

        std::string scheme;
        bool multithread;
        std::string project_filepath;

        std::unique_ptr<SingleChainAnnotatorCpp> light_chain_analyzer;
        std::unique_ptr<SingleChainAnnotatorCpp> heavy_chain_analyzer;
        std::unique_ptr<SingleChainAnnotatorCpp> analyzer;

        std::unique_ptr<CTermFinder> boundary_finder;
};

#endif
