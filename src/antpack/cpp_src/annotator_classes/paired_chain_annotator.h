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
#include "annotator_base_class.h"
#include "single_chain_annotator.h"
#include "../numbering_constants.h"


#define CTERMINAL_STANDARD_OFFSET 9


namespace py = pybind11;



class PairedChainAnnotatorCpp : public AnnotatorBaseClassCpp {
    public:
        PairedChainAnnotatorCpp(std::string scheme = "imgt",
                bool multithread = false,
                std::string consensus_filepath = "");

        std::pair<std::tuple<std::vector<std::string>, double, std::string, std::string>,
            std::tuple<std::vector<std::string>, double, std::string, std::string>>
            analyze_seq(std::string sequence);


    protected:

        std::string scheme;
        bool multithread;

        std::unique_ptr<SingleChainAnnotatorCpp> light_chain_analyzer;
        std::unique_ptr<SingleChainAnnotatorCpp> heavy_chain_analyzer;
        std::unique_ptr<SingleChainAnnotatorCpp> analyzer;

        std::unique_ptr<CTermFinder> boundary_finder;
};

#endif