#ifndef GENERIC_CHAIN_ANNOTATOR_HEADER_H
#define GENERIC_CHAIN_ANNOTATOR_HEADER_H

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




namespace py = pybind11;



class ChainAnnotatorCpp : public AnnotatorBaseClassCpp {
    public:
        ChainAnnotatorCpp(std::string scheme = "imgt",
                std::string consensus_filepath = "");

        std::vector<std::tuple<std::vector<std::string>, double, std::string, std::string>>
                        analyze_seq(std::string sequence, double pidentity_threshold);


    protected:
        std::string scheme;

        std::unique_ptr<SingleChainAnnotatorCpp> analyzer;
        std::unique_ptr<CTermFinder> boundary_finder;

        bool split_subregion(const std::pair<int, int> &init_subregion,
                int &split_point, const std::string &sequence);
};

#endif
