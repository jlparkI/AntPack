#ifndef PAIRED_CHAIN_ANNOTATOR_HEADER_H
#define PAIRED_CHAIN_ANNOTATOR_HEADER_H

// C++ headers
#include <vector>
#include <string>
#include <algorithm>
#include <tuple>
#include <set>
#include <array>

// Library headers


// Project headers
#include "cterm_finder.h"
#include "annotator_base_class.h"
#include "single_chain_annotator.h"
#include "numbering_constants.inc"



namespace NumberingTools {


class PairedChainAnnotatorCpp : public AnnotatorBaseClassCpp {
 public:
        PairedChainAnnotatorCpp(std::string scheme = "imgt",
                std::string consensus_filepath = "");

        std::pair<std::tuple<std::vector<std::string>, double, std::string, std::string>,
            std::tuple<std::vector<std::string>, double, std::string, std::string>>
            analyze_seq(std::string sequence);
        std::tuple<std::vector<std::tuple<std::vector<std::string>, double, std::string, std::string>>,
            std::vector<std::tuple<std::vector<std::string>, double, std::string, std::string>>>
            analyze_seqs(std::vector<std::string> sequences);


 protected:
        std::string scheme;

        std::unique_ptr<NumberingTools::SingleChainAnnotatorCpp> light_chain_analyzer;
        std::unique_ptr<NumberingTools::SingleChainAnnotatorCpp> heavy_chain_analyzer;
        std::unique_ptr<NumberingTools::SingleChainAnnotatorCpp> analyzer;
};

}  // namespace NumberingTools

#endif
