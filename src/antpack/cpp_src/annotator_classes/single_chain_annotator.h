#ifndef SINGLE_CHAIN_ANNOTATOR_HEADER_H
#define SINGLE_CHAIN_ANNOTATOR_HEADER_H

// C++ headers
#include <vector>
#include <string>
#include <tuple>
#include <unordered_map>

// Library headers
#include <nanobind/ndarray.h>

// Project headers
#include "annotator_base_class.h"
#include "../utilities/utilities.h"
#include "ig_aligner.h"
#include "numbering_constants.inc"



namespace nb = nanobind;


namespace NumberingTools{

    class SingleChainAnnotatorCpp : public AnnotatorBaseClassCpp {
        public:
            SingleChainAnnotatorCpp(std::vector<std::string> chains = {"H", "K", "L"},
                    std::string scheme = "imgt", bool compress_init_gaps = false,
                    std::string consensus_filepath = "");


            /// @brief Numbers an input sequence
            /// @param query_sequence Sequence to number
            /// @return A tuple containing the numbering, percent identity,
            ///         chain type and error message.
            std::tuple<std::vector<std::string>, double, std::string,
                std::string> analyze_seq(std::string);

            /// @brief Numbers a list of input sequences
            /// @param query_sequences A vector of sequences to number.
            /// @return A vector of tuples of the same length as the input vector
            ///         of sequences. Each tuple contains the numbering, percent identity,
            ///         chain type and error message for the corresponding sequence.
            std::vector<std::tuple<std::vector<std::string>, double, std::string,
                std::string>> analyze_seqs(std::vector<std::string> sequences);
        

        protected:
            std::vector<std::string> chains;
            std::string scheme;
            bool compress_init_gaps;

            std::vector<NumberingTools::IGAligner> scoring_tools;

            /// @brief Aligns a subregion of an input sequence (for situations where
            ///        a subregion needs to be extracted).
            int align_input_subregion(std::tuple<std::vector<std::string>, double,
                    std::string, std::string> &best_result, double &best_identity,
                    std::string &query_sequence);

    };

}

#endif
