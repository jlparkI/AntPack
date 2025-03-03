/* Copyright (C) 2025 Jonathan Parkinson
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef SINGLE_CHAIN_ANNOTATOR_HEADER_H
#define SINGLE_CHAIN_ANNOTATOR_HEADER_H

// C++ headers
#include <vector>
#include <string>
#include <tuple>
#include <memory>
#include <unordered_map>

// Library headers
#include <nanobind/ndarray.h>

// Project headers
#include "prefiltering_tool.h"
#include "annotator_base_class.h"
#include "../utilities/utilities.h"
#include "ig_aligner.h"
#include "numbering_constants.inc"



namespace nb = nanobind;


namespace NumberingTools {

class SingleChainAnnotatorCpp : public AnnotatorBaseClassCpp {
 public:
        SingleChainAnnotatorCpp(std::vector<std::string> chains,
                std::string scheme, std::string consensus_filepath,
                std::unordered_map<std::string, size_t> nterm_kmers);


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

        /// @brief Aligns a subregion of an input sequence (for situations where
        ///        a subregion needs to be extracted).
        int align_input_subregion(std::tuple<std::vector<std::string>, double,
                std::string, std::string> &best_result, std::string &query_sequence,
                std::string preferred_chain);

 protected:
        std::vector<std::string> chains;
        std::string scheme;
        std::unique_ptr<PrefilteringRoutines::PrefilteringTool> boundary_finder;

        std::vector<NumberingTools::IGAligner> scoring_tools;
};

}  // namespace NumberingTools

#endif
