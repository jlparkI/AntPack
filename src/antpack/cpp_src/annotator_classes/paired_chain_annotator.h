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
#ifndef PAIRED_CHAIN_ANNOTATOR_HEADER_H
#define PAIRED_CHAIN_ANNOTATOR_HEADER_H

// C++ headers
#include <vector>
#include <string>
#include <algorithm>
#include <tuple>
#include <set>
#include <array>
#include <memory>
#include <unordered_map>

// Library headers


// Project headers
#include "prefiltering_tool.h"
#include "annotator_base_class.h"
#include "single_chain_annotator.h"
#include "numbering_constants.inc"



namespace SequenceAnnotators {


class PairedChainAnnotatorCpp : public AnnotatorBaseClassCpp {
 public:
        PairedChainAnnotatorCpp(std::string scheme,
                std::string consensus_filepath,
                std::unordered_map<std::string, size_t> nterm_kmers);

        std::pair<std::tuple<std::vector<std::string>, double, std::string, std::string>,
            std::tuple<std::vector<std::string>, double, std::string, std::string>>
            analyze_seq(std::string sequence);
        std::tuple<std::vector<std::tuple<std::vector<std::string>, double, std::string, std::string>>,
            std::vector<std::tuple<std::vector<std::string>, double, std::string, std::string>>>
            analyze_seqs(std::vector<std::string> sequences);


 protected:
        std::string scheme;
        std::unique_ptr<SequenceAnnotators::SingleChainAnnotatorCpp>
            light_chain_analyzer;
        std::unique_ptr<SequenceAnnotators::SingleChainAnnotatorCpp>
            heavy_chain_analyzer;
        std::unique_ptr<SequenceAnnotators::SingleChainAnnotatorCpp>
            analyzer;

        std::unique_ptr<PrefilteringRoutines::PrefilteringTool> boundary_finder;

        /// @brief Pads the alignment represented by the input
        ///        so it is the same length as the sequence.
        std::vector<std::string> pad_alignment(
                const std::string &query_sequence,
                const std::vector<std::string> &alignment,
                const int &align_start, const int &align_end);
};

}  // namespace SequenceAnnotators

#endif
