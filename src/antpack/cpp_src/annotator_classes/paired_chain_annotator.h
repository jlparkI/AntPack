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
                std::unordered_map<std::string, size_t> nterm_kmers,
                std::string receptor_type);

        /// @brief Numbers an input sequence that contains a
        /// paired heavy and light chain (or alternatively just
        /// a heavy or light chain). Returns results for both
        /// heavy and light chains (the result will be blank for
        /// one of the two if only one is present).
        /// @param sequence The sequence to be numbered.
        /// @return A pair of tuples, each containing numbering,
        /// percent identity, chain name and error message. The
        /// first tuple is for heavy and the second for light.
        std::pair<std::tuple<std::vector<std::string>, double, std::string, std::string>,
            std::tuple<std::vector<std::string>, double, std::string, std::string>>
            analyze_seq(std::string sequence);

        /// @brief Numbers a list of input sequences, each of which
        /// contains a paired heavy and light chain (or just
        /// a heavy or light chain). Returns results for both
        /// heavy and light chains (the result will be blank for
        /// one of the two if only one is present).
        /// @param sequences The list of sequences to be numbered.
        /// @return A tuple of vectors. Each vector is a vector
        /// of tuples. Each sub-tuple contains numbering,
        /// percent identity, chain name and error message for
        /// the corresponding sequence. The first vector is a
        /// vector of results for heavy chains and the second
        /// is for light.
        std::tuple<std::vector<std::tuple<std::vector<std::string>, double, std::string, std::string>>,
            std::vector<std::tuple<std::vector<std::string>, double, std::string, std::string>>>
            analyze_seqs(std::vector<std::string> sequences);


 protected:
        std::string scheme;
        // Receptor type must be one of mab, tcr.
        bool receptor_type_is_mab;
        std::unique_ptr<SequenceAnnotators::SingleChainAnnotatorCpp>
            light_chain_analyzer;
        std::unique_ptr<SequenceAnnotators::SingleChainAnnotatorCpp>
            heavy_chain_analyzer;
        std::unique_ptr<SequenceAnnotators::SingleChainAnnotatorCpp>
            analyzer;

        // The boundary finder is only used if receptor_type is mab.
        std::unique_ptr<PrefilteringRoutines::PrefilteringTool> boundary_finder;

        /// @brief Numbers an input sequence that contains a
        /// paired heavy and light chain (or alternatively just
        /// a heavy or light chain), IF the receptor type is mab.
        /// Returns results for both heavy and light chains (the
        /// result will be blank for one of the two if only one is
        /// present).
        /// @param sequence The sequence to be numbered.
        /// @return A pair of tuples, each containing numbering,
        /// percent identity, chain name and error message. The
        /// first tuple is for heavy and the second for light.
        std::pair<std::tuple<std::vector<std::string>, double, std::string, std::string>,
            std::tuple<std::vector<std::string>, double, std::string, std::string>>
            mab_analyze_seq(std::string &sequence);

        /// @brief Numbers an input sequence that contains a
        /// paired heavy and light chain (or alternatively just
        /// a heavy or light chain), IF the receptor type is tcr.
        /// Returns results for both heavy and light chains (the
        /// result will be blank for one of the two if only one is
        /// present).
        /// @param sequence The sequence to be numbered.
        /// @return A pair of tuples, each containing numbering,
        /// percent identity, chain name and error message. The
        /// first tuple is for heavy and the second for light.
        std::pair<std::tuple<std::vector<std::string>, double, std::string, std::string>,
            std::tuple<std::vector<std::string>, double, std::string, std::string>>
            tcr_analyze_seq(std::string &sequence);


        /// @brief Pads the alignment represented by the input
        /// so it is the same length as the sequence.
        /// @param query_sequence The input sequence
        /// @param alignment The numbering, which is not the same
        /// length as query_sequence and therefore needs to be
        /// padded.
        /// @param align_start The position in query_sequence at
        /// which alignment starts.
        /// @param align_end The position in query_sequence at which
        /// alignment ends.
        /// @return The updated padded numbering.
        std::vector<std::string> pad_alignment(
                const std::string &query_sequence,
                const std::vector<std::string> &alignment,
                const int &align_start, const int &align_end);
};

}  // namespace SequenceAnnotators

#endif
