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
#include "vj_aligner.h"
#include "numbering_constants.inc"



namespace nb = nanobind;


namespace SequenceAnnotators {

class SingleChainAnnotatorCpp : public AnnotatorBaseClassCpp {
 public:
        SingleChainAnnotatorCpp(std::vector<std::string> chains,
                std::string scheme, std::string consensus_filepath,
                std::unordered_map<std::string, size_t> nterm_kmers);

        /// @brief Numbers a single input sequence by internally
        /// routing to mab_analyze_seq or tcr_analyze_seq depending
        /// on what the class was initialized to do.
        /// @param sequence The input sequence.
        /// @return A tuple of numbering, percent_identity, chain_type
        /// and error_message.
        std::tuple<std::vector<std::string>, double, std::string,
                std::string> analyze_seq(std::string sequence);

        /// @brief Numbers a single input sequence when the object has been
        /// initialized to annotate antibody sequences specifically.
        /// @param sequence The input sequence.
        /// @return A tuple of numbering, percent_identity, chain_type
        /// and error_message.
        std::tuple<std::vector<std::string>, double, std::string,
                std::string> mab_analyze_seq(std::string sequence);

        /// @brief Numbers a single input sequence when the object has been
        /// initialized to annotate TCR sequences specifically.
        /// @param sequence The input sequence.
        /// @return A tuple of numbering, percent_identity, chain_type
        /// and error_message.
        std::tuple<std::vector<std::string>, double, std::string,
                std::string> tcr_analyze_seq(std::string sequence);

        /// @brief Aligns a subregion of an input sequence (for situations where
        /// a subregion needs to be extracted) in cases where the object has
        /// been initialized to analyze mab sequences.
        /// @param best_result The tuple of numbering, percent_identity,
        /// chain_type and error_message in which output will be stored.
        /// @param query_sequence The input sequence.
        /// @param preferred_chain Indicates if a specific chain has already
        /// been found (based on kmer counting) to be a better match; if
        /// so, only this one needs to be aligned.
        int mab_align_input_subregion(std::tuple<std::vector<std::string>,
                    double, std::string, std::string> &best_result,
                std::string &query_sequence,
                std::string preferred_chain);


        /// @brief Numbers a list of input sequences
        /// @param query_sequences A vector of sequences to number.
        /// @return A vector of tuples of the same length as the input vector
        /// of sequences. Each tuple contains the numbering, percent identity,
        /// chain type and error message for the corresponding sequence.
        std::vector<std::tuple<std::vector<std::string>, double, std::string,
            std::string>> analyze_seqs(std::vector<std::string> sequences);

 protected:
        std::vector<std::string> chains;
        std::string scheme;

        // The following two variables are used ONLY if this annotator
        // handles mAbs. If it handles tcrs, see below.
        std::unique_ptr<PrefilteringRoutines::PrefilteringTool> boundary_finder;
        std::vector<NumberingTools::IGAligner> scoring_tools;

        // The following variables are used ONLY if this annotator
        // handles tcrs. If it handles mAbs, see above.
        std::vector<NumberingTools::VJAligner> tcr_scoring_tools;

        /// @brief Aligns a subregion of an input sequence (for situations where
        /// a subregion needs to be extracted) in cases where the object has
        /// been initialized to analyze TCR sequences.
        /// @param best_result The tuple of numbering, percent_identity,
        /// chain_type and error_message in which output will be stored.
        /// @param query_sequence The input sequence.
        /// @param queryAsIdx A pointer to an array containing the query
        /// sequence encoded as integers.
        /// @param preferred_chain Indicates which chain should be used. Note
        /// that in contrast to mab_align_input_subregion, for TCRs the
        /// preferred chain must be specified.
        /// @param preferred_vgene Indicates which vgene to use (notice that
        /// this is different from mab_align_input_subregion).
        /// @param preferred_jgene Indicates which jgene to use (notice that
        /// this is different from mab_align_input_subregion).
        int tcr_align_input_subregion(
            std::tuple<std::vector<std::string>, double,
                std::string, std::string> &best_result,
            std::string &query_sequence,
            int *queryAsIdx,
            const std::string &preferred_chain,
            const int &preferred_vgene,
            const int &preferred_jgene);
};

}  // namespace SequenceAnnotators

#endif
