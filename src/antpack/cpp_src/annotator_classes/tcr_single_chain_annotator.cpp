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
// C++ headers
#include <tuple>
#include <vector>
#include <string>
#include <memory>

// Library headers

// Project headers
#include "tcr_single_chain_annotator.h"
#include "../utilities/utilities.h"
#include "numbering_constants.inc"



namespace TCRNumberingTools {

TCRSingleChainAnnotatorCpp::TCRSingleChainAnnotatorCpp(
        std::vector<std::string> chains,
        std::string scheme,
        std::string consensus_filepath):
chains(chains),
scheme(scheme) {
}

std::tuple<std::vector<std::string>, double, std::string,
            std::string> TCRSingleChainAnnotatorCpp::analyze_seq(
                    std::string sequence) {
}


/// @brief Aligns the input sequence, which may be either the full sequence or a subregion of it.
///
int TCRSingleChainAnnotatorCpp::align_input_subregion(
        std::tuple<std::vector<std::string>, double,
        std::string, std::string> &best_result,
        std::string &query_sequence,
        std::string preferred_chain) {
    auto queryAsIdx = std::make_unique<int[]>(query_sequence.length());
    bool fwgxg_error = false;

    // Set initial percent identity to minimum.
    std::get<1>(best_result) = 0.0;

    if (!SequenceUtilities::convert_x_sequence_to_array(queryAsIdx.get(),
                query_sequence))
        return NumberingTools::INVALID_SEQUENCE;

    // If preferred chain is specified (i.e. is not ""), we can
    // only align to one scoring tool...neat! In that case, score
    // using that tool then return, AS LONG AS the preferred chain
    // matches one of the chains specified for this object.
    if (preferred_chain != "") {
        auto first_chain = this->chains.begin();
        auto last_chain = this->chains.end();
        auto selection = std::find(first_chain, last_chain,
                preferred_chain);

        if (selection != this->chains.end()) {
            size_t aligner_id = selection - this->chains.begin();
            std::vector<std::string> final_numbering;
            std::string error_message = "";
            double percent_identity = -1;

            this->scoring_tools[aligner_id].align(
                    query_sequence, queryAsIdx.get(), final_numbering,
                    percent_identity, error_message);

            std::get<0>(best_result) = final_numbering;
            std::get<1>(best_result) = percent_identity;
            std::get<2>(best_result) = this->scoring_tools[aligner_id].
                get_chain_name();
            std::get<3>(best_result) = error_message;
            return NumberingTools::VALID_SEQUENCE;
        }
    }



    // Otherwise, if no preferred chain was specified or if the preferred
    // chain is not one this object can handle (e.g. it was initialized
    // with chain H only but preferred is K), we need to
    // try all available aligners. Save the one that yields the highest
    // percent identity.

    for (size_t i=0; i < this->scoring_tools.size(); i++) {
        std::vector<std::string> final_numbering;
        std::string error_message = "";
        double percent_identity = -1;

        final_numbering.reserve(query_sequence.length() * 2);

        this->scoring_tools[i].align(query_sequence,
                    queryAsIdx.get(), final_numbering, percent_identity,
                    error_message);
        if (percent_identity > std::get<1>(best_result)) {
            std::get<0>(best_result) = final_numbering;
            std::get<1>(best_result) = percent_identity;
            std::get<2>(best_result) = this->scoring_tools[i].get_chain_name();
            std::get<3>(best_result) = error_message;
        }
        if (std::get<3>(best_result).length() > 2) {
            if (std::get<3>(best_result).substr(0, 2) == "> ")
                fwgxg_error = true;
        }
    }
    if (fwgxg_error && std::get<1>(best_result) < 0.75)
        return NumberingTools::POSSIBLE_FWGXG_ERROR_ON_ALIGNMENT;
    return NumberingTools::VALID_SEQUENCE;
}


/// @brief TCRSingleChainAnnotatorCpp function which numbers a list
///        of input sequences. Essentially a wrapper on analyze_seq
///        that is convenient if user wants to pass a list.
std::vector<std::tuple<std::vector<std::string>, double, std::string,
        std::string>> TCRSingleChainAnnotatorCpp::analyze_seqs(std::vector<std::string> sequences) {
    std::vector<std::tuple<std::vector<std::string>, double, std::string,
            std::string>> output_results;

    for (size_t i=0; i < sequences.size(); i++)
        output_results.push_back(this->analyze_seq(sequences[i]));

    return output_results;
}


}  // namespace TCRNumberingTools
