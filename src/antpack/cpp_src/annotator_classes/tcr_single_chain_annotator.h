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
#ifndef TCR_SINGLE_CHAIN_ANNOTATOR_HEADER_H
#define TCR_SINGLE_CHAIN_ANNOTATOR_HEADER_H
// C++ headers
#include <tuple>
#include <vector>
#include <string>

// Library headers

// Project headers


namespace TCRNumberingTools {

class TCRSingleChainAnnotatorCpp {
 public:
    TCRSingleChainAnnotatorCpp(std::vector<std::string> chains,
            std::string scheme, std::string consensus_filepath);

    /// @brief Numbers a subregion of an input sequence.
    /// @param best_result A reference to a tuple in which the results
    ///        will be stored.
    /// @param query_sequence The sequence to number.
    /// @param preferred_chain A string indicating which chain should be
    ///        used.
    /// @return An int indicating whether an error occurred.
    int align_input_subregion(
        std::tuple<std::vector<std::string>, double,
        std::string, std::string> &best_result,
        std::string &query_sequence,
        std::string preferred_chain);

    /// @brief Numbers an input sequence
    /// @param query_sequence Sequence to number
    /// @return A tuple containing the numbering, percent identity,
    ///         chain type and error message.
    std::tuple<std::vector<std::string>, double, std::string,
        std::string> analyze_seq(std::string sequence);

    /// @brief Numbers a list of input sequences
    /// @param query_sequences A vector of sequences to number.
    /// @return A vector of tuples of the same length as the input vector
    ///         of sequences. Each tuple contains the numbering, percent
    ///         identity, chain type and error message for the corresponding
    ///         sequence.
    std::vector<std::tuple<std::vector<std::string>, double, std::string,
        std::string>> analyze_seqs(std::vector<std::string> sequences);


 protected:
    std::vector<std::string> chains;
    std::string scheme;
};


}  // namespace TCRNumberingTools


#endif  // SRC_ANTPACK_CPP_SRC_ANNOTATOR_CLASSES_TCR_SINGLE_CHAIN_ANNOTATOR_H_
