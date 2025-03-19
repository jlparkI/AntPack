/* The IGAligner runs a highly efficient profile alignment for mAb sequences.
 * Copyright (C) 2025 Jonathan Parkinson
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
#ifndef SRC_ANTPACK_CPP_SRC_ANNOTATOR_CLASSES_IG_ALIGNER_H_
#define SRC_ANTPACK_CPP_SRC_ANNOTATOR_CLASSES_IG_ALIGNER_H_

// C++ headers
#include <vector>
#include <string>
#include <algorithm>
#include <filesystem>
#include <tuple>
#include <set>
#include <array>
#include <memory>

// Library headers
#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>

// Project headers
#include "../utilities/consensus_file_utilities.h"
#include "../utilities/utilities.h"
#include "numbering_constants.inc"


namespace nb = nanobind;

namespace NumberingTools {

class IGAligner {
 public:
    IGAligner(std::string consensus_filepath,
        std::string chain_name, std::string scheme);


    /// @brief Numbers an input sequence
    /// @param query_sequence Sequence to number
    /// @param encoded_sequence Pointer to array of same length
    ///        as query sequence with the encoded query seq.
    /// @param final_numbering vector that will store the generated
    ///        numbering
    /// @param percent_identity the percent identity to the template.
    /// @param error_message the error message (if none, "").
    void align(std::string query_sequence, int *encoded_sequence,
            std::vector<std::string> &final_numbering,
            double &percent_identity, std::string &error_message);

    /// @brief Convenience function to obtain the chain with which this
    ///        analyzer is associated.
    std::string get_chain_name();


    // Allowed error codes. These will be mapped to strings
    // which explain in more detail.
    enum allowedErrorCodes {noError = 0, invalidSequence = 1,
        fatalRuntimeError = 2,
        tooManyInsertions  = 3, alignmentWrongLength = 4,
        unacceptableConservedPositions = 5};


 protected:
    // Default gap penalties for gaps at the beginning and end of the sequence.
    // template is a weak penalty for placing gaps outside the numbering,
    // while query is a weak penalty for placing gaps in the numbering
    // (i.e. skipping numbers).

    const std::string chain_name;
    std::string scheme;

    std::unique_ptr<double[]> score_array;
    size_t score_arr_shape[2];
    int num_positions;
    int num_restricted_positions;

    std::vector<std::set<char>> consensus_map;
    std::vector<int> highly_conserved_positions;
    std::array<std::string, 6> error_code_to_message {{"",
        "Sequence contains invalid characters",
        "Fatal runtime error in IGAligner. Unusual. Please report",
        "> max allowed insertions -- suggests alignment error or unusual sequence",
        "Alignment length != length of input sequence. Unusual. Please report.",
        "Unexpected AA at conserved position."}};


    // Alphabet for numbering insertions. Our preference would
    // be to number insertions as _1, _2 etc, but most numbering
    // programs use letters, so we do the same here for consistency
    // and ease of comparison. This does limit the number of possible
    // insertions BUT if there is ever a case where more than 70
    // insertions have occurred at some point in the sequence, there
    // is almost certainly an error in the alignment, and
    // that's how we handle this eventuality if it occurs. (Indeed,
    // it is very rare to have more than a few at any given position,
    // so the size of alphabet here is overkill.)
    const std::vector<std::string> alphabet = {
                    "A", "B", "C", "D", "E", "F", "G", "H", "I",
                    "J", "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T", "U",
                    "V", "W", "X", "Y", "Z", "AA", "BB", "CC", "DD", "EE",
                    "FF", "GG", "HH", "II", "JJ", "KK", "LL", "MM", "NN",
                    "OO", "PP", "QQ", "RR", "SS", "TT", "UU", "VV",
                    "WW", "XX", "YY", "ZZ", "AAA", "BBB", "CCC", "DDD",
                    "EEE", "FFF", "GGG", "HHH", "III", "JJJ", "KKK", "LLL",
                    "MMM", "NNN", "OOO", "PPP", "QQQ", "RRR", "SSS", "TTT",
                    "UUU", "VVV", "WWW", "XXX", "YYY", "ZZZ"};


    /// @brief Fills in the scoring table to construct an alignment.
    /// @param path_trace The array that will store the best path found (aka
    /// the scoring table).
    /// @param query_seq_len The length of the query sequence.
    /// @param encoded_sequence A pointer to the array containing the encoded
    /// sequence.
    /// @param num_elements The number of elements in the scoring table.
    void fill_needle_scoring_table(uint8_t *path_trace,
            int query_seq_len, int row_size,
            const int *encoded_sequence,
            int &numElements);
};

}  // namespace NumberingTools

#endif
