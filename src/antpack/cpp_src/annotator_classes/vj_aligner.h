/* The VJAligner class handles alignment of an input sequence by
 * finding the most appropriate V and J genes from a list then
 * using custom preconstructed scoring matrices. This is different
 * from IGAligner, which aligns to a profile (not individual VJ
 * genes).
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


#ifndef SRC_ANTPACK_CPP_SRC_ANNOTATOR_CLASSES_VJ_ALIGNER_H_
#define SRC_ANTPACK_CPP_SRC_ANNOTATOR_CLASSES_VJ_ALIGNER_H_

// C++ headers
#include <vector>
#include <string>
#include <algorithm>
#include <unordered_map>
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



class VJAligner {
 public:
    VJAligner(std::string consensus_filepath,
        std::string chain_name, std::string scheme,
        std::string receptor_type);


    /// @brief Identifies the vgene which has the most kmers in
    /// common with an input sequence.
    /// @param query_sequence The input sequence.
    /// @param encoded_sequence Pointer to an array containing the
    /// input sequence encoded as integers; must be same size as
    /// query_sequence.length().
    /// @param identity The score for matching to the best vgene;
    /// the result is stored in this reference.
    /// @param best_vgene_number The id of the best vgene that is
    /// found; the result is stored in this reference.
    /// @return Returns 1 (VALID_SEQUENCE) or 0 for an error.
    int identify_best_vgene(std::string &query_sequence,
            const int *encoded_sequence, int &identity, int &best_vgene_number);


    /// @brief Identifies the jgene which is the best match for
    /// an input sequence based on a sliding window, and finds
    /// the location in the query sequence where the jgene should
    /// likely be aligned.
    /// @param query_sequence The input sequence.
    /// @param encoded_sequence Pointer to an array containing the
    /// input sequence encoded as integers; must be same size as
    /// query_sequence.length().
    /// @param optimal_position The position at which the jgene
    /// alignment should most likely occur; the result is stored
    /// in this reference.
    /// @param identity The score for matching to the best jgene;
    /// the result is stored in this reference.
    /// @param best_jgene_number The id of the best jgene that is
    /// found; the result is stored in this reference.
    /// @return Returns 1 (VALID_SEQUENCE) or 0 for an error.
    int identify_best_jgene(std::string &query_sequence,
            const int *encoded_sequence,
            int &optimal_position, int &identity,
            int &best_jgene_number);


    /// @brief Aligns an input sequence to a specified template V and J gene.
    /// @param query_sequence The sequence to be aligned.
    /// @param encoded_sequence A pointer to an array of ints containing the
    /// sequence encoded as numbers.
    /// @param final_number The vector in which the final numbering (as a vector
    /// of strings) will be stored.
    /// @param percent_identity The variable in which the percent identity to
    /// the template vgene will be stored.
    /// @param error_message The variable in which the error message for the
    /// alignment (if any) will be stored.
    /// @param vgene_number The number of the vgene. Determines which vgene in
    /// this class's list of vgenes will be used. If the input value is
    /// unacceptable this will cause an exception to be thrown.
    /// @param jgene_number The number of the jgene. Determines which jgene in
    /// this class's list of jgenes will be used. If the input value is
    /// unacceptable this will cause an exception to be thrown.
    void align(std::string query_sequence, int *encoded_sequence,
            std::vector<std::string> &final_numbering,
            double &percent_identity, std::string &error_message,
            int vgene_number, int jgene_number);

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

    size_t vgene_score_arr_shape[3];
    size_t jgene_score_arr_shape[3];
    std::unique_ptr<double[]> vgene_score_array;
    std::unique_ptr<double[]> jgene_score_array;
    std::unique_ptr<int32_t[]> blosum_array;

    std::unique_ptr<int[]> encoded_v_window1;
    std::unique_ptr<int[]> encoded_v_window2;
    std::unique_ptr<int[]> encoded_v_window3;
    std::unique_ptr<int[]> encoded_j_window1;
    std::vector<std::string> vgenes;
    std::vector<std::string> jgenes;

    // Since we always use IMGT, num_positions is set to 128.
    int num_positions = 128;

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
    /// @param vgene_number The vgene we should use.
    /// @param jgene_number The jgene we should use.
    void fill_needle_scoring_table(uint8_t *path_trace,
            int query_seq_len, int row_size,
            const int *encoded_sequence,
            const int &numElements,
            const int &vgene_number,
            const int &jgene_number);
};

}  // namespace NumberingTools

#endif
