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
#ifndef SRC_ANTPACK_CPP_SRC_ANNOTATOR_CLASSES_ANNOTATOR_BASE_CLASS_H_
#define SRC_ANTPACK_CPP_SRC_ANNOTATOR_CLASSES_ANNOTATOR_BASE_CLASS_H_

// C++ headers
#include <vector>
#include <array>
#include <string>
#include <utility>
#include <tuple>
#include <unordered_map>

// Project headers
#include "numbering_constants.inc"
#include "../utilities/utilities.h"




namespace SequenceAnnotators {


class AnnotatorBaseClassCpp {
 public:
    AnnotatorBaseClassCpp(std::string scheme);

    /// @brief Sorts a list of position codes, respecting the
    ///        rules of the numbering scheme.
    /// @param position_code_list A list of position codes to
    ///        sort. Each must be a string.
    /// @return A vector of position codes that have been sorted.
    std::vector<std::string> sort_position_codes(std::vector<std::string>
            position_code_list);

    /// @brief Takes a list of numbered sequences and converts this
    ///        to a multiple sequence alignment which can be easily
    ///        viewed or written to file.
    /// @param sequences A vector of strings each of which is one of
    ///        the numbered sequences.
    /// @param annotations A vector of tuples of the same length as
    ///        sequences. Each tuple has four elements: a vector of
    ///        numberings of the same length as that string, a percent
    ///        identity, a chain type and an error message. This tuple
    ///        is returned by the analyze_seq methods of child classes.
    /// @param add_unobserved_positions If True, positions expected
    ///        for a given scheme are added even if not observed for
    ///        any of the input alignments.
    /// @return A tuple of two vectors. The first contains a consensus
    ///         set of position codes. The second contains all of the
    ///         input sequences, gapped to be the same length.
    std::tuple<std::vector<std::string>, std::vector<std::string>>
        build_msa(std::vector<std::string> sequences,
        std::vector<std::tuple<std::vector<std::string>,
        double, std::string, std::string>> annotations,
        bool add_unobserved_positions = false);

    /// @brief Takes an input sequence and corresponding numbering
    ///        and trims it so that gaps at the c- and n-terminals
    ///        are removed.
    /// @param sequence A string which is the sequence to align.
    /// @param alignment A tuple containing four elements -- the
    ///        numbering, percent identity, chain type and error message.
    ///        This tuple is generated / returned by the analyze_seq methods
    ///        of child classes.
    /// @return Returns a tuple with four elements: the trimmed sequence,
    ///         the trimmed alignment,
    ///         the first position in the input sequence that was not
    ///         trimmed and the last position in the input sequence that
    ///         was not trimmed.
    std::tuple<std::string, std::vector<std::string>, int, int>
        trim_alignment(std::string sequence,
        std::tuple<std::vector<std::string>, double, std::string,
        std::string> alignment);

    /// @brief Assigns cdr labels to an input sequence.
    /// @param numbering A list of numbering codes valid for the
    ///        scheme specified when the class was created. This is
    ///        the first element of the tuple returned by 'analyze_seq'.
    /// @param chain A valid chain ('H', 'K', 'L'). 'K' and 'L' are treated
    ///        as equivalent.
    /// @param scheme Either "" or a valid scheme. If "", the scheme
    ///        used when the class was created is used. Otherwise
    ///        cross-scheme assignment occurs (the sequence was numbered
    ///        using one scheme but the CDRs are assigned using another).
    /// @return Returns a vector of strings '-', 'fmwk1', 'cdr1',
    ///         'fmwk2', 'cdr2' etc. to indicate the region to which
    ///         each position belongs.
    std::vector<std::string> assign_cdr_labels(
            std::vector<std::string> numbering,
                std::string chain,
                std::string scheme = "");


 protected:
        std::string scheme;
};

}  // namespace SequenceAnnotators

#endif  // SRC_ANTPACK_CPP_SRC_ANNOTATOR_CLASSES_ANNOTATOR_BASE_CLASS_H_
