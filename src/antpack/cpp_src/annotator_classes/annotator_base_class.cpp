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
#include <unordered_map>

// Project headers
#include "annotator_base_class.h"
#include "../utilities/cdr_assignment_utilities.h"

namespace SequenceAnnotators {


AnnotatorBaseClassCpp::AnnotatorBaseClassCpp(std::string scheme):
scheme(scheme) {
}





/// AnnotatorBaseClassCpp function which sorts a list of
/// position codes (essentially by calling the corresponding
/// function in utilities). Essentially a wrapper on the
/// corresponding utility function that can be accessed by
/// Python or outside callers.
std::vector<std::string> AnnotatorBaseClassCpp::sort_position_codes(
        std::vector<std::string> position_code_list) {
    std::vector<std::string> cleanedCodes, sortedCodes;
    int isValid;

    for (auto &code : position_code_list) {
        if (code != "-")
            cleanedCodes.push_back(code);
    }
    if (!SequenceUtilities::sort_position_codes_utility(cleanedCodes,
                this->scheme, sortedCodes)) {
        throw std::runtime_error(std::string("One or more of the "
                    "supplied position codes is invalid given the "
                    "specified scheme."));
    }
    return sortedCodes;
}



/// AnnotatorBaseClassCpp function which builds a multiple sequence alignment from
/// a list of previously numbered sequences. The sequences must all be of the
/// same chain type. Essentially a wrapper on the corresponding utility function
/// that can be accessed by Python or outside callers.
std::tuple<std::vector<std::string>, std::vector<std::string>>
AnnotatorBaseClassCpp::build_msa(std::vector<std::string> sequences,
        std::vector<std::tuple<std::vector<std::string>, double,
        std::string, std::string>> annotations,
        bool add_unobserved_positions) {
    std::vector<std::string> positionCodes, alignedSeqs;

    if (!SequenceUtilities::build_msa_utility(sequences, annotations, positionCodes,
            alignedSeqs, this->scheme, add_unobserved_positions)) {
        throw std::runtime_error(std::string("A fatal error occured when "
                    "building an MSA -- please report."));
    }

    return std::tuple<std::vector<std::string>,
           std::vector<std::string>>{positionCodes, alignedSeqs};
}


/// AnnotatorBaseClassCpp function which trims an alignment to
/// remove gap positions at either end. Essentially a wrapper
/// on the corresponding utility function that can
/// be accessed by Python or outside callers.
std::tuple<std::string, std::vector<std::string>, int, int>
    AnnotatorBaseClassCpp::trim_alignment(std::string sequence,
        std::tuple<std::vector<std::string>, double,
        std::string, std::string> alignment) {
    int exstart, exend;
    std::vector<char> trimmedSeq;
    std::vector<std::string> trimmedAlignment;

    if (!SequenceUtilities::trim_alignment_utility(sequence,
                alignment, trimmedAlignment, exstart, exend,
                trimmedSeq)) {
        throw std::runtime_error(std::string("Invalid sequence / "
                    "alignment pairing supplied."));
    }

    std::string trimmedSeqStr(trimmedSeq.begin(), trimmedSeq.end());

    return std::tuple<std::string, std::vector<std::string>,
           int, int>{trimmedSeqStr, trimmedAlignment, exstart, exend};
}


/// Wraps the assign_cdr_labels function in SequenceUtilities
/// for access by Python callers.
std::vector<std::string> AnnotatorBaseClassCpp::assign_cdr_labels(
        std::vector<std::string> numbering, std::string chain,
        std::string scheme) {
    std::vector<std::string> cdr_labeling;
    if (scheme == "") {
        CDRConversionUtilities::assign_cdr_labels(numbering, chain,
            cdr_labeling, this->scheme, this->scheme);
    } else {
        CDRConversionUtilities::assign_cdr_labels(numbering, chain,
            cdr_labeling, this->scheme, scheme);
    }
    return cdr_labeling;
}

}  // namespace SequenceAnnotators
