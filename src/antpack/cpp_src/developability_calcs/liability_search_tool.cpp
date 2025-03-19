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
#include <vector>
#include <stdexcept>
#include <string>
#include <utility>

// Library headers

// Project headers
#include "liability_search_tool.h"
#include "../utilities/cdr_assignment_utilities.h"
#include "../annotator_classes/annotator_base_class.h"


namespace LiabilitySearch {


LiabilitySearchToolCpp::LiabilitySearchToolCpp() {
}


/// @brief Finds potential liabilities in the sequence.
/// @param alignment A tuple containing the output of AntPack's single
/// or paired chain annotator.
/// @param sequence The sequence as a string. Should not contain gaps.
/// @param scheme One of 'imgt', 'aho', 'martin', 'kabat'. This must be
/// the scheme that was used to number the input sequence.
/// @param cdr_scheme One of 'imgt', 'aho', 'martin', 'kabat', 'north.
/// This is the cdr definitions that are used. Note that you can use
/// a different set of cdr definitions than numbering scheme (e.g. number
/// with IMGT and define CDRs using Kabat) although usually you will
/// want this to be the same as 'scheme'.
/// @return A vector of tuples. Each tuple contains a two-tuple of ints
/// (the start and end of the liability) and a string describing
/// the type of liability.
std::vector<std::pair<std::pair<int, int>, std::string>>
LiabilitySearchToolCpp::analyze_seq(std::string sequence,
                std::tuple<std::vector<std::string>,
                double, std::string, std::string> alignment,
                std::string scheme,
                std::string cdr_scheme) {
    // As always, exceptions thrown here get sent back to Python
    // if this is used within the Nanobind wrapper.

    if (std::get<2>(alignment) != "H" && std::get<2>(alignment) != "K" &&
            std::get<2>(alignment) != "L") {
        throw std::runtime_error(std::string("Only antibody chains are "
                    "supported."));
    }

    if (sequence.length() != std::get<0>(alignment).size()) {
        throw std::runtime_error(std::string("Alignment and sequence "
                    "must have the same size."));
    }
    // Note that the cdr_assignment_utilities functions check that the input
    // alignment and scheme are valid. If they are not they will throw an
    // exception, which if this is used within the Nanobind wrapper will get
    // passed back to Python. For these purposes we assume
    std::vector<std::string> cdr_labeling;
    CDRConversionUtilities::assign_cdr_labels(std::get<0>(alignment),
            std::get<2>(alignment), cdr_labeling,
            scheme, cdr_scheme);

    if (cdr_labeling.size() != sequence.length()) {
        throw std::runtime_error(std::string("Alignment and cdr labeling are "
                    "not the same length! This is a serious internal error; "
                    "please report."));
    }

    size_t sequence_length = sequence.length();
    std::vector<std::pair<std::pair<int, int>, std::string>> output;

    for (size_t i=0; i < sequence_length; i++) {
        if (cdr_labeling.at(i) == "-" || cdr_labeling.at(i).length() < 4)
            continue;

        if (cdr_labeling.at(i).substr(0, 3) == "cdr") {
            // (re.compile(r"[MW]"), True,
            // "Methionine / Tryptophan oxidation moderate risk"),
            if (sequence.at(i) == 'M' || sequence.at(i) == 'W') {
                std::pair<std::pair<int, int>, std::string> error =
                    { {i, i+1}, "Methionine / Tryptophan oxidation "
                        "moderate risk"};
                output.push_back(error);
                // In this case, if the letter is MW, no other liabilities
                // are applicable, so continue.
                continue;
            }
            if (i+2 <= sequence_length) {
                // (re.compile(r"TS"), True, "pH-dependent hydrolysis risk"),
                std::string next2 = sequence.substr(i, 2);
                if (next2 == "TS") {
                    std::pair<std::pair<int, int>, std::string> error =
                        { {i, i+2}, "pH-dependent hydrolysis risk"};
                    output.push_back(error);
                    // In this case, if the letters are TS, no other liabilities
                    // are applicable, so continue.
                    continue;
                }
                if (next2.at(1) == 'N') {
                    if (next2.at(0) == 'S' || next2.at(0) == 'T' ||
                            next2.at(0) == 'K') {
                        // (re.compile(r"[STK]N"), True,
                        // "Deamidation (low risk)"),
                        std::pair<std::pair<int, int>, std::string> error =
                            { {i, i+2}, "Deamidation (low risk)"};
                        output.push_back(error);
                    }
                }
                if (next2.at(0) == 'N' || next2.at(0) == 'D') {
                    if (next2.at(1) == 'P') {
                        // (re.compile(r"[ND]P"), True,
                        // "pH-dependent hydrolysis (moderate risk)"),
                        std::pair<std::pair<int, int>, std::string> error =
                            { {i, i+2}, "pH-dependent hydrolysis "
                                "(moderate risk)"};
                        output.push_back(error);
                    }
                }
                if (next2.at(0) == 'D') {
                    if (next2.at(1) == 'D' || next2.at(1) == 'G' ||
                            next2.at(1) == 'H' || next2.at(1) == 'S' ||
                            next2.at(1) == 'T') {
                        // (re.compile(r"D[DGHST]"), True,
                        // "Isomerization (elevated risk)"),
                        std::pair<std::pair<int, int>, std::string> error =
                            { {i, i+2}, "Isomerization (elevated risk)"};
                        output.push_back(error);
                    }
                }
                if (next2.at(0) == 'N') {
                    // (re.compile(r"N[GS]"), True,
                    // "Deamidation (elevated risk)"),
                    // (re.compile(r"N[AHNT]"), True,
                    // "Deamidation (moderate risk)"),
                    if (next2.at(1) == 'G' || next2.at(1) == 'S') {
                        std::pair<std::pair<int, int>, std::string> error =
                            { {i, i+2}, "Deamidation (elevated risk)"};
                        output.push_back(error);
                    }
                    if (next2.at(1) == 'A' || next2.at(1) == 'H' ||
                            next2.at(1) == 'N' || next2.at(1) == 'T') {
                        std::pair<std::pair<int, int>, std::string> error =
                            { {i, i+2}, "Deamidation (elevated risk)"};
                        output.push_back(error);
                    }
                }
            }
        }
        if (i+3 <= sequence_length) {
            std::string next3 = sequence.substr(i, 3);
            // (re.compile(r"N[^P][ST]"), False, "N-glycosylation"),
            if (next3.at(0) == 'N' && next3.at(1) != 'P') {
                if (next3.at(2) == 'S' || next3.at(2) == 'T') {
                    std::pair<std::pair<int, int>, std::string> error =
                        { {i, i+3}, "N-glycosylation risk"};
                    output.push_back(error);
                }
            }
        }
        if (sequence.at(i) == 'C') {
            if (std::get<0>(alignment).at(i) != "23" &&
                    std::get<0>(alignment).at(i) != "104") {
                std::pair<std::pair<int, int>, std::string> error =
                        { {i, i+3}, "Unusual cysteine"};
                output.push_back(error);
            }
        }
    }
    return output;
}


}  // namespace LiabilitySearch
