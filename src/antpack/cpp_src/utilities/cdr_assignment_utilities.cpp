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
#include <string>
#include <vector>
#include <array>
#include <stdexcept>

// Library headers

// Project headers
#include "cdr_assignment_utilities.h"
#include "../annotator_classes/numbering_constants.inc"


namespace SequenceUtilities {



/// @brief Assigns cdr labels to an input sequence.
/// @param numbering A list of numbering codes valid for the
///        scheme specified when the class was created. This is
///        the first element of the tuple returned by 'analyze_seq'.
/// @param chain A valid chain ('H', 'K', 'L'). 'K' and 'L' are treated
///        as equivalent.
/// @param cdr_labels A vector in which the output will be stored.
/// @param scheme The scheme -- one of 'imgt', 'aho', 'martin', 'kabat'.
void assign_cdr_labels(const std::vector<std::string> &numbering,
        const std::string &chain,
        std::vector<std::string> &cdr_labeling,
        const std::string &scheme) {
    std::array<int, 6> current_breakpoints;
    int numeric_portion;
    std::string current_label;
    size_t current_token = 0;
    int next_breakpoint;

    std::array<std::string, 7> cdr_region_labels = {"fmwk1", "cdr1",
            "fmwk2", "cdr2", "fmwk3", "cdr3", "fmwk4"};

    if (scheme == "imgt") {
        if (chain == "H" || chain == "L" || chain == "K")
            current_breakpoints = NumberingTools::IMGT_CDR_BREAKPOINTS;
        else
            throw std::runtime_error(std::string("Unrecognized chain or "
                        "scheme supplied."));
    } else if (scheme == "martin") {
        if (chain == "H")
            current_breakpoints = NumberingTools::MARTIN_HEAVY_CDR_BREAKPOINTS;
        else if (chain == "L" || chain == "K")
            current_breakpoints = NumberingTools::MARTIN_LIGHT_CDR_BREAKPOINTS;
        else
            throw std::runtime_error(std::string("Unrecognized chain or "
                        "scheme supplied."));
    } else if (scheme == "kabat") {
        if (chain == "H")
            current_breakpoints = NumberingTools::KABAT_HEAVY_CDR_BREAKPOINTS;
        else if (chain == "L" || chain == "K")
            current_breakpoints = NumberingTools::KABAT_LIGHT_CDR_BREAKPOINTS;
        else
            throw std::runtime_error(std::string("Unrecognized chain or "
                        "scheme supplied."));
    } else if (scheme == "aho") {
        if (chain == "H" || chain == "L" || chain == "K")
            current_breakpoints = NumberingTools::AHO_CDR_BREAKPOINTS;
        else
            throw std::runtime_error(std::string("Unrecognized chain or "
                        "scheme supplied."));
    }

    next_breakpoint = current_breakpoints.at(current_token);
    current_label = cdr_region_labels[current_token];

    for (size_t i=0; i < numbering.size(); i++) {
        if (numbering.at(i) == "-") {
            cdr_labeling.push_back("-");
            continue;
        }
        // This will throw if the string starts with a letter but will
        // otherwise extract the integer piece.
        // AntPack never places a letter at the start of the code,
        // so this will not happen unless
        // the user has passed some altered / corrupted input.
        try {
            numeric_portion = std::stoi(numbering[i]);
        }
        catch (...) {
            throw std::runtime_error(std::string("An invalid position "
                        "code was passed. The alignment passed to this "
                        "function should be unaltered output from "
                        "analyze_seq and not some other procedure."));
        }
        if (numeric_portion >= next_breakpoint) {
            if ( current_token < (current_breakpoints.size() - 1) ) {
                while (current_token < (current_breakpoints.size() - 1) &&
                        numeric_portion >= next_breakpoint) {
                    current_token += 1;
                    next_breakpoint = current_breakpoints.at(current_token);
                    current_label =
                        cdr_region_labels[current_token];
                }
                if (numeric_portion >= next_breakpoint) {
                    // Set next_breakpoint to an arbitrarily
                    // high, unachievable number.
                    next_breakpoint = 10000;
                    current_token += 1;
                    current_label = cdr_region_labels.at(
                        cdr_region_labels.size() - 1);
                }
            } else {
                // Set next_breakpoint to an arbitrarily high,
                // unachievable number.
                next_breakpoint = 10000;
                current_token += 1;
                current_label = cdr_region_labels.at(
                    cdr_region_labels.size() - 1);
            }
        }
        cdr_labeling.push_back(current_label);
    }
}

}  // namespace SequenceUtilities
