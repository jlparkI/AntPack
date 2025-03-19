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


namespace CDRConversionUtilities {



/// @brief Assigns cdr labels to an input sequence.
/// @param numbering A list of numbering codes valid for the
///        scheme specified when the class was created. This is
///        the first element of the tuple returned by 'analyze_seq'.
/// @param chain A valid chain ('H', 'K', 'L'). 'K' and 'L' are treated
///        as equivalent.
/// @param cdr_labels A vector in which the output will be stored.
/// @param scheme The scheme -- one of 'imgt', 'aho', 'martin', 'kabat'.
/// @param cdr_scheme The scheme -- one of 'imgt', 'aho', 'martin', 'kabat',
///        'north' -- that should be used to construct CDRs. This may
///        or may not be the same as scheme.
void assign_cdr_labels(const std::vector<std::string> &numbering,
        const std::string &chain,
        std::vector<std::string> &cdr_labeling,
        const std::string &scheme,
        const std::string &cdr_scheme) {
    if (chain == "A" || chain == "D" || chain == "B" ||
            chain == "G") {
        if (scheme != "imgt" || cdr_scheme != "imgt") {
            throw std::runtime_error(std::string("For TCRs, only 'imgt' "
                        "is supported."));
        }
    } else if (chain == "H" || chain == "K" || chain == "L") {
        if (scheme != "martin" && scheme != "kabat" && scheme != "imgt" &&
                scheme != "aho") {
            throw std::runtime_error(std::string("Unsupported numbering "
                        "scheme."));
        }
        if (cdr_scheme != "martin" && cdr_scheme != "kabat" &&
                cdr_scheme != "imgt" && cdr_scheme != "aho"
                && cdr_scheme != "north") {
            throw std::runtime_error(std::string("Unsupported cdr "
                        "assignment scheme."));
        }
    } else {
        throw std::runtime_error(std::string("Unrecognized chain or "
                        "scheme supplied."));
    }

    std::array<int, 6> current_breakpoints;

    // TODO: Replace this with a constexpr map -- the following has very
    // poor readability. There are no other scheme / cdr scheme conversions
    // we will ever need to add to this, unless someone were to invent
    // a new antibody numbering scheme (hopefully not...)
    //
    // Notice incidentally that for IMGT and AHO, we do not need to know
    // chain even when cross-converting, and that these are the ONLY
    // two schemes that make sense for TCRs. Currently only IMGT is
    // supported for TCRs (see if statements above).

    if (scheme == "imgt") {
        if (cdr_scheme == "imgt") {
            current_breakpoints = IMGT_CDR_BREAKPOINTS;
        } else if (cdr_scheme == "kabat") {
            if (chain == "H")
                current_breakpoints = IMGT_NMBR_KABAT_H_CDR_BREAKPOINTS;
            else
                current_breakpoints = IMGT_NMBR_KABAT_L_CDR_BREAKPOINTS;
        } else if (cdr_scheme == "martin") {
            if (chain == "H")
                current_breakpoints = IMGT_NMBR_MARTIN_H_CDR_BREAKPOINTS;
            else
                current_breakpoints = IMGT_NMBR_MARTIN_L_CDR_BREAKPOINTS;
        } else if (cdr_scheme == "aho") {
            current_breakpoints = IMGT_NMBR_AHO_CDR_BREAKPOINTS;
        } else if (cdr_scheme == "north") {
            if (chain == "H")
                current_breakpoints = IMGT_NMBR_NORTH_H_CDR_BREAKPOINTS;
            else
                current_breakpoints = IMGT_NMBR_NORTH_L_CDR_BREAKPOINTS;
        } else {
            throw std::runtime_error(std::string("Unrecognized chain or scheme."));
        }
    }
    else if (scheme == "martin") {
        if (cdr_scheme == "martin") {
            if (chain == "H")
                current_breakpoints = MARTIN_HEAVY_CDR_BREAKPOINTS;
            else
                current_breakpoints = MARTIN_LIGHT_CDR_BREAKPOINTS;
        } else if (cdr_scheme == "kabat") {
            if (chain == "H")
                current_breakpoints = MARTIN_NMBR_KABAT_H_CDR_BREAKPOINTS;
            else
                current_breakpoints = MARTIN_NMBR_KABAT_L_CDR_BREAKPOINTS;
        } else if (cdr_scheme == "imgt") {
            if (chain == "H")
                current_breakpoints = MARTIN_NMBR_IMGT_H_CDR_BREAKPOINTS;
            else
                current_breakpoints = MARTIN_NMBR_IMGT_L_CDR_BREAKPOINTS;
        } else if (cdr_scheme == "aho") {
            if (chain == "H")
                current_breakpoints = MARTIN_NMBR_AHO_H_CDR_BREAKPOINTS;
            else
                current_breakpoints = MARTIN_NMBR_AHO_L_CDR_BREAKPOINTS;
        } else if (cdr_scheme == "north") {
            if (chain == "H")
                current_breakpoints = MARTIN_NMBR_NORTH_H_CDR_BREAKPOINTS;
            else
                current_breakpoints = MARTIN_NMBR_NORTH_L_CDR_BREAKPOINTS;
        } else {
            throw std::runtime_error(std::string("Unrecognized chain or scheme."));
        }
    }
    else if (scheme == "kabat") {
        if (cdr_scheme == "kabat") {
            if (chain == "H")
                current_breakpoints = KABAT_HEAVY_CDR_BREAKPOINTS;
            else
                current_breakpoints = KABAT_LIGHT_CDR_BREAKPOINTS;
        } else if (cdr_scheme == "martin") {
            if (chain == "H")
                current_breakpoints = KABAT_NMBR_MARTIN_H_CDR_BREAKPOINTS;
            else
                current_breakpoints = KABAT_NMBR_MARTIN_L_CDR_BREAKPOINTS;
        } else if (cdr_scheme == "imgt") {
            if (chain == "H")
                current_breakpoints = KABAT_NMBR_IMGT_H_CDR_BREAKPOINTS;
            else
                current_breakpoints = KABAT_NMBR_IMGT_L_CDR_BREAKPOINTS;
        } else if (cdr_scheme == "aho") {
            if (chain == "H")
                current_breakpoints = KABAT_NMBR_AHO_H_CDR_BREAKPOINTS;
            else
                current_breakpoints = KABAT_NMBR_AHO_L_CDR_BREAKPOINTS;
        } else if (cdr_scheme == "north") {
            if (chain == "H")
                current_breakpoints = KABAT_NMBR_NORTH_H_CDR_BREAKPOINTS;
            else
                current_breakpoints = KABAT_NMBR_NORTH_L_CDR_BREAKPOINTS;
        } else {
            throw std::runtime_error(std::string("Unrecognized chain or "
                        "scheme supplied."));
        }
    }
    else if (scheme == "aho") {
        if (cdr_scheme == "aho") {
            current_breakpoints = AHO_CDR_BREAKPOINTS;
        } else if (cdr_scheme == "imgt") {
            current_breakpoints = AHO_NMBR_IMGT_CDR_BREAKPOINTS;
        } else if (cdr_scheme == "kabat") {
            if (chain == "H")
                current_breakpoints = AHO_NMBR_KABAT_H_CDR_BREAKPOINTS;
            else
                current_breakpoints = AHO_NMBR_KABAT_L_CDR_BREAKPOINTS;
        } else if (cdr_scheme == "martin") {
            if (chain == "H")
                current_breakpoints = AHO_NMBR_MARTIN_H_CDR_BREAKPOINTS;
            else
                current_breakpoints = AHO_NMBR_MARTIN_L_CDR_BREAKPOINTS;
        } else if (cdr_scheme == "north") {
            if (chain == "H")
                current_breakpoints = AHO_NMBR_NORTH_H_CDR_BREAKPOINTS;
            else
                current_breakpoints = AHO_NMBR_NORTH_L_CDR_BREAKPOINTS;
        } else {
            throw std::runtime_error(std::string("Unrecognized chain or "
                        "scheme supplied."));
        }
    }

    int numeric_portion;
    std::array<std::string, 7> cdr_region_labels = {"fmwk1", "cdr1",
            "fmwk2", "cdr2", "fmwk3", "cdr3", "fmwk4"};
    size_t current_token = 0;
    std::string current_label = cdr_region_labels.at(current_token);
    int next_breakpoint = current_breakpoints.at(current_token);


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
