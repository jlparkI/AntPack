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



/// @brief Assigns cdr labels to a previously constructed alignment.
/// @param alignment An alignment output by SingleChainAnnotator or
///        PairedChainAnnotator.
/// @param cdr_labels A vector in which the output will be stored.
/// @param scheme The scheme -- one of 'imgt', 'aho', 'martin', 'kabat'.
void assign_cdr_labels(const std::tuple<std::vector<std::string>,
        double, std::string, std::string> &alignment,
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
        if (std::get<2>(alignment) == "H" || std::get<2>(alignment) == "L" ||
                std::get<2>(alignment) == "K")
            current_breakpoints = NumberingTools::IMGT_CDR_BREAKPOINTS;
        else
            throw std::runtime_error(std::string("Unrecognized chain or "
                        "scheme supplied."));
    } else if (scheme == "martin") {
        if (std::get<2>(alignment) == "H")
            current_breakpoints = NumberingTools::MARTIN_HEAVY_CDR_BREAKPOINTS;
        else if (std::get<2>(alignment) == "L" || std::get<2>(alignment) == "K")
            current_breakpoints = NumberingTools::MARTIN_LIGHT_CDR_BREAKPOINTS;
        else
            throw std::runtime_error(std::string("Unrecognized chain or "
                        "scheme supplied."));
    } else if (scheme == "kabat") {
        if (std::get<2>(alignment) == "H")
            current_breakpoints = NumberingTools::KABAT_HEAVY_CDR_BREAKPOINTS;
        else if (std::get<2>(alignment) == "L" || std::get<2>(alignment) == "K")
            current_breakpoints = NumberingTools::KABAT_LIGHT_CDR_BREAKPOINTS;
        else
            throw std::runtime_error(std::string("Unrecognized chain or "
                        "scheme supplied."));
    } else if (scheme == "aho") {
        if (std::get<2>(alignment) == "H" || std::get<2>(alignment) == "L" ||
                std::get<2>(alignment) == "K")
            current_breakpoints = NumberingTools::AHO_CDR_BREAKPOINTS;
        else
            throw std::runtime_error(std::string("Unrecognized chain or "
                        "scheme supplied."));
    }

    next_breakpoint = current_breakpoints.at(current_token);
    current_label = cdr_region_labels[current_token];

    for (size_t i=0; i < std::get<0>(alignment).size(); i++) {
        if (std::get<0>(alignment).at(i) == "-") {
            cdr_labeling.push_back("-");
            continue;
        }
        // This will throw if the string starts with a letter but will
        // otherwise extract the integer piece.
        // AntPack never places a letter at the start of the code,
        // so this will not happen unless
        // the user has passed some altered / corrupted input.
        try {
            numeric_portion = std::stoi(std::get<0>(alignment)[i]);
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
