/* Contains tools needed for liability search.*/

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


std::vector<std::pair<std::pair<int, int>, std::string>>
LiabilitySearchToolCpp::analyze_seq(std::string sequence,
                std::tuple<std::vector<std::string>,
                double, std::string, std::string> alignment,
                std::string scheme) {
    // As always, exceptions thrown here get sent back to Python
    // if this is used within the Nanobind wrapper.

    if (sequence.length() != std::get<0>(alignment).size()) {
        throw std::runtime_error(std::string("Alignment and sequence "
                    "must have the same size."));
    }
    // Note that the cdr_assignment_utilities functions check that the input
    // alignment and scheme are valid. If they are not they will throw an
    // exception, which if this is used within the Nanobind wrapper will get
    // passed back to Python.
    std::vector<std::string> cdr_labeling;
    SequenceUtilities::assign_cdr_labels(alignment, cdr_labeling,
            scheme);

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
