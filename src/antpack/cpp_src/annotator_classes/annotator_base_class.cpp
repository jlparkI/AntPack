/* Contains the parent class for other annotator classes, which provides methods
 * that all annotator tools should expose to their callers.*/
// C++ headers
#include <tuple>
#include <vector>
#include <string>
#include <memory>
#include <unordered_map>

// Project headers
#include "annotator_base_class.h"
#include "../utilities/cdr_assignment_utilities.h"

namespace NumberingTools {


AnnotatorBaseClassCpp::AnnotatorBaseClassCpp(std::string scheme,
        std::string consensus_filepath,
        std::unordered_map<std::string, size_t> nterm_kmers):
    scheme(scheme) {
    this->boundary_finder = std::make_unique
        <PrefilteringRoutines::PrefilteringTool>(consensus_filepath,
                nterm_kmers);
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
    if (!SequenceUtilities::sort_position_codes_utility(cleanedCodes, this->scheme,
            sortedCodes)) {
        throw std::runtime_error(std::string("One or more of the supplied position "
                    "codes is invalid given the specified scheme."));
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
        std::tuple<std::vector<std::string>, 
               double, std::string, std::string> alignment) {

    std::vector<std::string> cdr_labeling;
    SequenceUtilities::assign_cdr_labels(alignment,
            cdr_labeling, this->scheme);
    return cdr_labeling;
}

}  // namespace NumberingTools
