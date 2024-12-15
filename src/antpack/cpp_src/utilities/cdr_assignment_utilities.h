#ifndef CDR_ASSIGNMENT_UTILITIES_HEADER_H
#define CDR_ASSIGNMENT_UTILITIES_HEADER_H

// C++ headers
#include <vector>
#include <tuple>
#include <string>

// Library headers

// Project headers
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
        const std::string &scheme);

}  // namespace SequenceUtilities


#endif
