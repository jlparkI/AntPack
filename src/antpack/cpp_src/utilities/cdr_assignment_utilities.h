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
        const std::string &scheme);

}  // namespace SequenceUtilities


#endif
