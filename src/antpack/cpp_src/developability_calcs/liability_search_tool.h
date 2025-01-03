#ifndef LIABILITY_SEARCH_TOOL_HEADER_H
#define LIABILITY_SEARCH_TOOL_HEADER_H

// C++ headers
#include <vector>
#include <tuple>
#include <string>

// Library headers


// Project headers


namespace LiabilitySearch {

class LiabilitySearchToolCpp {
 public:
        LiabilitySearchToolCpp();

        /// @brief Finds potential liabilities in the sequence.
        /// @param alignment A tuple containing the output of AntPack's single
        ///                  or paired chain annotator.
        /// @param sequence The sequence as a string. Should not contain gaps.
        /// @param scheme One of 'imgt', 'aho', 'martin', 'kabat'.
        /// @return A vector of tuples. Each tuple contains a two-tuple of ints
        ///         (the start and end of the liability) and a string describing
        ///         the type of liability.
        std::vector<std::pair<std::pair<int, int>, std::string>>
            analyze_seq(std::string sequence,
                std::tuple<std::vector<std::string>,
                double, std::string, std::string> alignment,
                std::string scheme);
};

}  // namespace LiabilitySearch


#endif
