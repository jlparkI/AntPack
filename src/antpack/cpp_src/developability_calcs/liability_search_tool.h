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
