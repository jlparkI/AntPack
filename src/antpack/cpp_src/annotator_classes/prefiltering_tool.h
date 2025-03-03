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
#ifndef C_TERMINAL_FINDER_HEADER_H
#define C_TERMINAL_FINDER_HEADER_H

// C++ headers
#include <vector>
#include <array>
#include <string>
#include <memory>
#include <filesystem>
#include <unordered_map>

// Library headers
#include <nanobind/ndarray.h>
#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>

// Project headers
#include "../utilities/consensus_file_utilities.h"
#include "../utilities/utilities.h"
#include "numbering_constants.inc"


namespace nb = nanobind;




namespace PrefilteringRoutines {


static constexpr int CTERM_FINDER_ERROR = 0;
static constexpr int CTERM_FINDER_SUCCESS = 1;
static constexpr int KMER_SIZE = 9;
static constexpr int KMER_WINDOW_WIDTH = 14;

// Threshold score below which it is unlikely the retrieved position
// is in fact a cterminal or a start. This is based on analysis
// of a large (>10 million) sequence dataset.
static constexpr double CTERM_THRESHOLD_SCORE = 33.;
static constexpr int START_THRESHOLD_SCORE = 4;


class PrefilteringTool {
 public:
        PrefilteringTool(std::string consensus_filepath,
                std::unordered_map<std::string, size_t> nterm_kmers);


        /// @brief Convenience function to retrieve chain indexing used
        ///        by this class.
        std::vector<std::string> get_chain_list(void);

        /// @brief Convenience function to retrieve number of template
        ///        positions.
        int get_num_positions(void);


        /// @brief Convenience function to retrieve the threshold score
        ///        below which it is unlikely the retrieved position
        ///        is in fact a cterminal. This is based on analysis
        ///        of a large (>10 million) sequence dataset.
        double get_cterm_threshold_score(void);

        /// @brief Convenience function to retrieve the threshold score
        ///        below which it is unlikely the retrieved position
        ///        is in fact a start region. This is based on analysis
        ///        of a large sequence dataset.
        int get_nterm_threshold_score(void);


        /// @brief Convenience function to retrieve the kmer window
        ///        size.
        int get_kmer_window_size(void);


        /// @brief Takes the input scores and positions and merges them
        ///        so that the greater of the two light chain scores (K,L)
        ///        is retained.
        void merge_light_chain_cterm_scores(const std::array<double, 3> &scores,
            const std::array<int, 3> &positions, std::array<double, 2> &merged_scores,
            std::array<int, 2> &merged_positions);


        /// @brief Takes the input scores and positions and merges them
        ///        so that the greater of the two light chain scores (K,L)
        ///        is retained.
        void merge_light_chain_nterm_scores(const std::array<int, 3> &scores,
            const std::array<int, 3> &positions, std::array<int, 2> &merged_scores,
            std::array<int, 2> &merged_positions);

        /// @brief Finds the most likely H, K and L c-terminal
        ///        regions in the input sequence.
        int find_start_end_zones(std::string &query_sequence,
            std::array<double, 3> &best_cterm_scores,
            std::array<int, 3> &best_cterm_positions,
            std::array<int, 3> &best_nterm_scores,
            std::array<int, 3> &best_nterm_positions);


        /// @brief A wrapper on find_best_cterminal that can be called from
        ///        Python. Used for testing only.
        int pyfind_c_terminals(std::string query_sequence,
            nb::ndarray<double, nb::shape<3>, nb::device::cpu, nb::c_contig> best_cterm_scores,
            nb::ndarray<int32_t, nb::shape<3>, nb::device::cpu, nb::c_contig> best_cterm_positions,
            nb::ndarray<int32_t, nb::shape<3>, nb::device::cpu, nb::c_contig> best_nterm_scores,
            nb::ndarray<int32_t, nb::shape<3>, nb::device::cpu, nb::c_contig> best_nterm_positions);

 protected:
        std::unordered_map<std::string, size_t> ig_nterm_kmers;
        std::unique_ptr<double[]> score_array;
        size_t score_arr_shape[3];
        int num_positions;
        std::vector<std::string> boundary_chains = {"H", "K", "L"};
};

}  // namespace NumberingTools

#endif
