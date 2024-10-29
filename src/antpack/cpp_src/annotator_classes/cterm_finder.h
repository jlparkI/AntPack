#ifndef C_TERMINAL_FINDER_HEADER_H
#define C_TERMINAL_FINDER_HEADER_H

// C++ headers
#include <vector>
#include <array>
#include <string>
#include <memory>
#include <filesystem>

// Library headers
#include <nanobind/ndarray.h>
#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>

// Project headers
#include "../utilities/consensus_file_utilities.h"
#include "../utilities/utilities.h"
#include "numbering_constants.inc"


namespace nb = nanobind;




namespace NumberingTools {


static constexpr int CTERM_FINDER_ERROR = 0;
static constexpr int CTERM_FINDER_SUCCESS = 1;
// Threshold score below which it is unlikely the retrieved position
// is in fact a cterminal. This is based on analysis
// of a large (>10 million) sequence dataset.
static constexpr double THRESHOLD_SCORE = 33.;


class CTermFinder {
 public:
        explicit CTermFinder(std::string consensus_filepath);


        /// @brief Convenience function to retrieve number of template
        ///        positions.
        int get_num_positions(void);


        /// @brief Convenience function to retrieve the threshold score
        ///        below which it is unlikely the retrieved position
        ///        is in fact a cterminal. This is based on analysis
        ///        of a large (>10 million) sequence dataset.
        double get_threshold_score(void);


        /// @brief Takes the input scores and positions and merges them
        ///        so that the greater of the two light chain scores (K,L)
        ///        is retained.
        void merge_light_chain_scores(const std::array<double, 3> &scores,
            const std::array<int, 3> &positions, std::array<double, 2> &merged_scores,
            std::array<int, 2> &merged_positions);

        /// @brief Finds the most likely H, K and L c-terminal
        ///        regions in the input sequence.
        int find_best_cterminal(std::string &query_sequence,
                std::array<double, 3> &best_scores,
                std::array<int, 3> &best_positions);


        /// @brief A wrapper on find_best_cterminal that can be called from
        ///        Python. Used for testing only.
        int pyfind_c_terminals(std::string query_sequence,
                nb::ndarray<double, nb::shape<3>, nb::device::cpu, nb::c_contig> best_scores,
                nb::ndarray<int32_t, nb::shape<3>, nb::device::cpu, nb::c_contig> best_positions);

 protected:
        std::unique_ptr<double[]> score_array;
        size_t score_arr_shape[3];
        int num_positions;
        std::vector<std::string> boundary_chains = {"H", "K", "L"};
};

}  // namespace NumberingTools

#endif
