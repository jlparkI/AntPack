#ifndef C_TERMINAL_FINDER_HEADER_H
#define C_TERMINAL_FINDER_HEADER_H

// C++ headers
#include <vector>
#include <array>
#include <string>
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




namespace NumberingTools{


    static constexpr int CTERM_FINDER_ERROR = 0;
    static constexpr int CTERM_FINDER_SUCCESS = 1;


    class CTermFinder {
        public:
            CTermFinder(std::string consensus_filepath);

            /// @brief Finds the most likely H, K and L c-terminal regions in the input sequence.
            int find_best_cterminal(std::string &query_sequence,
                    std::array<double, 3> &best_scores, std::array<int, 3> &best_positions);

            /// @brief A wrapper on find_best_cterminal that can be called from Python. Used
            ///        for testing only.
            int pyfind_c_terminals(std::string query_sequence,
                    nb::ndarray<double, nb::shape<3>, nb::device::cpu, nb::c_contig> best_scores,
                    nb::ndarray<int32_t, nb::shape<3>, nb::device::cpu, nb::c_contig> best_positions);

        protected:
            std::unique_ptr<double[]> score_array;
            size_t score_arr_shape[3];
            int num_positions;

    };
}

#endif
