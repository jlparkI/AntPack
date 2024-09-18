#ifndef C_TERMINAL_FINDER_HEADER_H
#define C_TERMINAL_FINDER_HEADER_H

#include <nanobind/ndarray.h>
#include <nanobind/nanobind.h>
#include <vector>
#include <array>
#include <string>
#include <filesystem>
#include "../utilities/consensus_file_utilities.h"
#include "../utilities/utilities.h"
#include "../numbering_constants.h"


#define INVALID_SEQUENCE 0
#define VALID_SEQUENCE 1



namespace nb = nanobind;






class CTermFinder {
    public:
        CTermFinder(std::string consensus_filepath);

        int find_c_terminals(std::string &query_sequence,
                std::array<double, 3> &best_scores, std::array<int, 3> &best_positions);

    protected:
        std::unique_ptr<double[]> score_array;
        size_t score_arr_shape[3];
        int num_positions;

};

#endif
