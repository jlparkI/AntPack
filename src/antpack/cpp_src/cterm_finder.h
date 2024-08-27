#ifndef C_TERMINAL_FINDER_HEADER_H
#define C_TERMINAL_FINDER_HEADER_H

#include <pybind11/numpy.h>
#include <vector>
#include <string>
#include <filesystem>
#include "consensus_file_utilities.h"
#include "utilities.h"
#include "numbering_constants.h"


#define INVALID_SEQUENCE 0
#define VALID_SEQUENCE 1



namespace py = pybind11;






class CTermFinder {
    public:
        CTermFinder(py::array_t<double, py::array::c_style> score_array);

        int find_c_terminals(std::string query_sequence,
                std::array<double, 3> &best_scores, std::array<int, 3> &best_positions);

    protected:

        int num_positions;
        py::array_t<double, py::array::c_style> score_array;

};

#endif
