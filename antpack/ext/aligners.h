#ifndef ALIGNERS_HEADER_H
#define ALIGNERS_HEADER_H

#include <pybind11/numpy.h>
#include <vector>
#include <string>
#include <tuple>
#include <set>
#include <array>
#include "utilities.h"

namespace py = pybind11;

// The IMGT numbering system will always have (at least) 128 positions.
#define NUM_IMGT_POSITIONS 128


class IMGTAligner {
    public:
        IMGTAligner(py::array_t<double> scoreArray,
                std::vector<std::vector<std::string>> consensus);

        std::tuple<std::vector<std::string>, double, int> 
                align(std::string query_sequence);
    protected:
        py::array_t<double> scoreArray;
        int numPositions;
        int numRestrictedPositions;
        const int numAAs=21;
        std::array<std::set<char>, NUM_IMGT_POSITIONS> consensusMap;
};

#endif