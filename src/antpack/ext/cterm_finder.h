#ifndef C_TERMINAL_FINDER_HEADER_H
#define C_TERMINAL_FINDER_HEADER_H

#include <pybind11/numpy.h>
#include <vector>
#include <string>
#include "utilities.h"
#include "numbering_constants.h"



namespace py = pybind11;






class CTermFinder {
    public:
        CTermFinder(py::array_t<double, py::array::c_style> scoreArray);
        std::string find_c_terminals(std::string query_sequence,
                py::array_t<double, py::array::c_style> bestScores,
                py::array_t<int32_t, py::array::c_style> bestPositions);
        // Allowed error codes. These will be mapped to strings which explain in more detail.
        enum cTermAllowedErrorCodes {noError = 0, invalidSequence = 1};


    protected:

        int numPositions;
        py::array_t<double, py::array::c_style> scoreArray;

        std::array<std::string, 2> errorCodeToMessage {{"",
                "Sequence contains invalid characters"}};

};

#endif
