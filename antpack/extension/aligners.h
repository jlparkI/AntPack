#ifndef ALIGNERS_HEADER_H
#define ALIGNERS_HEADER_H

#include <pybind11/numpy.h>
#include <string>
#include "needle.h"

namespace py = pybind11;

class IMGTAligner {
    public:
        IMGTAligner(
                py::array_t<double> scoreArray
                );

        std::vector<std::string> align(std::string query_sequence);

    protected:
        py::array_t<double> scoreArray;
        int numPositions;
        const int numAAs=21;
};

#endif