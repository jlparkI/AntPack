#ifndef ALIGNERS_HEADER_H
#define ALIGNERS_HEADER_H

#include <pybind11/numpy.h>
#include <string>
#include <tuple>

namespace py = pybind11;

std::vector<std::string> alphabet = {'A', 'B', 'C', 'D', 'E', 'F',
                        'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N',
                        'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V',
                        'W', 'X', 'Y', 'Z'}

class IMGTAligner {
    public:
        IMGTAligner(
                py::array_t<double> scoreArray
                );

        std::tuple<std::vector<std::string>, int> align(std::string query_sequence);

    protected:
        py::array_t<double> scoreArray;
        int numPositions;
        const int numAAs=21;
};

#endif