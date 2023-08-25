#ifndef ALIGNERS_HEADER_H
#define ALIGNERS_HEADER_H

#include <pybind11/numpy.h>
#include <vector>
#include <string>
#include <tuple>

namespace py = pybind11;


class IMGTAligner {
    public:
        IMGTAligner(py::array_t<double> scoreArray);

        std::tuple<std::vector<std::string>, int> align(std::string query_sequence);

    protected:
        py::array_t<double> scoreArray;
        int numPositions;
        const int numAAs=21;
};

#endif