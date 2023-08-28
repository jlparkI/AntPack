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
        // Alphabet for numbering insertions. Our preference would be to number insertions
        // as _1, _2 etc, but most numbering programs use letters, so we do the same here
        // for consistency and ease of comparison. This does limit the number of possible
        // insertions BUT if there is ever a case where more than 78 insertions have occurred
        // at some point in the sequence, there is almost certainly an error in the alignment, and
        // that's how we handle this eventuality if it occurs. (Indeed, it is very rare to have
        // more than a few at any given position, so the size of alphabet here is overkill.)
        std::vector<std::string> alphabet = {"A", "B", "C", "D", "E", "F", "G", "H", "I",
                                    "J", "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T", "U",
                                    "V", "W", "X", "Y", "Z", "AA", "BB", "CC", "DD", "EE",
                                    "FF", "GG", "HH", "II", "JJ", "KK", "LL", "MM", "NN",
                                    "OO", "PP", "QQ", "RR", "SS", "TT", "UU", "VV",
                                    "WW", "XX", "YY", "ZZ", "AAA", "BBB", "CCC", "DDD", 
                                    "EEE", "FFF", "GGG", "HHH", "III", "JJJ", "KKK", "LLL",
                                    "MMM", "NNN", "OOO", "PPP", "QQQ", "RRR", "SSS", "TTT",
                                    "UUU", "VVV", "WWW", "XXX", "YYY", "ZZZ"};
};

#endif