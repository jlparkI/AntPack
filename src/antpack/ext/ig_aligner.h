#ifndef IG_ALIGNERS_HEADER_H
#define IG_ALIGNERS_HEADER_H

#include <pybind11/numpy.h>
#include <vector>
#include <string>
#include <algorithm>
#include <tuple>
#include <set>
#include <array>
#include "utilities.h"
#include "numbering_constants.h"

#include <iostream>


namespace py = pybind11;

// We have 22 AAs -- the 20 standard aminos then a gap in query penalty
// vs gap in template -- thus, 22 expected values.
#define NUM_AAS 22

// Codes for the pathways that can link a score
// to the best-scoring parent.
#define LEFT_TRANSFER 0
#define DIAGONAL_TRANSFER 1
#define UP_TRANSFER 2

// The columns of the score matrix that are accessed for gap penalties.
#define QUERY_GAP_COLUMN 20
#define TEMPLATE_GAP_COLUMN 21





class IGAligner {
    public:
        IGAligner(py::array_t<double, py::array::c_style> scoreArray,
                std::vector<std::vector<std::string>> consensus,
                std::string chainName, std::string scheme,
                double terminalTemplateGapPenalty,
                double CterminalQueryGapPenalty,
                bool compressInitialGaps);
        std::tuple<std::vector<std::string>, double, std::string,
                std::string, std::vector<std::string>> align_test_only(std::string query_sequence,
                bool retrieve_cdr_labeling, py::array_t<double> scoreMatrix,
                py::array_t<uint8_t> pathTrace);
        std::tuple<std::vector<std::string>, double, std::string,
                std::string, std::vector<std::string>> align(std::string query_sequence,
                bool retrieve_cdr_labeling);
        // Allowed error codes. These will be mapped to strings which explain in more detail.
        enum allowedErrorCodes {noError = 0, invalidSequence = 1, fatalRuntimeError = 2,
            tooManyInsertions  = 3, alignmentWrongLength = 4,
            unacceptableConservedPositions = 5};


    protected:
        void fillNeedleScoringTable(double *needleScores, uint8_t *pathTrace,
                    int querySeqLen, int rowSize, int *queryAsIdx);
        double core_align_test_only(std::string const &query_sequence,
                bool retrieve_cdr_labeling, std::vector<std::string> &finalNumbering,
                std::vector<std::string> &cdrLabeling,
                allowedErrorCodes &errorCode,
                py::array_t<double> scoreMatrix,
                py::array_t<uint8_t> pathTrace);

        // Default gap penalties for gaps at the beginning and end of the sequence.
        // template is a weak penalty for placing gaps outside the numbering,
        // while query is a weak penalty for placing gaps in the numbering
        // (i.e. skipping numbers).

        int numPositions;
        int numRestrictedPositions;
        py::array_t<double, py::array::c_style> scoreArray;
        const std::string chainName;
        std::string scheme;
        double terminalTemplateGapPenalty;
        double CterminalQueryGapPenalty;
        bool compressInitialGaps;

        std::vector<std::set<char>> consensusMap;
        std::vector<int> highlyConservedPositions;
        std::array<std::string, 6> errorCodeToMessage {{"",
                "Sequence contains invalid characters",
                "Fatal runtime error in IGAligner. Unusual. Please report",
                "> 72 insertions. Suggests a problem with this sequence",
                "Alignment length != length of input sequence. Unusual. Please report.",
                "Unexpected AA at conserved position."}};

        std::array<std::string, 7> cdrRegionLabels {{"fmwk1", "cdr1", "fmwk2", "cdr2",
                                "fmwk3", "cdr3", "fmwk4"}};

        std::vector<int> cdrBreakpoints;

        // Alphabet for numbering insertions. Our preference would be to number insertions
        // as _1, _2 etc, but most numbering programs use letters, so we do the same here
        // for consistency and ease of comparison. This does limit the number of possible
        // insertions BUT if there is ever a case where more than 78 insertions have occurred
        // at some point in the sequence, there is almost certainly an error in the alignment, and
        // that's how we handle this eventuality if it occurs. (Indeed, it is very rare to have
        // more than a few at any given position, so the size of alphabet here is overkill.)
        const std::vector<std::string> alphabet = {"A", "B", "C", "D", "E", "F", "G", "H", "I",
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
