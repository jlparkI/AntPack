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
// Technically light chains also have 128, but due to another weird quirk
// position 128 is never used for light chains.
#define NUM_HEAVY_IMGT_POSITIONS 128
#define NUM_LIGHT_IMGT_POSITIONS 127

// We have 22 AAs -- the 20 standard aminos then a gap in query penalty
// vs gap in template -- thus, 22 expected values.
#define NUM_AAS 22

// Codes for the pathways that can link a score
// to the best-scoring parent.
#define LEFT_TRANSFER 0
#define DIAGONAL_TRANSFER 1
#define UP_TRANSFER 2

// These are "magic number" positions in the IMGT framework at
// which "forwards-backwards" insertion numbering must be applied.
// This is a nuisance, but is out of our control -- the IMGT #ing
// system has this quirk built-in... Note that because IMGT numbers
// from 1, these positions are the actual IMGT position - 1.
#define CDR1_INSERTION_PT 32
#define CDR2_INSERTION_PT 60
#define CDR3_INSERTION_PT 110

// The columns of the score matrix that are accessed for gap penalties.
#define QUERY_GAP_COLUMN 20
#define TEMPLATE_GAP_COLUMN 21

// Highly conserved positions in the IMGT scheme. These are the IMGT #s - 1.
#define HIGHLY_CONSERVED_POSITION_1 22
#define HIGHLY_CONSERVED_POSITION_2 40
#define HIGHLY_CONSERVED_POSITION_3 103
#define HIGHLY_CONSERVED_POSITION_4 117
#define HIGHLY_CONSERVED_POSITION_5 118
#define HIGHLY_CONSERVED_POSITION_6 120

// A default gap penalty for gaps at the beginning and end of the sequence.
#define DEFAULT_GAP_PENALTY -1

// Codes for heavy and light chains.
#define HEAVY_CHAIN 0
#define LIGHT_CHAIN 1


class IMGTAligner {
    public:
        IMGTAligner(py::array_t<double> scoreArray,
                std::vector<std::vector<std::string>> consensus,
                std::string chainName);

        std::tuple<std::vector<std::string>, double,
                std::string, std::string> align(std::string query_sequence);

    protected:
        void fillNeedleScoringTable(double *needleScores, int *pathTrace,
                    int querySeqLen, int rowSize, int *queryAsIdx);

        int numPositions;
        int numRestrictedPositions;
        py::array_t<double> scoreArray;
        const std::string chainName;
        std::vector<std::set<char>> consensusMap;
        std::array<std::string, 6> errorCodeToMessage {"",
                "The sequence contains invalid characters",
                "A fatal runtime error occurred in IMGTAligner. This is very unusual. "
                        "Please report",
                "More than 72 insertions were found. This is highly unlikely for an "
                    "antibody and suggests a problem with this sequence",
                "The alignment length does not match the length of the input sequence. "
                    "This is an unusual error; please report.",
                "The sequence does not have an expected amino acid at one or more "
                    "conserved positions 23, 41, 104, 118, 119, 121. Either it is "
                    "not an antibody, OR there was a serious alignment error, "
                    "OR it contains a large deletion."};

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