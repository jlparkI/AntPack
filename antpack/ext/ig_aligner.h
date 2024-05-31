#ifndef IG_ALIGNERS_HEADER_H
#define IG_ALIGNERS_HEADER_H

//#include <pybind11/numpy.h>
#include <vector>
#include <string>
#include <algorithm>
#include <tuple>
#include <set>
#include <array>
#include "utilities.h"


#include <iostream>


namespace py = pybind11;

// The minimum number of amino acids in a sequence to try to align it.
// Less than this and it will be immediately rejected. This is fairly
// arbitrary, we didn't put much thought into the selection of 25 --
// a typical chain is > 100 AAs, so anything MUCH less than that is
// clearly a fragment that probably can't be reliably numbered.
#define MINIMUM_SEQUENCE_LENGTH 25


// The IMGT numbering system will always have (at least) 128 positions.
// Technically light chains also have 128, but due to another weird quirk
// position 128 is never used for light chains.
#define NUM_HEAVY_IMGT_POSITIONS 128
#define NUM_LIGHT_IMGT_POSITIONS 127

// Expected number of positions for Martin, Kabat.
#define NUM_HEAVY_MARTIN_KABAT_POSITIONS 113
#define NUM_LIGHT_MARTIN_KABAT_POSITIONS 107

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
#define HIGHLY_CONSERVED_IMGT_1 22
#define HIGHLY_CONSERVED_IMGT_2 40
#define HIGHLY_CONSERVED_IMGT_3 103
#define HIGHLY_CONSERVED_IMGT_4 117
#define HIGHLY_CONSERVED_IMGT_5 118
#define HIGHLY_CONSERVED_IMGT_6 120


// Highly conserved positions in the Martin / Kabat schemes for heavy chains.
// These are the #s - 1.
#define HIGHLY_CONSERVED_KABAT_HEAVY_1 21
#define HIGHLY_CONSERVED_KABAT_HEAVY_2 35
#define HIGHLY_CONSERVED_KABAT_HEAVY_3 91
#define HIGHLY_CONSERVED_KABAT_HEAVY_4 102
#define HIGHLY_CONSERVED_KABAT_HEAVY_5 103
#define HIGHLY_CONSERVED_KABAT_HEAVY_6 105

// Highly conserved positions in the Martin / Kabat schemes for light chains.
// These are the #s - 1.
#define HIGHLY_CONSERVED_KABAT_LIGHT_1 22
#define HIGHLY_CONSERVED_KABAT_LIGHT_2 34
#define HIGHLY_CONSERVED_KABAT_LIGHT_3 87
#define HIGHLY_CONSERVED_KABAT_LIGHT_4 97
#define HIGHLY_CONSERVED_KABAT_LIGHT_5 98
#define HIGHLY_CONSERVED_KABAT_LIGHT_6 100


// The start and end of CDRs in IMGT (ignoring insertions).
// These are the IMGTs - 1. The breakpoints are set on the last
// amino acid before the start of the CDR and the last amino acid
// in the CDR.
#define IMGT_CDR_BREAKPOINT_1 25
#define IMGT_CDR_BREAKPOINT_2 37
#define IMGT_CDR_BREAKPOINT_3 54
#define IMGT_CDR_BREAKPOINT_4 64
#define IMGT_CDR_BREAKPOINT_5 103
#define IMGT_CDR_BREAKPOINT_6 116


// The start and end of CDRs in Kabat (ignoring insertions).
// These are the #s - 1. The breakpoints are set on the last amino
// acid before the start of the CDR and the last amino acid in the
// cdr.
#define KABAT_HEAVY_CDR_BREAKPOINT_1 29
#define KABAT_HEAVY_CDR_BREAKPOINT_2 34
#define KABAT_HEAVY_CDR_BREAKPOINT_3 48
#define KABAT_HEAVY_CDR_BREAKPOINT_4 64
#define KABAT_HEAVY_CDR_BREAKPOINT_5 93
#define KABAT_HEAVY_CDR_BREAKPOINT_6 101

#define KABAT_LIGHT_CDR_BREAKPOINT_1 22
#define KABAT_LIGHT_CDR_BREAKPOINT_2 33
#define KABAT_LIGHT_CDR_BREAKPOINT_3 48
#define KABAT_LIGHT_CDR_BREAKPOINT_4 55
#define KABAT_LIGHT_CDR_BREAKPOINT_5 87
#define KABAT_LIGHT_CDR_BREAKPOINT_6 96

// The start and end of CDRs in Martin (ignoring insertions).
// These are the #s - 1. The breakpoints are set on the last amino
// acid before the start of the CDR and the last amino acid in the
// cdr.
#define MARTIN_HEAVY_CDR_BREAKPOINT_1 24
#define MARTIN_HEAVY_CDR_BREAKPOINT_2 31
#define MARTIN_HEAVY_CDR_BREAKPOINT_3 50
#define MARTIN_HEAVY_CDR_BREAKPOINT_4 55
#define MARTIN_HEAVY_CDR_BREAKPOINT_5 93
#define MARTIN_HEAVY_CDR_BREAKPOINT_6 101

#define MARTIN_LIGHT_CDR_BREAKPOINT_1 24
#define MARTIN_LIGHT_CDR_BREAKPOINT_2 31
#define MARTIN_LIGHT_CDR_BREAKPOINT_3 48
#define MARTIN_LIGHT_CDR_BREAKPOINT_4 51
#define MARTIN_LIGHT_CDR_BREAKPOINT_5 89
#define MARTIN_LIGHT_CDR_BREAKPOINT_6 95




class IGAligner {
    public:
        IGAligner(py::array_t<double> scoreArray,
                std::vector<std::vector<std::string>> consensus,
                std::string chainName, std::string scheme,
                double terminalTemplateGapPenalty,
                double CterminalQueryGapPenalty,
                bool compressInitialGaps);

        std::tuple<std::vector<std::string>, double,
           std::string, std::string, std::vector<std::string>> align(std::string query_sequence,
                        bool retrieve_cdr_labeling);

    protected:
        void fillNeedleScoringTable(double *needleScores, int *pathTrace,
                    int querySeqLen, int rowSize, int *queryAsIdx);

        // Default gap penalties for gaps at the beginning and end of the sequence.
        // template is a weak penalty for placing gaps outside the numbering,
        // while query is a weak penalty for placing gaps in the numbering
        // (i.e. skipping numbers).

        int numPositions;
        int numRestrictedPositions;
        py::array_t<double> scoreArray;
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
