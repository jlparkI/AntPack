#ifndef ANNOTATOR_BASE_CLASS_HEADER_H
#define ANNOTATOR_BASE_CLASS_HEADER_H

// C++ headers
#include <vector>
#include <array>
#include <string>
#include <utility>
#include <unordered_map>

// Project headers
#include "cterm_finder.h"
#include "numbering_constants.inc"
#include "../utilities/utilities.h"




namespace NumberingTools{


// The minimum length a segment must have in order to be worth
// splitting. Anything less than this and there's no point
// trying to split it into subregions.
static constexpr int MINIMUM_SEGMENT_LENGTH = 61;

// Offset on identified c-terminals if breaking up an input sequence.
static constexpr int CTERMINAL_OFFSET_ALIGNMENT_SHIFT = 11;

// The smallest score allowed for a postulated cterminal to
// be recognized as a possible cterminal.
static constexpr int MINIMUM_TOLERABLE_CTERMINAL_SCORE = 80;


class AnnotatorBaseClassCpp {
 public:
        AnnotatorBaseClassCpp(std::string scheme,
                std::string consensus_filepath);

        /// @brief Sorts a list of position codes, respecting the rules of the numbering scheme.
        /// @param position_code_list A list of position codes to sort. Each must be a string.
        /// @return A vector of position codes that have been sorted.
        std::vector<std::string> sort_position_codes(std::vector<std::string> position_code_list);

        /// @brief Takes a list of numbered sequences and converts this to a multiple sequence alignment
        ///        which can be easily viewed or written to file.
        /// @param sequences A vector of strings each of which is one of the numbered sequences.
        /// @param annotations A vector of tuples of the same length as sequences. Each tuple has four
        ///        elements: a vector of numberings of the same length as that string, a percent identity,
        ///        a chain type and an error message. This tuple is returned by the analyze_seq methods of
        ///        child classes.
        /// @return A tuple of two vectors. The first contains a consensus set of position codes. The
        ///         second contains all of the input sequences, gapped to be the same length.
        std::tuple<std::vector<std::string>, std::vector<std::string>> build_msa(std::vector<std::string> sequences,
            std::vector<std::tuple<std::vector<std::string>, double, std::string, std::string>> annotations);

        /// @brief Takes an input sequence and corresponding numbering and trims it so that gaps
        ///        at the c- and n-terminals are removed.
        /// @param sequence A string which is the sequence to align.
        /// @param alignment A tuple containing four elements -- the numbering, percent identity,
        ///        chain type and error message. This tuple is generated / returned by the
        ///        analyze_seq methods of child classes.
        /// @return Returns a tuple with four elements: the trimmed sequence, the trimmed alignment,
        ///         the first position in the input sequence that was not trimmed and the last position
        ///         in the input sequence that was not trimmed.
        std::tuple<std::string, std::vector<std::string>, int, int> trim_alignment(std::string sequence,
            std::tuple<std::vector<std::string>, double, std::string, std::string> alignment);

        /// @brief Assigns cdr labels to an input sequence.
        /// @param alignment A four-element tuple containing the numbering, the percent identity,
        ///        the chain type and the error message. This tuple is generated / returned by the
        ///        analyze_seq methods of child classes.
        /// @param cdr_scheme One of 'aho', 'imgt', 'kabat' or 'martin'.
        /// @return Returns a vector of strings '-', 'fmwk1', 'cdr1', 'fmwk2', 'cdr2' etc.
        ///         to indicate the region to which each position belongs.
        std::vector<std::string> assign_cdr_labels(std::tuple<std::vector<std::string>, 
               double, std::string, std::string> alignment, std::string cdr_scheme = "");

        /// @brief Looks for likely c-terminals and uses these to split the input sequence
        ///        into subregions.
        void split_sequence_into_subregions(std::vector<std::pair<size_t,size_t>>
                    &subregions, std::string &sequence,
                    size_t minimum_region_size = MINIMUM_SEGMENT_LENGTH,
                    size_t maximum_iterations = 2);

 protected:
        std::string scheme;

        std::unique_ptr<CTermFinder> boundary_finder;
        std::unordered_map<std::string, std::vector<int>> cdr_breakpoints;
        const std::array<std::string, 7> cdr_region_labels {{"fmwk1", "cdr1", "fmwk2", "cdr2",
                                "fmwk3", "cdr3", "fmwk4"}};

        /// @brief Splits a subregion into smaller subregions, or indicates that
        ///        this is not possible.
        bool split_subregion(const std::pair<size_t, size_t> &init_subregion,
                size_t &split_point, const std::string &sequence,
                size_t minimum_region_size = MINIMUM_SEGMENT_LENGTH);
};

}  // namespace NumberingTools

#endif
