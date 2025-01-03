#ifndef ANNOTATOR_BASE_CLASS_HEADER_H
#define ANNOTATOR_BASE_CLASS_HEADER_H

// C++ headers
#include <vector>
#include <array>
#include <string>
#include <utility>
#include <memory>
#include <tuple>
#include <unordered_map>

// Project headers
#include "prefiltering_tool.h"
#include "numbering_constants.inc"
#include "../utilities/utilities.h"




namespace NumberingTools {


class AnnotatorBaseClassCpp {
 public:
        AnnotatorBaseClassCpp(std::string scheme,
                std::string consensus_filepath,
                std::unordered_map<std::string, size_t> nterm_kmers);

        /// @brief Sorts a list of position codes, respecting the
        ///        rules of the numbering scheme.
        /// @param position_code_list A list of position codes to
        ///        sort. Each must be a string.
        /// @return A vector of position codes that have been sorted.
        std::vector<std::string> sort_position_codes(std::vector<std::string>
                position_code_list);

        /// @brief Takes a list of numbered sequences and converts this
        ///        to a multiple sequence alignment which can be easily
        ///        viewed or written to file.
        /// @param sequences A vector of strings each of which is one of
        ///        the numbered sequences.
        /// @param annotations A vector of tuples of the same length as
        ///        sequences. Each tuple has four elements: a vector of
        ///        numberings of the same length as that string, a percent
        ///        identity, a chain type and an error message. This tuple
        ///        is returned by the analyze_seq methods of child classes.
        /// @param add_unobserved_positions If True, positions expected
        ///        for a given scheme are added even if not observed for
        ///        any of the input alignments.
        /// @return A tuple of two vectors. The first contains a consensus
        ///         set of position codes. The second contains all of the
        ///         input sequences, gapped to be the same length.
        std::tuple<std::vector<std::string>, std::vector<std::string>>
            build_msa(std::vector<std::string> sequences,
            std::vector<std::tuple<std::vector<std::string>,
            double, std::string, std::string>> annotations,
            bool add_unobserved_positions = false);

        /// @brief Takes an input sequence and corresponding numbering
        ///        and trims it so that gaps at the c- and n-terminals
        ///        are removed.
        /// @param sequence A string which is the sequence to align.
        /// @param alignment A tuple containing four elements -- the
        ///        numbering, percent identity, chain type and error message.
        ///        This tuple is generated / returned by the analyze_seq methods
        ///        of child classes.
        /// @return Returns a tuple with four elements: the trimmed sequence,
        ///         the trimmed alignment,
        ///         the first position in the input sequence that was not
        ///         trimmed and the last position in the input sequence that
        ///         was not trimmed.
        std::tuple<std::string, std::vector<std::string>, int, int>
            trim_alignment(std::string sequence,
            std::tuple<std::vector<std::string>, double, std::string,
            std::string> alignment);

        /// @brief Assigns cdr labels to an input sequence.
        /// @param alignment A four-element tuple containing
        ///        the numbering, the percent identity, the chain
        ///        type and the error message. This tuple is generated
        ///        / returned by the analyze_seq methods of child classes.
        /// @param cdr_scheme One of 'aho', 'imgt', 'kabat' or 'martin'.
        /// @return Returns a vector of strings '-', 'fmwk1', 'cdr1',
        ///         'fmwk2', 'cdr2' etc. to indicate the region to which
        ///         each position belongs.
        std::vector<std::string>
            assign_cdr_labels(std::tuple<std::vector<std::string>,
               double, std::string, std::string> alignment);


 protected:
        std::string scheme;
        std::unique_ptr<PrefilteringRoutines::PrefilteringTool> boundary_finder;
};

}  // namespace NumberingTools

#endif
