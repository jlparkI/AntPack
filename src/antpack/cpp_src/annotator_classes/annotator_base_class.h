#ifndef ANNOTATOR_BASE_CLASS_HEADER_H
#define ANNOTATOR_BASE_CLASS_HEADER_H

#include <vector>
#include <string>
#include <unordered_map>
#include "cterm_finder.h"
#include "../numbering_constants.h"
#include "../utilities/utilities.h"


// The minimum length a segment must have in order to be worth
// splitting. Anything less than this and there's no point
// trying to split it into subregions.
#define MINIMUM_SEGMENT_LENGTH 61

// Offset on identified c-terminals if breaking up an input sequence.
#define CTERMINAL_OFFSET_ALIGNMENT_SHIFT 11

// The smallest score allowed for a postulated cterminal to
// be recognized as a possible cterminal.
#define MINIMUM_TOLERABLE_CTERMINAL_SCORE 80




class AnnotatorBaseClassCpp {
    public:
        AnnotatorBaseClassCpp(std::string scheme, std::string consensus_filepath);

        std::vector<std::string> sort_position_codes(std::vector<std::string> position_code_list);
        std::tuple<std::vector<std::string>, std::vector<std::string>> build_msa(std::vector<std::string> sequences,
            std::vector<std::tuple<std::vector<std::string>, double, std::string, std::string>> annotations);

        std::tuple<std::string, std::vector<std::string>, int, int> trim_alignment(std::string sequence,
            std::tuple<std::vector<std::string>, double, std::string, std::string> alignment);
        std::vector<std::string> assign_cdr_labels(std::tuple<std::vector<std::string>, 
               double, std::string, std::string> alignment, std::string cdr_scheme = "");

        void split_sequence_into_subregions(std::vector<std::pair<size_t,size_t>>
                    &subregions, std::string &sequence,
                    size_t minimum_region_size = MINIMUM_SEGMENT_LENGTH,
                    size_t maximum_iterations = 2);

    protected:
        std::string scheme;

        std::unique_ptr<CTermFinder> boundary_finder;
        std::unordered_map<std::string, std::vector<int>> cdr_breakpoints;
        const std::array<std::string, 7> cdr_region_labels {"fmwk1", "cdr1", "fmwk2", "cdr2",
                                "fmwk3", "cdr3", "fmwk4"};

        bool split_subregion(const std::pair<size_t, size_t> &init_subregion,
                size_t &split_point, const std::string &sequence,
                size_t minimum_region_size = MINIMUM_SEGMENT_LENGTH);

};

#endif
