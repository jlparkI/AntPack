/* Contains the parent class for other annotator classes, which provides methods
 * that all annotator tools should expose to their callers.*/

#include <vector>
#include <string>
#include <unordered_map>
#include <iostream>
#include "numbering_constants.h"
#include "utilities.h"


class AnnotatorBaseClassCpp {
    public:
        AnnotatorBaseClassCpp(std::string scheme);

        std::vector<std::string> sort_position_codes(std::vector<std::string> position_code_list);
        std::tuple<std::vector<std::string>, std::vector<std::string>> build_msa(std::vector<std::string> sequences,
            std::vector<std::tuple<std::vector<std::string>, double, std::string, std::string>> annotations);
        std::tuple<std::string, std::vector<std::string>, int, int> trim_alignment(std::string sequence,
            std::tuple<std::vector<std::string>, double, std::string, std::string> alignment);
        std::vector<std::string> assign_cdr_labels(std::tuple<std::vector<std::string>, 
               double, std::string, std::string> alignment, std::string cdr_scheme = "");


    protected:
        std::string scheme;
        std::unordered_map<std::string, std::vector<int>> cdr_breakpoints;
        const std::array<std::string, 7> cdr_region_labels {{"fmwk1", "cdr1", "fmwk2", "cdr2",
                                "fmwk3", "cdr3", "fmwk4"}};

};
