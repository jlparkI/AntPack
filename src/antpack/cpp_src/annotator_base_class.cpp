/* Contains the parent class for other annotator classes, which provides methods
 * that all annotator tools should expose to their callers.*/
#include "annotator_base_class.h"



AnnotatorBaseClassCpp::AnnotatorBaseClassCpp(std::string scheme):
    scheme(scheme)
{

    // Set up a list of breakpoints that mark the dividing lines between
    // framework and CDR regions for each possible scheme. This way the
    // user can extract CDRs for any scheme, not just the currently selected
    // one.
    this->cdr_breakpoints["imgt_H"] = {IMGT_CDR_BREAKPOINT_1, IMGT_CDR_BREAKPOINT_2,
                        IMGT_CDR_BREAKPOINT_3, IMGT_CDR_BREAKPOINT_4,
                        IMGT_CDR_BREAKPOINT_5, IMGT_CDR_BREAKPOINT_6};
    this->cdr_breakpoints["imgt_L"] = {IMGT_CDR_BREAKPOINT_1, IMGT_CDR_BREAKPOINT_2,
                        IMGT_CDR_BREAKPOINT_3, IMGT_CDR_BREAKPOINT_4,
                        IMGT_CDR_BREAKPOINT_5, IMGT_CDR_BREAKPOINT_6};
    this->cdr_breakpoints["kabat_L"] = {KABAT_LIGHT_CDR_BREAKPOINT_1,
                        KABAT_LIGHT_CDR_BREAKPOINT_2, KABAT_LIGHT_CDR_BREAKPOINT_3,
                        KABAT_LIGHT_CDR_BREAKPOINT_4, KABAT_LIGHT_CDR_BREAKPOINT_5,
                        KABAT_LIGHT_CDR_BREAKPOINT_6};
    this->cdr_breakpoints["kabat_H"] = {KABAT_HEAVY_CDR_BREAKPOINT_1,
                        KABAT_HEAVY_CDR_BREAKPOINT_2, KABAT_HEAVY_CDR_BREAKPOINT_3,
                        KABAT_HEAVY_CDR_BREAKPOINT_4, KABAT_HEAVY_CDR_BREAKPOINT_5,
                        KABAT_HEAVY_CDR_BREAKPOINT_6};
    this->cdr_breakpoints["martin_L"] = {MARTIN_LIGHT_CDR_BREAKPOINT_1,
                        MARTIN_LIGHT_CDR_BREAKPOINT_2, MARTIN_LIGHT_CDR_BREAKPOINT_3,
                        MARTIN_LIGHT_CDR_BREAKPOINT_4, MARTIN_LIGHT_CDR_BREAKPOINT_5,
                        MARTIN_LIGHT_CDR_BREAKPOINT_6};
    this->cdr_breakpoints["martin_H"] = {MARTIN_HEAVY_CDR_BREAKPOINT_1,
                        MARTIN_HEAVY_CDR_BREAKPOINT_2, MARTIN_HEAVY_CDR_BREAKPOINT_3,
                        MARTIN_HEAVY_CDR_BREAKPOINT_4, MARTIN_HEAVY_CDR_BREAKPOINT_5,
                        MARTIN_HEAVY_CDR_BREAKPOINT_6};
}





// AnnotatorBaseClassCpp function which sorts a list of position codes (essentially
// by calling the corresponding function in utilities). Essentially a wrapper on the
// corresponding utility function that can be accessed by Python or outside callers.
std::vector<std::string> AnnotatorBaseClassCpp::sort_position_codes(std::vector<std::string> position_code_list){
    std::vector<std::string> cleanedCodes, sortedCodes;
    int isValid;

    for (auto &code : position_code_list){
        if (code != "-")
            cleanedCodes.push_back(code);
    }
    isValid = sort_position_codes_utility(cleanedCodes, this->scheme,
            sortedCodes);
    if (isValid == 0){
        throw std::runtime_error(std::string("One or more of the supplied position "
                    "codes is invalid given the specified scheme."));
    }
    return sortedCodes;
}



// AnnotatorBaseClassCpp function which builds a multiple sequence alignment from
// a list of previously numbered sequences. The sequences must all be of the
// same chain type. Essentially a wrapper on the corresponding utility function
// that can be accessed by Python or outside callers.
std::tuple<std::vector<std::string>, std::vector<std::string>> AnnotatorBaseClassCpp::build_msa(std::vector<std::string> sequences,
        std::vector<std::tuple<std::vector<std::string>, double, std::string, std::string>> annotations){
    std::vector<std::string> positionCodes, alignedSeqs;
    int errCode;

    errCode = build_msa_utility(sequences, annotations, positionCodes,
            alignedSeqs, this->scheme);
    if (errCode != 1)
        throw std::runtime_error(std::string("A fatal error occured when building an MSA -- please report."));
    return std::tuple<std::vector<std::string>, std::vector<std::string>>{positionCodes, alignedSeqs};
}


// AnnotatorBaseClassCpp function which trims an alignment to remove gap positions at
// either end. Essentially a wrapper on the corresponding utility function that can
// be accessed by Python or outside callers.
std::tuple<std::string, std::vector<std::string>, int, int> AnnotatorBaseClassCpp::trim_alignment(std::string sequence,
        std::tuple<std::vector<std::string>, double, std::string, std::string> alignment){
    int exstart, exend;
    std::vector<char> trimmedSeq;
    std::vector<std::string> trimmedAlignment;

    int errCode = trim_alignment_utility(sequence, alignment, trimmedAlignment,
            exstart, exend, trimmedSeq);
    if (errCode != 1)
        throw std::runtime_error(std::string("Invalid sequence / alignment pairing supplied."));

    std::string trimmedSeqStr(trimmedSeq.begin(), trimmedSeq.end());

    return std::tuple<std::string, std::vector<std::string>, int, int>{trimmedSeqStr,
        trimmedAlignment, exstart, exend};
}


// AnnotatorBaseClassCpp function that assigns cdr and framework labels based on an
// alignment performed by analyze_seq (or similar). It can use a scheme other than
// the one used by the object, so that the user can assign
// cdr labels for an IMGT-numbered scheme using Kabat CDR definitions if desired.
std::vector<std::string> AnnotatorBaseClassCpp::assign_cdr_labels(std::tuple<std::vector<std::string>, 
               double, std::string, std::string> alignment, std::string cdr_scheme){
    std::string *current_scheme;
    std::vector<int> *current_breakpoints;
    std::vector<std::string> cdr_labeling;
    int numeric_portion;
    std::string current_label;
    size_t current_token = 0;
    int next_breakpoint;

    if (cdr_scheme == "")
        current_scheme = &this->scheme;
    else if (cdr_scheme == "imgt" || cdr_scheme == "kabat" || cdr_scheme == "martin" ||
            cdr_scheme == "aho")
        current_scheme = &cdr_scheme;
    else{
        throw std::runtime_error(std::string("Unrecognized scheme supplied. Use '' to use "
                    "the scheme currently selected for this Annotator, or use a currently "
                    "accepted scheme."));
    }

    if (*current_scheme == "imgt"){
        if (std::get<2>(alignment) == "H")
            current_breakpoints = &this->cdr_breakpoints.at("imgt_H");
        else if (std::get<2>(alignment) == "L" || std::get<2>(alignment) == "K")
            current_breakpoints = &this->cdr_breakpoints.at("imgt_L");
        else
            throw std::runtime_error(std::string("Unrecognized chain or scheme supplied."));
    }
    else if (*current_scheme == "martin"){
        if (std::get<2>(alignment) == "H")
            current_breakpoints = &this->cdr_breakpoints.at("martin_H");
        else if (std::get<2>(alignment) == "L" || std::get<2>(alignment) == "K")
            current_breakpoints = &this->cdr_breakpoints.at("martin_L");
        else
            throw std::runtime_error(std::string("Unrecognized chain or scheme supplied."));
    }
    else if (*current_scheme == "kabat"){
        if (std::get<2>(alignment) == "H")
            current_breakpoints = &this->cdr_breakpoints.at("kabat_H");
        else if (std::get<2>(alignment) == "L" || std::get<2>(alignment) == "K")
            current_breakpoints = &this->cdr_breakpoints.at("kabat_L");
        else
            throw std::runtime_error(std::string("Unrecognized chain or scheme supplied."));
    }
    else if (*current_scheme == "aho"){
        if (std::get<2>(alignment) == "H")
            current_breakpoints = &this->cdr_breakpoints.at("aho_H");
        else if (std::get<2>(alignment) == "L" || std::get<2>(alignment) == "K")
            current_breakpoints = &this->cdr_breakpoints.at("aho_L");
        else
            throw std::runtime_error(std::string("Unrecognized chain or scheme supplied."));
    }

    next_breakpoint = current_breakpoints->at(current_token);
    current_label = this->cdr_region_labels[current_token];

    for (size_t i=0; i < std::get<0>(alignment).size(); i++){
        if (std::get<0>(alignment).at(i) == "-"){
            cdr_labeling.push_back("-");
            continue;
        }
        // This will throw if the string starts with a letter but will otherwise extract the integer piece.
        // AntPack never places a letter at the start of the code, so this will not happen unless
        // the user has passed some altered / corrupted input.
        try{
            numeric_portion = std::stoi(std::get<0>(alignment)[i]);
        }
        catch (...){
            throw std::runtime_error(std::string("An invalid position code was passed. The alignment "
                        "passed to this function should be unaltered output from analyze_seq and not "
                        "some other procedure."));
        }
        if (numeric_portion >= next_breakpoint){
            if (current_token < (current_breakpoints->size() - 1) ){
                while (current_token < (current_breakpoints->size() - 1) &&
                        numeric_portion >= next_breakpoint){
                    current_token += 1;
                    next_breakpoint = current_breakpoints->at(current_token);
                    current_label = this->cdr_region_labels[current_token];
                }
                if (numeric_portion >= next_breakpoint){
                    // Set next_breakpoint to an arbitrarily high, unachievable number.
                    next_breakpoint = 10000;
                    current_token += 1;
                    current_label = this->cdr_region_labels[this->cdr_region_labels.size() - 1];
                }
            }
            else{
                // Set next_breakpoint to an arbitrarily high, unachievable number.
                next_breakpoint = 10000;
                current_token += 1;
                current_label = this->cdr_region_labels[this->cdr_region_labels.size() - 1];
            }
        }
        cdr_labeling.push_back(current_label);
    }

    return cdr_labeling;
}
