// C++ headers
#include <utility>
#include <tuple>
#include <vector>
#include <memory>
#include <string>


// Project headers
#include "single_chain_annotator.h"






namespace NumberingTools {


// Special error code for cases where a rare but specific type
// of alignment error may have occurred.
static const int POSSIBLE_FWGXG_ERROR_ON_ALIGNMENT = 2;


SingleChainAnnotatorCpp::SingleChainAnnotatorCpp(
    std::vector<std::string> chains,
            std::string scheme, bool compress_init_gaps,
            std::string consensus_filepath
):
    AnnotatorBaseClassCpp(scheme, consensus_filepath),
    chains(chains),
    scheme(scheme),
    compress_init_gaps(compress_init_gaps)
{
    // Note that exceptions thrown here go back to Python via
    // PyBind as long as this constructor is used within the wrapper.
    if (chains.size() < 1 || chains.size() > 3) {
        throw std::runtime_error(std::string("There must be at least one chain specified, "
                "and no more than 3."));
    }


    for (auto & chain : this->chains) {
        if (chain != "H" && chain != "K" && chain != "L") {
            throw std::runtime_error(std::string("Valid chains must be "
                        "one of 'H', 'K', 'L'."));
        } else {
            if (scheme == "imgt") {
                this->scoring_tools.push_back(IGAligner(consensus_filepath,
                        chain, scheme,
                        IMGT_DEFAULT_TERMINAL_TEMPLATE_GAP_PENALTY,
                        IMGT_DEFAULT_C_TERMINAL_QUERY_GAP_PENALTY,
                        compress_init_gaps) );
            } else if (scheme == "aho") {
                this->scoring_tools.push_back(IGAligner(consensus_filepath,
                        chain, scheme,
                        AHO_DEFAULT_TERMINAL_TEMPLATE_GAP_PENALTY,
                        AHO_DEFAULT_C_TERMINAL_QUERY_GAP_PENALTY,
                        compress_init_gaps) );
            } else if (scheme == "martin") {
                this->scoring_tools.push_back(IGAligner(consensus_filepath,
                        chain, scheme,
                        MARTIN_DEFAULT_TERMINAL_TEMPLATE_GAP_PENALTY,
                        MARTIN_DEFAULT_C_TERMINAL_QUERY_GAP_PENALTY,
                        compress_init_gaps) );
            } else if (scheme == "kabat") {
                this->scoring_tools.push_back(IGAligner(consensus_filepath,
                        chain, scheme,
                        KABAT_DEFAULT_TERMINAL_TEMPLATE_GAP_PENALTY,
                        KABAT_DEFAULT_C_TERMINAL_QUERY_GAP_PENALTY,
                        compress_init_gaps) );
            } else {
                throw std::runtime_error(std::string("Invalid scheme "
                        "specified. Please use one of 'imgt', 'kabat', "
                        "'martin', 'aho'."));
            }
        }
    }
}



/// SingleChainAnnotatorCpp function which numbers a single input sequence.
std::tuple<std::vector<std::string>, double, std::string,
        std::string> SingleChainAnnotatorCpp::analyze_seq(std::string sequence) {

    double best_identity = -1;
    std::vector<std::string> empty_numbering;
    std::tuple<std::vector<std::string>, double, std::string,
                    std::string> best_result{ empty_numbering,
                        0, "", "Invalid sequence supplied -- nonstandard AAs"};


    if (!SequenceUtilities::validate_x_sequence(sequence))
        return best_result;

    int err_code = this->align_input_subregion(best_result, best_identity,
            sequence);
    if (err_code !=  POSSIBLE_FWGXG_ERROR_ON_ALIGNMENT)
        return best_result;

    // It can (rarely) happen that we have a chain (usually light)
    // where the expected FGxG motif in the J-gene is altered AND
    // there is significant additional sequence beyond the end of
    // the J-gene, e.g. in a paired chain, AND this results in an alignment
    // error. In the event this unusual combination of circumstances
    // occurs, it usually manifests in the form of an unusually long
    // CDR that exceeds the maximum number of allowed insertion codes.

    // If this may have occurred, find the most likely c-terminal
    // region and re-align the portion of the sequence that precedes
    // this. This is not ideal if the sequence contains
    // multiple chains because the other chains will go unreported,
    // but a SingleChainAnnotator by definition is not ideal for
    // this scenario. MultiChainAnnotator or PairedChainAnnotator
    // should be used instead.

    std::array<double, 3> scores;
    std::array<int, 3> positions;
    if (!this->boundary_finder->find_best_cterminal(sequence,
                scores, positions)) {
        empty_numbering.clear();
        best_result = { empty_numbering, 0, "",
            "Alignment error -- no c-terminal region found."};
        return best_result;
    }


    // Now find the most likely cterminal match, use that to cut the
    // sequence and align the region prior to this.
    size_t best_position = 0;
    double best_score = 0;

    for (size_t i = 0; i < scores.size(); i++) {
        if (scores[i] > best_score) {
            best_position = positions[i];
            best_score = scores[i];
        }
    }

    // Add the appropriate number of letters to best_position to
    // use it as a cutoff. Make sure it is within the sequence.
    // Also make sure it is not unreasonably short.
    best_position = best_position + this->boundary_finder->
        get_num_positions() + 5;
    best_position = best_position > MINIMUM_SEQUENCE_LENGTH ?
        best_position : MINIMUM_SEQUENCE_LENGTH;
    best_position = best_position < sequence.length() ? best_position :
        sequence.length();

    // If the best position is the same as the sequence length,
    // realigning will make no difference; return existing
    // result.
    if (best_position == sequence.length())
        return best_result;



    std::string subsequence = sequence.substr(0, best_position);
    best_identity = -1;

    err_code = this->align_input_subregion(best_result, best_identity,
            subsequence);

    for (size_t i = 0; i < (sequence.length() - subsequence.length());
            i++) {
        std::get<0>(best_result).push_back("-");
    }

    return best_result;
}




/// Aligns the input sequence, which may be either the full sequence or a subregion of it.
int SingleChainAnnotatorCpp::align_input_subregion(std::tuple<std::vector<std::string>, double,
                std::string, std::string> &best_result, double &best_identity,
                std::string &query_sequence) {
    auto queryAsIdx = std::make_unique<int[]>(query_sequence.length());
    bool fwgxg_error = false;

    if (!SequenceUtilities::convert_x_sequence_to_array(queryAsIdx.get(),
                query_sequence))
        return INVALID_SEQUENCE;

    for (size_t i=0; i < this->scoring_tools.size(); i++) {
        std::vector<std::string> final_numbering;
        std::string error_message = "";
        double percent_identity = -1;

        final_numbering.reserve(query_sequence.length() * 2);

        this->scoring_tools[i].align(query_sequence,
                    queryAsIdx.get(), final_numbering, percent_identity,
                    error_message);
        if (percent_identity > best_identity) {
            std::get<0>(best_result) = final_numbering;
            std::get<1>(best_result) = percent_identity;
            std::get<2>(best_result) = this->scoring_tools[i].get_chain_name();
            std::get<3>(best_result) = error_message;
            best_identity = percent_identity;
        }
        if (error_message.length() > 2) {
            if (error_message.substr(0, 2) == "> ")
                fwgxg_error = true;
        }
    }
    if (fwgxg_error && best_identity < 0.8)
        return POSSIBLE_FWGXG_ERROR_ON_ALIGNMENT;
    return VALID_SEQUENCE;
}






// SingleChainAnnotatorCpp function which numbers a list of input sequences.
std::vector<std::tuple<std::vector<std::string>, double, std::string,
            std::string>> SingleChainAnnotatorCpp::analyze_seqs(std::vector<std::string> sequences) {
    std::vector<std::tuple<std::vector<std::string>, double, std::string,
            std::string>> output_results;

    for (size_t i=0; i < sequences.size(); i++)
        output_results.push_back(this->analyze_seq(sequences[i]));

    return output_results;
}

}  // namespace NumberingTools
