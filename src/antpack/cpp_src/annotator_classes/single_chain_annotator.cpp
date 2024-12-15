// C++ headers
#include <utility>
#include <tuple>
#include <vector>
#include <memory>
#include <string>
#include <algorithm>


// Project headers
#include "single_chain_annotator.h"






namespace NumberingTools {


// Special error code for cases where a rare but specific type
// of alignment error may have occurred.
static const int POSSIBLE_FWGXG_ERROR_ON_ALIGNMENT = 2;


SingleChainAnnotatorCpp::SingleChainAnnotatorCpp(
    std::vector<std::string> chains,
            std::string scheme, std::string consensus_filepath,
            std::unordered_map<std::string, size_t> nterm_kmers
):
    AnnotatorBaseClassCpp(scheme, consensus_filepath, nterm_kmers),
    chains(chains),
    scheme(scheme)
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
                        chain, scheme) );
            } else if (scheme == "aho") {
                this->scoring_tools.push_back(IGAligner(consensus_filepath,
                        chain, scheme) );
            } else if (scheme == "martin") {
                this->scoring_tools.push_back(IGAligner(consensus_filepath,
                        chain, scheme) );
            } else if (scheme == "kabat") {
                this->scoring_tools.push_back(IGAligner(consensus_filepath,
                        chain, scheme) );
            } else {
                throw std::runtime_error(std::string("Invalid scheme "
                        "specified. Please use one of 'imgt', 'kabat', "
                        "'martin', 'aho'."));
            }
        }
    }
}



/// SingleChainAnnotatorCpp function which numbers a single input sequence.
/// This is essentially a wrapper on align_input_subregion which decides
/// what chain(s) are worth aligning. This can save a considerable amount
/// of time since alignment is the most expensive step. In the event we
/// are unable to determine, we can of course try all possibilities.
std::tuple<std::vector<std::string>, double, std::string,
        std::string> SingleChainAnnotatorCpp::analyze_seq(std::string sequence) {

    std::string preferred_chain = "";
    std::vector<std::string> empty_numbering;
    std::tuple<std::vector<std::string>, double, std::string,
                    std::string> best_result{ empty_numbering,
                        0, "", "Sequence contains invalid characters"};


    if (!SequenceUtilities::validate_x_sequence(sequence))
        return best_result;

    // Check the sequence for c-terminal and nterm regions that
    // clearly belong to a particular chain type. If we find these,
    // we can align to one chain type only and save time. This
    // should never returns an error code, because we already
    // validated the sequence.

    std::array<double, 3> cterm_scores;
    std::array<int, 3> cterm_positions, nterm_positions, nterm_scores;
    if (this->boundary_finder->find_start_end_zones(sequence,
                cterm_scores, cterm_positions,
                nterm_scores, nterm_positions)) {

        // If no error code, find the most plausible cterminal.
        // If the most plausible n-terminal region matches,
        // we know the chain type with high confidence. Make
        // sure the scores meet at least some minimum threshold.
        int best_nterm_score = this->boundary_finder->
            get_nterm_threshold_score();
        double best_cterm_score = this->boundary_finder->
            get_cterm_threshold_score();
        std::vector<std::string> chain_list = this->boundary_finder->
            get_chain_list();

        for (size_t i=0; i < 3; i++) {
            if (cterm_scores[i] > best_cterm_score) {
                best_cterm_score = cterm_scores[i];

                if (nterm_scores[i] > best_nterm_score) {
                    preferred_chain = chain_list.at(i);
                    best_nterm_score = nterm_scores[i];
                } else {
                    preferred_chain = "";
                }
            } else if (nterm_scores[i] > best_nterm_score) {
                best_nterm_score = nterm_scores[i];
            }
        }
    }

    // Perform the alignment, passing the preferred_chain variable.
    // Since preferred_chain is initialized to "", if no preferred
    // chain was found, the aligner will try all available chains.

    int err_code = this->align_input_subregion(best_result, sequence,
            preferred_chain);

    // It is unlikely but possible to have a rare error when the j-region
    // is highly unusual and CDR3 is very long. Return if there is no
    // error code indicating this has happened.
    if (err_code !=  POSSIBLE_FWGXG_ERROR_ON_ALIGNMENT)
        return best_result;

    // If error code indicates this MAY have happened, we will proceed
    // by trying to re-align with the most likely c-terminal we found
    // as a cutoff and aligning the region prior to that cutoff.

    size_t best_seq_cutoff = 0;
    double best_score = 0;

    for (size_t i = 0; i < cterm_scores.size(); i++) {
        if (cterm_scores[i] > best_score) {
            best_seq_cutoff = cterm_positions[i];
            best_score = cterm_scores[i];
        }
    }

    // Add the appropriate number of letters to best_seq_cutoff to
    // use it as a cutoff. Make sure it is within the sequence.
    // Also make sure it is not unreasonably short.
    best_seq_cutoff = best_seq_cutoff + this->boundary_finder->
        get_num_positions() + 5;
    best_seq_cutoff = best_seq_cutoff > MINIMUM_SEQUENCE_LENGTH ?
        best_seq_cutoff : MINIMUM_SEQUENCE_LENGTH;
    best_seq_cutoff = best_seq_cutoff < sequence.length() ? best_seq_cutoff :
        sequence.length();

    std::string subsequence = sequence.substr(0, best_seq_cutoff);

    // Set to preferred chain to default so that we try all three
    // alignments (since we encountered an alignment error, best to
    // be safe).
    preferred_chain = "";
    err_code = this->align_input_subregion(best_result,
            subsequence, preferred_chain);

    for (size_t i = 0; i < (sequence.length() - subsequence.length());
            i++) {
        std::get<0>(best_result).push_back("-");
    }

    return best_result;
}




/// Aligns the input sequence, which may be either the full sequence or a subregion of it.
int SingleChainAnnotatorCpp::align_input_subregion(std::tuple<std::vector<std::string>, double,
                std::string, std::string> &best_result, std::string &query_sequence,
                std::string preferred_chain) {
    auto queryAsIdx = std::make_unique<int[]>(query_sequence.length());
    bool fwgxg_error = false;

    // Set initial percent identity to minimum.
    std::get<1>(best_result) = 0.0;

    if (!SequenceUtilities::convert_x_sequence_to_array(queryAsIdx.get(),
                query_sequence))
        return INVALID_SEQUENCE;

    // If preferred chain is specified (i.e. is not ""), we can
    // only align to one scoring tool...neat! In that case, score
    // using that tool then return, AS LONG AS the preferred chain
    // matches one of the chains specified for this object.
    if (preferred_chain != "") {
        auto first_chain = this->chains.begin();
        auto last_chain = this->chains.end();
        auto selection = std::find(first_chain, last_chain,
                preferred_chain);

        if (selection != this->chains.end()) {
            size_t aligner_id = selection - this->chains.begin();
            std::vector<std::string> final_numbering;
            std::string error_message = "";
            double percent_identity = -1;

            this->scoring_tools[aligner_id].align(
                    query_sequence, queryAsIdx.get(), final_numbering,
                    percent_identity, error_message);

            std::get<0>(best_result) = final_numbering;
            std::get<1>(best_result) = percent_identity;
            std::get<2>(best_result) = this->scoring_tools[aligner_id].
                get_chain_name();
            std::get<3>(best_result) = error_message;
            return VALID_SEQUENCE;
        }
    }



    // Otherwise, if no preferred chain was specified or if the preferred
    // chain is not one this object can handle (e.g. it was initialized
    // with chain H only but preferred is K), we need to
    // try all available aligners. Save the one that yields the highest
    // percent identity.

    for (size_t i=0; i < this->scoring_tools.size(); i++) {
        std::vector<std::string> final_numbering;
        std::string error_message = "";
        double percent_identity = -1;

        final_numbering.reserve(query_sequence.length() * 2);

        this->scoring_tools[i].align(query_sequence,
                    queryAsIdx.get(), final_numbering, percent_identity,
                    error_message);
        if (percent_identity > std::get<1>(best_result)) {
            std::get<0>(best_result) = final_numbering;
            std::get<1>(best_result) = percent_identity;
            std::get<2>(best_result) = this->scoring_tools[i].get_chain_name();
            std::get<3>(best_result) = error_message;
        }
        if (std::get<3>(best_result).length() > 2) {
            if (std::get<3>(best_result).substr(0, 2) == "> ")
                fwgxg_error = true;
        }
    }
    if (fwgxg_error && std::get<1>(best_result) < 0.75)
        return POSSIBLE_FWGXG_ERROR_ON_ALIGNMENT;
    return VALID_SEQUENCE;
}






/// SingleChainAnnotatorCpp function which numbers a list of input sequences.
/// Essentially a wrapper on analyze_seq that is convenient if user wants
/// to pass a list.
std::vector<std::tuple<std::vector<std::string>, double, std::string,
            std::string>> SingleChainAnnotatorCpp::analyze_seqs(std::vector<std::string> sequences) {
    std::vector<std::tuple<std::vector<std::string>, double, std::string,
            std::string>> output_results;

    for (size_t i=0; i < sequences.size(); i++)
        output_results.push_back(this->analyze_seq(sequences[i]));

    return output_results;
}

}  // namespace NumberingTools
