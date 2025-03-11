/* Copyright (C) 2025 Jonathan Parkinson
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
// C++ headers
#include <utility>
#include <tuple>
#include <vector>
#include <memory>
#include <string>
#include <algorithm>


// Project headers
#include "single_chain_annotator.h"
#include "prefiltering_tool.h"





namespace SequenceAnnotators {


// Special error code for cases where a rare but specific type
// of alignment error may have occurred.
static const int POSSIBLE_FWGXG_ERROR_ON_ALIGNMENT = 2;


SingleChainAnnotatorCpp::SingleChainAnnotatorCpp(
    std::vector<std::string> chains,
    std::string scheme, std::string consensus_filepath,
    std::unordered_map<std::string, size_t> nterm_kmers):
AnnotatorBaseClassCpp(scheme),
chains(chains),
scheme(scheme) {
    // First, determine whether the chains that the user passed are TCR or
    // mAb. If the user supplied a combination, generate an exception
    // (they should never do this).
    int tcr_chain_count = 0, mab_chain_count = 0;
    for (const auto & chain : chains) {
        if (chain == "H" || chain == "K" || chain == "L") {
            mab_chain_count += 1;
        } else if (chain == "A" || chain == "B" || chain == "D" ||
                chain == "G") {
            tcr_chain_count += 1;
        } else {
            throw std::runtime_error(std::string("Unrecognized "
                        "chain supplied."));
        }
    }

    if (tcr_chain_count > 0 && mab_chain_count > 0) {
        throw std::runtime_error(std::string("You have supplied chains "
                    "for TCRs and chains for mAbs to the same annotator. "
                    "You cannot have an annotator for both mAbs and TCRs "
                    "-- if you need to analyze both types of sequences, "
                    "you should create an annotator for each. When "
                    "creating an annotator, supply either TCR chains "
                    "(A, B, D, G) or mAb chains (H, K, L) -- not both."));
    }

    if (mab_chain_count > 0) {
        this->boundary_finder = std::make_unique
            <PrefilteringRoutines::PrefilteringTool>(consensus_filepath,
                nterm_kmers);

        // Note that exceptions thrown here go back to Python via
        // PyBind as long as this constructor is used within the wrapper.
        if (chains.size() < 1 || chains.size() > 3) {
            throw std::runtime_error(std::string("There must be at least "
                        "one chain specified, and no more than 3."));
        }

        if (scheme != "aho" && scheme != "imgt" && scheme != "kabat"
                && scheme != "martin") {
            throw std::runtime_error(std::string("Invalid scheme "
                            "specified. Please use one of 'imgt', 'kabat', "
                            "'martin', 'aho'."));
        }

        for (auto & chain : this->chains) {
            if (chain != "H" && chain != "K" && chain != "L") {
                throw std::runtime_error(std::string("For mAbs, valid "
                            "chains must be one of 'H', 'K', 'L'."));
            } else {
                this->scoring_tools.push_back(
                        NumberingTools::IGAligner(consensus_filepath,
                        chain, scheme) );
            }
        }
    }

    if (tcr_chain_count > 0) {
        if (chains.size() < 1 || chains.size() > 4) {
            throw std::runtime_error(std::string("There must be at least "
                        "one chain specified, and no more than 4."));
        }
        if (scheme != "imgt") {
            throw std::runtime_error(std::string("Invalid scheme "
                            "specified. For TCRs, please use IMGT."));
        }
        for (auto & chain : this->chains) {
            if (chain != "A" && chain != "B" && chain != "D" &&
                    chain != "G") {
                throw std::runtime_error(std::string("For TCRs, "
                            "valid chains must be one of 'A', 'B', "
                            "'D', 'G'."));
            } else {
                this->tcr_scoring_tools.push_back(
                        NumberingTools::VJAligner(consensus_filepath,
                        chain, scheme, "TCR") );
            }
        }
    }
}

/// @brief Numbers a single input sequence by internally
/// routing to mab_analyze_seq or tcr_analyze_seq depending
/// on what the class was initialized to do.
/// @param sequence The input sequence.
/// @return A tuple of numbering, percent_identity, chain_type
/// and error_message.
std::tuple<std::vector<std::string>, double, std::string,
    std::string> SingleChainAnnotatorCpp::analyze_seq(
    std::string sequence) {
    if (this->scoring_tools.size() > 0) {
        return this->mab_analyze_seq(sequence);
    } else if (this->tcr_scoring_tools.size() > 0) {
        return this->tcr_analyze_seq(sequence);
    } else {
        // This should never happen -- class is always initialized with
        // one or the other and if not throws an exception during
        // construction.
        std::vector<std::string> empty_numbering;
        std::tuple<std::vector<std::string>, double, std::string,
                    std::string> best_result{ empty_numbering,
                        0, "", "Annotator not initialized properly."};
        return best_result;
    }
}


/// @brief Numbers a single input sequence when the object has been
/// initialized to annotate antibody sequences specifically.
/// @param sequence The input sequence.
/// @return A tuple of numbering, percent_identity, chain_type
/// and error_message.
std::tuple<std::vector<std::string>, double, std::string,
        std::string> SingleChainAnnotatorCpp::mab_analyze_seq(
                std::string sequence) {
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

    int err_code = this->mab_align_input_subregion(best_result, sequence,
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
    best_seq_cutoff = best_seq_cutoff >
        NumberingTools::MINIMUM_SEQUENCE_LENGTH ?
        best_seq_cutoff : NumberingTools::MINIMUM_SEQUENCE_LENGTH;
    best_seq_cutoff = best_seq_cutoff < sequence.length() ? best_seq_cutoff :
        sequence.length();

    std::string subsequence = sequence.substr(0, best_seq_cutoff);

    // Set to preferred chain to default so that we try all three
    // alignments (since we encountered an alignment error, best to
    // be safe).
    preferred_chain = "";
    err_code = this->mab_align_input_subregion(best_result,
            subsequence, preferred_chain);

    for (size_t i = 0; i < (sequence.length() - subsequence.length());
            i++) {
        std::get<0>(best_result).push_back("-");
    }

    return best_result;
}



/// @brief Numbers a single input sequence when the object has been
/// initialized to annotate TCR sequences specifically.
/// @param sequence The input sequence.
/// @return A tuple of numbering, percent_identity, chain_type
/// and error_message.
std::tuple<std::vector<std::string>, double, std::string,
        std::string> SingleChainAnnotatorCpp::tcr_analyze_seq(
                std::string sequence) {
    std::string preferred_chain = "";
    int preferred_vgene = 0, preferred_jgene = 0;
    std::vector<std::string> empty_numbering;
    std::tuple<std::vector<std::string>, double, std::string,
                    std::string> best_result{ empty_numbering,
                        0, "", "Sequence contains invalid characters"};


    if (!SequenceUtilities::validate_x_sequence(sequence))
        return best_result;

    auto queryAsIdx = std::make_unique<int[]>(sequence.length());
    if (!SequenceUtilities::convert_x_sequence_to_array(queryAsIdx.get(),
                sequence))
        return best_result;

    // First, determine using the available VJAligners which chain
    // is the best match. If both V and J genes concur, this is
    // straightforward, and we can set a preferred_chain.
    int best_jg_position = 0, best_vg_score = -1, best_jg_score = -1;
    int best_vg_idx = -1, best_jg_idx = -1;
    std::string best_vg_chain, best_jg_chain;
    std::vector<int> best_vgenes, best_jgenes;

    for (auto & aligner : this->tcr_scoring_tools) {
        int vg_score, jg_score, vg_idx, jg_idx, jg_position;
        if (!aligner.identify_best_vgene(sequence,
                queryAsIdx.get(), vg_score, vg_idx)) {
            return best_result;
        }
        if (!aligner.identify_best_jgene(sequence,
                    queryAsIdx.get(), jg_position, jg_score, jg_idx)) {
            return best_result;
        }

        best_vgenes.push_back(vg_idx);
        best_jgenes.push_back(jg_idx);

        if (vg_score > best_vg_score) {
            best_vg_score = vg_score;
            best_vg_idx = vg_idx;
            best_vg_chain = aligner.get_chain_name();
        }
        if (jg_score > best_jg_score) {
            best_jg_position = jg_position;
            best_jg_score = jg_score;
            best_jg_idx = jg_idx;
            best_jg_chain = aligner.get_chain_name();
        }
    }

    // If we identified different chains using jgene and vgene,
    // or if no chain was identified, try aligning to all possible
    // chains.
    if (best_jg_chain != best_vg_chain || best_vg_idx < 0 ||
            best_jg_idx < 0) {
        for (size_t i=0; i < this->tcr_scoring_tools.size(); i++) {
            std::tuple<std::vector<std::string>, double, std::string,
                    std::string> current_result;

            if (!this->tcr_align_input_subregion(current_result,
                    sequence, queryAsIdx.get(), this->chains[i],
                    best_vgenes.at(i), best_jgenes.at(i))) {
                std::tuple<std::vector<std::string>, double, std::string,
                    std::string> null_result{ empty_numbering,
                        0, "", "Error in sequence processing"};
                return null_result;
            }
            if (std::get<1>(current_result) > std::get<1>(best_result))
                best_result = current_result;
        }
        return best_result;
    }

    // Otherwise, we have agreement on the v and jgene, which
    // strongly suggests a particular chain. Align only to that
    // chain.
    if (!this->tcr_align_input_subregion(best_result,
            sequence, queryAsIdx.get(), best_vg_chain,
            best_vg_idx, best_jg_idx)) {
        std::tuple<std::vector<std::string>, double, std::string,
            std::string> null_result{ empty_numbering,
                0, "", "Error in sequence processing"};
        return null_result;
    }

    return best_result;
}



/// @brief Aligns a subregion of an input sequence (for situations where
/// a subregion needs to be extracted) in cases where the object has
/// been initialized to analyze mab sequences.
/// @param best_result The tuple of numbering, percent_identity,
/// chain_type and error_message in which output will be stored.
/// @param query_sequence The input sequence.
/// @param preferred_chain Indicates if a specific chain has already
/// been found (based on kmer counting) to be a better match; if
/// so, only this one needs to be aligned.
int SingleChainAnnotatorCpp::mab_align_input_subregion(
    std::tuple<std::vector<std::string>, double,
    std::string, std::string> &best_result, std::string &query_sequence,
    std::string preferred_chain) {
    bool fwgxg_error = false;

    // Set initial percent identity to minimum.
    std::get<1>(best_result) = 0.0;

    auto queryAsIdx = std::make_unique<int[]>(query_sequence.length());
    if (!SequenceUtilities::convert_x_sequence_to_array(queryAsIdx.get(),
                query_sequence))
        return NumberingTools::INVALID_SEQUENCE;

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
            return NumberingTools::VALID_SEQUENCE;
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
    return NumberingTools::VALID_SEQUENCE;
}



/// @brief Aligns a subregion of an input sequence (for situations where
/// a subregion needs to be extracted) in cases where the object has
/// been initialized to analyze TCR sequences.
/// @param best_result The tuple of numbering, percent_identity,
/// chain_type and error_message in which output will be stored.
/// @param query_sequence The input sequence.
/// @param preferred_chain Indicates which chain should be used. Note
/// that in contrast to mab_align_input_subregion, for TCRs the
/// preferred chain must be specified.
/// @param preferred_vgene Indicates which vgene to use (notice that
/// this is different from mab_align_input_subregion).
/// @param preferred_jgene Indicates which jgene to use (notice that
/// this is different from mab_align_input_subregion).
int SingleChainAnnotatorCpp::tcr_align_input_subregion(
    std::tuple<std::vector<std::string>, double,
        std::string, std::string> &best_result,
    std::string &query_sequence,
    int *queryAsIdx,
    const std::string &preferred_chain,
    const int &preferred_vgene,
    const int &preferred_jgene) {

    // Set initial percent identity to minimum.
    std::get<1>(best_result) = 0.0;

    // Notice that in contrast to mab_align_input_subregion,
    // a preferred chain MUST be specified (since otherwise this
    // function would have to again find the preferred vgene and
    // jgene which have already been found by caller).
    auto first_chain = this->chains.begin();
    auto last_chain = this->chains.end();
    auto selection = std::find(first_chain, last_chain,
            preferred_chain);

    if (selection != this->chains.end()) {
        size_t aligner_id = selection - this->chains.begin();
        std::vector<std::string> final_numbering;
        std::string error_message = "";
        double percent_identity = -1;

        this->tcr_scoring_tools[aligner_id].align(
                query_sequence, queryAsIdx, final_numbering,
                percent_identity, error_message,
                preferred_vgene, preferred_jgene);

        std::get<0>(best_result) = final_numbering;
        std::get<1>(best_result) = percent_identity;
        std::get<2>(best_result) = this->tcr_scoring_tools[aligner_id].
            get_chain_name();
        std::get<3>(best_result) = error_message;
        return NumberingTools::VALID_SEQUENCE;
    }
    // If we could not find the specified chain in the list of
    // aligners, something went wrong (highly unlikely, because
    // that would mean the object was finding in its list of
    // VJAligners an object that did not exist in its list of
    // VJAligners, but add this error return here just in case).
    return NumberingTools::INVALID_SEQUENCE;
}




/// @brief Numbers a list of input sequences
/// @param query_sequences A vector of sequences to number.
/// @return A vector of tuples of the same length as the input vector
/// of sequences. Each tuple contains the numbering, percent identity,
/// chain type and error message for the corresponding sequence.
std::vector<std::tuple<std::vector<std::string>, double, std::string,
    std::string>> SingleChainAnnotatorCpp::analyze_seqs(std::vector<std::string> sequences) {
    std::vector<std::tuple<std::vector<std::string>, double, std::string,
            std::string>> output_results;

    if (this->scoring_tools.size() > 0) {
        for (size_t i=0; i < sequences.size(); i++)
            output_results.push_back(this->mab_analyze_seq(sequences[i]));
    } else if (this->tcr_scoring_tools.size() > 0) {
        for (size_t i=0; i < sequences.size(); i++)
            output_results.push_back(this->tcr_analyze_seq(sequences[i]));
    } else {
        // This should never happen -- class is always initialized with
        // one or the other and if not throws an exception during
        // construction.
        for (size_t i=0; i < sequences.size(); i++)
            output_results.push_back(this->analyze_seq(sequences[i]));
    }

    return output_results;
}

}  // namespace SequenceAnnotators
