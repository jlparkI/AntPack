/* Contains the implementation of the paired chain annotator tool for
 * sequences containing a paired heavy and light chain. */

// C++ headers
#include <memory>

// Project headers
#include "paired_chain_annotator.h"



namespace NumberingTools{

PairedChainAnnotatorCpp::PairedChainAnnotatorCpp(
        std::string scheme, std::string consensus_filepath,
        std::unordered_map<std::string, size_t> nterm_kmers
        ):
    AnnotatorBaseClassCpp(scheme, consensus_filepath, nterm_kmers),
    scheme(scheme) {

    std::vector<std::string> chains = {"K", "L"};
    this->light_chain_analyzer = std::make_unique<SingleChainAnnotatorCpp>
        (chains, scheme, consensus_filepath, nterm_kmers);

    chains = {"H"};
    this->heavy_chain_analyzer = std::make_unique<SingleChainAnnotatorCpp>
        (chains, scheme, consensus_filepath, nterm_kmers);

    chains = {"H", "K", "L"};
    this->analyzer = std::make_unique<SingleChainAnnotatorCpp>
        (chains, scheme, consensus_filepath, nterm_kmers);


    std::filesystem::path extensionPath = consensus_filepath;

    std::string npyFName = "CTERMFINDER_CONSENSUS_H.npy";
    std::filesystem::path npyFPath = extensionPath / npyFName;
    cnpy::NpyArray raw_score_arr;
    try {
        raw_score_arr = cnpy::npy_load(npyFPath.string());
    }
    catch (...) {
        throw std::runtime_error(std::string("The consensus file / "
                    "library installation has an issue."));
    }
}



// PairedChainAnnotatorCpp function which numbers an input sequence.
// An alignment is returned for the heavy chain and light chain that
// are found.
std::pair<std::tuple<std::vector<std::string>, double, std::string, std::string>,
    std::tuple<std::vector<std::string>, double, std::string, std::string>>
    PairedChainAnnotatorCpp::analyze_seq(std::string sequence) {
    
    if (!SequenceUtilities::validate_x_sequence(sequence)) {
        std::tuple<std::vector<std::string>, double,
            std::string, std::string> vecres =
            {{}, 0., "", "Sequence contains invalid characters"};
        std::pair<std::tuple<std::vector<std::string>, double,
            std::string, std::string>,     std::tuple<std::vector<std::string>,
            double, std::string, std::string>> output_result = {vecres, vecres};
        return output_result;
    }

    // There are three possible chains.
    std::array<double, 3> cterm_scores;
    std::array<int, 3> cterm_positions, nterm_positions, nterm_scores;

    // Find likely c- and n-terminal regions in the input sequence.
    if (!this->boundary_finder->find_start_end_zones(sequence, cterm_scores,
            cterm_positions, nterm_scores, nterm_positions)) {
        std::tuple<std::vector<std::string>, double,
            std::string, std::string> vecres =
            {{}, 0., "", "Invalid sequence supplied -- either too short or "
                "contains nonstandard AAs"};
        std::pair<std::tuple<std::vector<std::string>, double,
            std::string, std::string>,     std::tuple<std::vector<std::string>,
            double, std::string, std::string>> output_result = {vecres, vecres};
        return output_result;
    }

    // Cterminal and nterminal scores are not always sufficiently reliable
    // to distinguish between K and L. We therefore merge the light chain
    // scores and if K or L is indicated use the light chain aligner.
    std::array<int, 2> merged_nterm_scores, merged_nterm_positions,
        merged_cterm_positions;
    std::array<double, 2> merged_cterm_scores;
    this->boundary_finder->merge_light_chain_nterm_scores(nterm_scores,
            nterm_positions, merged_nterm_scores, merged_nterm_positions);
    this->boundary_finder->merge_light_chain_cterm_scores(cterm_scores,
            cterm_positions, merged_cterm_scores, merged_cterm_positions);

    // Ask the c-terminal finder for a threshold below which a score likely does
    // not indicate the presence of a c-terminal and store this.
    double best_cterm_score = this->boundary_finder->
        get_cterm_threshold_score();
    int best_nterm_score = this->boundary_finder->
        get_nterm_threshold_score();

    int start_cutoff = 0, end_cutoff = sequence.length();

    // Since this is a paired chain, there should be at least one and probably
    // 2 c-terminal regions. Start by finding the highest-confidence c-terminal
    // region. If we do not find anything, we keep the length of the sequence
    // as the end cutoff.
    std::vector<std::string> chain_list = this->boundary_finder->
            get_chain_list();
    int selected_idx = -1;

    for (size_t i=0; i < merged_cterm_scores.size(); i++) {
        if (merged_cterm_scores[i] > best_cterm_score) {
            best_cterm_score = merged_cterm_scores[i];
            selected_idx = i;
            // Add onto the end to make sure all of cterminal is
            // included plus a safety margin.
            end_cutoff = merged_cterm_positions[i] + this->boundary_finder->
                get_num_positions() + 5;
            end_cutoff = end_cutoff < sequence.length() ?
                end_cutoff : sequence.length();
        }
    }


    // Try to find something looking like an nterminal. If we
    // found a specific cterminal, we should look for the
    // matching nterminal. If we did not, we should look for
    // an n-terminal.
    if (selected_idx != -1) {
        if (merged_nterm_scores[selected_idx] > best_nterm_score) {
            start_cutoff = merged_nterm_positions[selected_idx] -
                this->boundary_finder->get_kmer_window_size() - 5;
        }
    } else {
        for (size_t i=0; i < merged_nterm_scores.size(); i++) {
            if (merged_nterm_scores[i] > best_nterm_score) {
                best_nterm_score = merged_nterm_scores[i];
                // Subtract from the end to make sure all of nterminal is
                // included plus a safety margin.
                start_cutoff = merged_nterm_positions[i] -
                    this->boundary_finder->get_kmer_window_size() - 5;
                selected_idx = i;
            }
        }
    }

    // If the start cutoff is too close to the end cutoff, reset it to
    // a reasonable margin.
    if (end_cutoff - start_cutoff < 100) {
        start_cutoff = end_cutoff - 200;
        if (start_cutoff < 0) {
            start_cutoff = 0;
            if (end_cutoff - start_cutoff < 100)
                end_cutoff = start_cutoff + 200;
        }
    }

    // Of course, always double check that everything is within bounds.
    start_cutoff = start_cutoff > 0 ? start_cutoff : 0;
    end_cutoff = end_cutoff < sequence.length() ? end_cutoff :
        sequence.length();
    std::string subsequence = sequence.substr(start_cutoff,
            end_cutoff - start_cutoff);

    std::tuple<std::vector<std::string>, double,
            std::string, std::string> init_results =
                this->analyzer->analyze_seq(subsequence);


    // If the initial alignment returned an error, the numbering
    // vector will be empty and chain assignment may be arbitrary.
    // This may happen if user passed a single chain by mistake.
    // We can try realigning with no cutoffs and no preferred chain
    // and see what we get.
    int err_code = VALID_SEQUENCE;
    if (std::get<0>(init_results).size() == 0) {
        start_cutoff = 0;
        end_cutoff = sequence.length();
        err_code = this->analyzer->align_input_subregion(init_results,
                sequence, "");
    }

    // If we still have the same problem, abort.
    if (std::get<0>(init_results).size() == 0 || err_code != VALID_SEQUENCE) {
        std::tuple<std::vector<std::string>, double,
            std::string, std::string> blank =
            {{}, 0., "", "Alignment error; cterminal, nterminal not found."};
        std::pair<std::tuple<std::vector<std::string>, double,
            std::string, std::string>,     std::tuple<std::vector<std::string>,
            double, std::string, std::string>> output_result;

        if (std::get<2>(init_results) == "H")
            output_result = {init_results, blank};
        else
            output_result = {blank, init_results};
        return output_result;
    }

    // Otherwise, we can proceed. Add gaps onto the initial
    // alignment so it is the same length as the sequence.
    std::get<0>(init_results) = this->pad_alignment(sequence,
            std::get<0>(init_results), start_cutoff,
            end_cutoff);

    // Now we need to align either the region before init_results or
    // the region after. Find where the alignment starts and
    // ends.

    size_t postalign_start = 0, postalign_end =
        std::get<0>(init_results).size();

    // Now number the region that is gaps in the first alignment.
    for (size_t i=0; i < std::get<0>(init_results).size(); i++) {
        if (std::get<0>(init_results)[i] != "-") {
            postalign_start = i;
            break;
        }
    }
    for (int i=std::get<0>(init_results).size();
            i > 0; i--) {
        if (std::get<0>(init_results)[i-1] != "-") {
            postalign_end = i;
            break;
        }
    }

    // If exstart or exend is too close to the ends of the
    // sequence to do an alignment, skip it. If both have
    // to be skipped, then user passed a single chain by
    // mistake, abort and return results. For now we are
    // aligning both the start and end region -- TODO:
    // for greater efficiency identify the most promising
    // region using the nterm and cterm scores.
    std::tuple<std::vector<std::string>, double, std::string,
        std::string> second_result;
    std::tuple<std::vector<std::string>, double, std::string,
        std::string> third_result;

    // Initialize to pessimistic values.
    std::get<1>(second_result) = -1;
    std::get<1>(third_result) = -1;
    std::get<3>(second_result) = "Alignment error; chain not found.";
    std::get<3>(third_result) = "Alignment error; chain not found.";
    std::get<2>(second_result) = "";
    std::get<2>(third_result) = "";

    if (sequence.length() - postalign_end > MINIMUM_SEQUENCE_LENGTH) {
        subsequence = sequence.substr(postalign_end, sequence.length() -
                postalign_end);
        if (std::get<2>(init_results) == "H") {
            err_code = this->light_chain_analyzer->align_input_subregion(
                    second_result, subsequence, "");
        } else {
            err_code = this->heavy_chain_analyzer->align_input_subregion(
                    second_result, subsequence, "");
        }
    }
    if (postalign_start > MINIMUM_SEQUENCE_LENGTH) {
        subsequence = sequence.substr(0, postalign_start);
        if (std::get<2>(init_results) == "H") {
            err_code = this->light_chain_analyzer->align_input_subregion(
                    third_result, subsequence, "");
        } else {
            err_code = this->heavy_chain_analyzer->align_input_subregion(
                    third_result, subsequence, "");
        }
    }


    std::pair<std::tuple<std::vector<std::string>, double,
        std::string, std::string>,     std::tuple<std::vector<std::string>,
        double, std::string, std::string>> output_result;

    if (std::get<1>(second_result) > std::get<1>(third_result)) {
        std::get<0>(second_result) = this->pad_alignment(sequence,
            std::get<0>(second_result), postalign_end,
            sequence.length());

        if (std::get<2>(init_results) == "H")
            output_result = {init_results, second_result};
        else
            output_result = {second_result, init_results};
    } else {
        std::get<0>(third_result) = this->pad_alignment(sequence,
            std::get<0>(third_result), 0, postalign_start);

        if (std::get<2>(init_results) == "H")
            output_result = {init_results, third_result};
        else
            output_result = {third_result, init_results};
    }

    return output_result;
}



/// PairedChainAnnotatorCpp function which numbers a list of input sequences.
/// Basically a wrapper on analyze_seq for convenience if user passes a list.
std::tuple<std::vector<std::tuple<std::vector<std::string>, double, std::string, std::string>>,
           std::vector<std::tuple<std::vector<std::string>, double, std::string, std::string>>>
             PairedChainAnnotatorCpp::analyze_seqs(std::vector<std::string> sequences) {
    std::tuple<std::vector<std::tuple<std::vector<std::string>, double, std::string, std::string>>,
           std::vector<std::tuple<std::vector<std::string>, double, std::string, std::string>>>
               output_results;

    std::pair<std::tuple<std::vector<std::string>, double,
        std::string, std::string>,
            std::tuple<std::vector<std::string>, double,
            std::string, std::string>> seq_result;

    for (size_t i=0; i < sequences.size(); i++) {
        seq_result = this->analyze_seq(sequences[i]);
        std::get<0>(output_results).push_back(std::move(seq_result.first));
        std::get<1>(output_results).push_back(std::move(seq_result.second));
    }
    return output_results;
}



/// Pads the input alignment so it is the same length as the
/// query sequence.
std::vector<std::string> PairedChainAnnotatorCpp::pad_alignment(
        const std::string &query_sequence,
                const std::vector<std::string> &alignment,
                const int &align_start, const int &align_end) {
    std::vector<std::string> padded_numbering;

    for (size_t i=0; i < align_start; i++)
        padded_numbering.push_back("-");
    for (size_t i=0; i < alignment.size(); i++)
        padded_numbering.push_back(alignment.at(i));
    for (size_t i=align_end; i < query_sequence.length(); i++)
        padded_numbering.push_back("-");

    return padded_numbering;
}



}  // namespace NumberingTools
