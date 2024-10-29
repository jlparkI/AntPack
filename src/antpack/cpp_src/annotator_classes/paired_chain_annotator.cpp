/* Contains the implementation of the paired chain annotator tool for
 * sequences containing a paired heavy and light chain. */

// C++ headers
#include <memory>

// Project headers
#include "paired_chain_annotator.h"



namespace NumberingTools{

PairedChainAnnotatorCpp::PairedChainAnnotatorCpp(
        std::string scheme, std::string consensus_filepath):
    AnnotatorBaseClassCpp(scheme, consensus_filepath),
    scheme(scheme) {

    std::vector<std::string> chains = {"K", "L"};
    this->light_chain_analyzer = std::make_unique<SingleChainAnnotatorCpp>
        (chains, scheme, false, consensus_filepath);

    chains = {"H"};
    this->heavy_chain_analyzer = std::make_unique<SingleChainAnnotatorCpp>
        (chains, scheme, false, consensus_filepath);

    chains = {"H", "K", "L"};
    this->analyzer = std::make_unique<SingleChainAnnotatorCpp>
        (chains, scheme, false, consensus_filepath);


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
            {{}, 0., "", "Invalid sequence supplied -- nonstandard AAs"};
        std::pair<std::tuple<std::vector<std::string>, double,
            std::string, std::string>,     std::tuple<std::vector<std::string>,
            double, std::string, std::string>> output_result = {vecres, vecres};
        return output_result;
    }

    // There are three possible chains, so we will initially generate
    // three scores, then merge the two light chain scores.
    std::array<double, 3> scores;
    std::array<int, 3> positions;
    std::array<double, 2> merged_scores;
    std::array<int, 2> merged_positions;

    // Find likely c-terminal regions in the input sequence.
    if (!this->boundary_finder->find_best_cterminal(sequence, scores,
            positions)) {
        std::tuple<std::vector<std::string>, double,
            std::string, std::string> vecres =
            {{}, 0., "", "Alignment error -- no c-terminal region found."};
        std::pair<std::tuple<std::vector<std::string>, double,
            std::string, std::string>,     std::tuple<std::vector<std::string>,
            double, std::string, std::string>> output_result = {vecres, vecres};
        return output_result;
    }

    this->boundary_finder->merge_light_chain_scores(scores, positions,
            merged_scores, merged_positions);

    // Ask the c-terminal finder for a threshold below which a score likely does
    // not indicate the presence of a c-terminal and store this.
    double score_threshold = this->boundary_finder->get_threshold_score();

    int cutoff_position = 0;

    // Since this is a paired chain, there should be at least one and probably
    // 2 c-terminal regions. Which one occurs first? We ensure the score is
    // > cutoff to rule out regions that are likely not "real". If no region
    // meets the cutoff, input sequence is either not an antibody or a severely
    // truncated one. In that case, try to analyze it with the
    // SingleChainAnnotator and return that result, whatever it is.
    if (merged_scores[0] <= score_threshold &&
            merged_scores[1] <= score_threshold) {
        std::tuple<std::vector<std::string>, double,
            std::string, std::string> init_results =
            this->analyzer->analyze_seq(sequence);
        std::pair<std::tuple<std::vector<std::string>, double,
            std::string, std::string>,     std::tuple<std::vector<std::string>,
            double, std::string, std::string>> output_result;
        std::tuple<std::vector<std::string>, double,
            std::string, std::string> second_result =
            {{}, 0., "", "Alignment error; second chain not numbered "
                "correctly / not found."};

        if (std::get<2>(init_results) == "H")
            output_result = {init_results, second_result};
        else
            output_result = {second_result, init_results};

        return output_result;
    }

    // Otherwise, if at least one region meets the cutoff, try
    // to find whether light or heavy chain occurs first.
    if (merged_positions[0] < merged_positions[1] &&
            merged_scores[0] > score_threshold) {
        cutoff_position = merged_positions[0];
    } else if (merged_scores[1] > score_threshold) {
        cutoff_position = merged_positions[1];
    } else {
        cutoff_position = merged_positions[0];
    }

    // Make sure cutoff position obeys reasonable constraints.
    cutoff_position = cutoff_position > 0 ? cutoff_position : 0;
    cutoff_position = cutoff_position < sequence.length() ?
        cutoff_position : sequence.length();


    std::string subsequence = sequence.substr(cutoff_position,
            sequence.length() - cutoff_position);

    std::tuple<std::vector<std::string>, double,
        std::string, std::string> init_results =
        this->analyzer->analyze_seq(subsequence);

    // If the initial alignment returned an error, the numbering
    // vector will be empty and chain assignment may be arbitrary.
    // This may happen if user passed a single chain by mistake.
    // If so, abort and use a simpler procedure for the remaining
    // analysis.
    if (std::get<0>(init_results).size() == 0) {
        cutoff_position += this->boundary_finder->get_num_positions() + 5;
        cutoff_position = cutoff_position < sequence.length() ?
            cutoff_position : sequence.length();
        subsequence = sequence.substr(0, cutoff_position);

        // Change the error message in init_result to something
        // more accurate.
        std::get<3>(init_results) = "Alignment error; second chain "
            "not numbered correctly / not found.";

        std::tuple<std::vector<std::string>, double, std::string,
                std::string> second_result;
        second_result = this->analyzer->analyze_seq(subsequence);

        size_t padding_amount = (sequence.length() -
                std::get<0>(second_result).size());

        for (size_t i=0; i < padding_amount; i++)
            std::get<0>(second_result).push_back("-");

        std::pair<std::tuple<std::vector<std::string>, double,
            std::string, std::string>,     std::tuple<std::vector<std::string>,
            double, std::string, std::string>> output_result;

        if (std::get<2>(second_result) == "H") {
            std::get<2>(init_results) = "L";
            output_result = {second_result, init_results};
        } else {
            std::get<2>(init_results) = "H";
            output_result = {init_results, second_result};
        }

        return output_result;
    }

    // Otherwise, we can proceed. Add gaps onto the initial
    // alignment so it is the same length as the sequence.
    std::vector<std::string> padded_init_numbering((sequence.length() -
                std::get<0>(init_results).size()), "-");

    for (size_t i=0; i < std::get<0>(init_results).size(); i++)
        padded_init_numbering.push_back(std::get<0>(init_results).at(i));

    std::get<0>(init_results) = std::move(padded_init_numbering);

    size_t exstart = 0;

    // Now number the region that is gaps in the first alignment.
    for (size_t i=0; i < std::get<0>(init_results).size(); i++) {
        if (std::get<0>(init_results)[i] != "-") {
            exstart = i;
            break;
        }
    }

    // It is possible that exstart is 0 (e.g. if user passed a
    // single chain sequence by mistake). In this case the aligner
    // will immediately reject with an appropriate error message
    // so this is ok.

    subsequence = sequence.substr(0, exstart);
    std::tuple<std::vector<std::string>, double, std::string,
        std::string> second_result;


    if (std::get<2>(init_results) == "H")
        second_result = this->light_chain_analyzer->analyze_seq(subsequence);
    else
        second_result = this->heavy_chain_analyzer->analyze_seq(subsequence);

    size_t padding_amount = (sequence.length() -
            std::get<0>(second_result).size());

    for (size_t i=0; i < padding_amount; i++)
        std::get<0>(second_result).push_back("-");

    std::pair<std::tuple<std::vector<std::string>, double,
        std::string, std::string>,     std::tuple<std::vector<std::string>,
        double, std::string, std::string>> output_result;

    if (std::get<2>(init_results) == "H")
        output_result = {init_results, second_result};

    else
        output_result = {second_result, init_results};

    return output_result;
}


// PairedChainAnnotatorCpp function which numbers a list of input sequences.
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

}  // namespace NumberingTools
