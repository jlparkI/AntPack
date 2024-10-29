#include "cterm_finder.h"

// C++ headers
#include <vector>
#include <memory>
#include <string>


namespace NumberingTools {


CTermFinder::CTermFinder(std::string consensus_filepath) {
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

    this->score_arr_shape[0] = raw_score_arr.shape[0];
    this->score_arr_shape[1] = raw_score_arr.shape[1];
    this->score_arr_shape[2] = boundary_chains.size();
    if (this->score_arr_shape[1] != 21) {
        throw std::runtime_error(std::string("The consensus file / "
                    "library installation has an issue."));
    }
    this->score_array = std::make_unique<double[]>(this->score_arr_shape[0] *
           this->score_arr_shape[1] * this->score_arr_shape[2]);

    for (size_t i=0; i < boundary_chains.size(); i++) {
        std::string npyFName = "CTERMFINDER_CONSENSUS_" +
            boundary_chains[i] + ".npy";
        std::filesystem::path npyFPath = extensionPath / npyFName;

        try {
            raw_score_arr = cnpy::npy_load(npyFPath.string());
        }
        catch (...) {
            throw std::runtime_error(std::string("The consensus file / "
                        "library installation has an issue."));
        }

        int ld_shape0 = raw_score_arr.shape[0];
        int ld_shape1 = raw_score_arr.shape[1];

        if (raw_score_arr.word_size != 8 ||
                ld_shape0 != this->score_arr_shape[0] ||
                ld_shape1 != this->score_arr_shape[1]) {
            throw std::runtime_error("The consensus file / "
                        "library installation has an issue.");
        }

        double *raw_score_ptr = raw_score_arr.data<double>();

        for (int j=0; j < this->score_arr_shape[0]; j++) {
            for (int k=0; k < this->score_arr_shape[1]; k++) {
                score_array[j * this->score_arr_shape[1] *
                    this->score_arr_shape[2] + k *
                    this->score_arr_shape[2] + i] = *raw_score_ptr;
                raw_score_ptr++;
            }
        }
    }
    this->num_positions = this->score_arr_shape[0];
}



/// Convenience function to retrieve the number of positions in the cterm
/// template.
int CTermFinder::get_num_positions(void) {
    return this->num_positions;
}


/// Convenience function to retrieve the threshold score below which
/// it is very unlikely the retrieved position is actually a cterminal.
double CTermFinder::get_threshold_score(void) {
    return THRESHOLD_SCORE;
}



/// Takes the input scores and positions and merges them
/// so that the greater of the two light chain scores (K,L)
/// is retained. The light chain value will always go in position
/// 1, the heavy chain in 0.
void CTermFinder::merge_light_chain_scores(const std::array<double, 3> &scores,
        const std::array<int, 3> &positions,
        std::array<double, 2> &merged_scores,
        std::array<int, 2> &merged_positions) {
    for (size_t i = 0; i < 2; i ++) {
        merged_scores[i] = 0;
        merged_positions[i] = 0;
    }

    for (size_t i = 0; i < this->boundary_chains.size(); i++) {
        if (this->boundary_chains[i] == "H") {
            merged_positions[0] = positions[i];
            merged_scores[0] = scores[i];
        } else {
            if (scores[i] > merged_scores[1]) {
                merged_scores[1] = scores[i];
                merged_positions[1] = positions[i];
            }
        }
    }
}




// Does a fast (relative to typical alignments) search for templates
// which indicate the c-terminal of a kappa, lambda or heavy chain.
int CTermFinder::find_best_cterminal(std::string &query_sequence,
        std::array<double, 3> &best_scores, std::array<int, 3> &best_positions) {

    std::string error_message;

    for (int k=0; k < 3; k++) {
        best_scores[k] = 0;
        best_positions[k] = 0;
    }

    if (query_sequence.length() < MINIMUM_SEQUENCE_LENGTH ||
            (query_sequence.length() - this->num_positions) <= 0)
        return CTERM_FINDER_ERROR;

    auto encoded_sequence = std::make_unique<int[]>(query_sequence.length());

    if (!SequenceUtilities::convert_x_sequence_to_array(encoded_sequence.get(),
                query_sequence))
        return CTERM_FINDER_ERROR;

    // Look for the first two cysteines in the sequence,
    // then start searching for template matches which occur
    // after this.
    size_t startPosition = 0, ncysteines = 0;

    for (size_t i=0; i < query_sequence.length(); i++) {
        if (encoded_sequence[i] == 1) {
            ncysteines += 1;
            if (ncysteines >= 2) {
                startPosition = i + 1;
                break;
            }
        }
    }

    // If we could not find two cysteines, there is
    // something wrong; abort.
    if ( ncysteines < 2 )
        return INVALID_SEQUENCE;


    // We loop up to query_sequence.length() - 5, which is <
    // this->num_positions, in case there is a truncated c-
    // terminal at the end of the sequence. 5 is arbitrary
    // but past a certain point looking for say the first 3
    // letters of the c-terminal would be prone to false
    // positives.
    for (size_t i = 0; i < query_sequence.length() - 5; i++) {
        double match_score[3];
        size_t score_arr_row = 0;

        for (int k=0; k < 3; k++)
            match_score[k] = 0;

        for (int j=0; j < this->num_positions; j++) {
            // This is a little clunky and introduces unnecessary
            // branching. TODO: Break into two loops so we don't
            // have to do this.
            if (i + j >= query_sequence.length())
                break;

            for (int k=0; k < 3; k++) {
                size_t position = score_arr_row + encoded_sequence[i+j] *
                    this->score_arr_shape[2] + k;
                match_score[k] += this->score_array[position];
            }
            score_arr_row += (this->score_arr_shape[1] *
                    this->score_arr_shape[2]);
        }

        for (int k=0; k < 3; k++) {
            if (match_score[k] > best_scores[k]) {
                best_scores[k] = match_score[k];
                best_positions[k] = i;
            }
        }
    }

    return CTERM_FINDER_SUCCESS;
}






// A light wrapper on find_c_terminals for use by Python callers.
int CTermFinder::pyfind_c_terminals(std::string query_sequence,
        nb::ndarray<double, nb::shape<3>, nb::device::cpu, nb::c_contig> best_scores,
        nb::ndarray<int32_t, nb::shape<3>, nb::device::cpu, nb::c_contig> best_positions) {
    std::array<double, 3> scores;
    std::array<int, 3> positions;

    int err_code = this->find_best_cterminal(query_sequence, scores, positions);

    for (size_t i=0; i < 3; i++) {
        best_scores(i) = scores[i];
        best_positions(i) = positions[i];
    }
    return err_code;
}

}  // namespace NumberingTools
