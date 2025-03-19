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
#include "prefiltering_tool.h"

// C++ headers
#include <vector>
#include <memory>
#include <string>
#include <unordered_map>


namespace PrefilteringRoutines {


PrefilteringTool::PrefilteringTool(std::string consensus_filepath,
        std::unordered_map<std::string, size_t> nterm_kmers) {

    std::filesystem::path extension_path = consensus_filepath;
    extension_path = extension_path / "mabs";
    std::string npyFName = "CTERMFINDER_CONSENSUS_H.npy";
    std::filesystem::path npyFPath = extension_path / npyFName;
    cnpy::NpyArray raw_score_arr;

    try {
        raw_score_arr = cnpy::npy_load(npyFPath.string());
    }
    catch (...) {
        throw std::runtime_error(std::string("The consensus file / "
                    "library installation has an issue."));
    }
    this->ig_nterm_kmers = nterm_kmers;

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
        std::filesystem::path npyFPath = extension_path / npyFName;

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

    // Check the k-mer map to make sure all kmers are expected size
    // and map to an int in [0,1,2].
    for (auto & map_pair : nterm_kmers) {
        if (map_pair.first.length() != KMER_SIZE) {
            throw std::runtime_error(std::string("The consensus file / "
                        "library installation has an issue."));
        }
        if (map_pair.second > 2) {
            throw std::runtime_error(std::string("The consensus file / "
                        "library installation has an issue."));
        }
    }
}

/// Convenience function to retrieve the ordering of chains used by
/// this class.
std::vector<std::string> PrefilteringTool::get_chain_list(void) {
    return this->boundary_chains;
}


/// Convenience function to retrieve the number of positions in the cterm
/// template.
int PrefilteringTool::get_num_positions(void) {
    return this->num_positions;
}


/// Convenience function to retrieve the threshold score below which
/// it is very unlikely the retrieved position is actually a cterminal.
double PrefilteringTool::get_cterm_threshold_score(void) {
    return CTERM_THRESHOLD_SCORE;
}


/// Convenience function to retrieve the threshold score below which
/// it is very unlikely the retrieved position is actually a
/// nterm region.
int PrefilteringTool::get_nterm_threshold_score(void) {
    return START_THRESHOLD_SCORE;
}


/// Convenience function to retrieve the kmer window size.
int PrefilteringTool::get_kmer_window_size(void) {
    return KMER_WINDOW_WIDTH;
}



/// Takes the input scores and positions and merges them
/// so that the greater of the two light chain scores (K,L)
/// is retained. The light chain value will always go in position
/// 1, the heavy chain in 0.
void PrefilteringTool::merge_light_chain_cterm_scores(
        const std::array<double, 3> &scores,
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


/// Takes the input scores and positions and merges them
/// so that the greater of the two light chain scores (K,L)
/// is retained. The light chain value will always go in position
/// 1, the heavy chain in 0.
void PrefilteringTool::merge_light_chain_nterm_scores(
        const std::array<int, 3> &scores,
        const std::array<int, 3> &positions,
        std::array<int, 2> &merged_scores,
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
// which indicate the c-terminal or n-terminal of a light or heavy
// chain.
int PrefilteringTool::find_start_end_zones(std::string &query_sequence,
        std::array<double, 3> &best_cterm_scores,
        std::array<int, 3> &best_cterm_positions,
        std::array<int, 3> &best_nterm_scores,
        std::array<int, 3> &best_nterm_positions) {

    std::string error_message;

    for (int k=0; k < 3; k++) {
        best_cterm_scores[k] = 0;
        best_cterm_positions[k] = 0;
        best_nterm_scores[k] = 0;
        best_nterm_positions[k] = 0;
    }

    // Possible overkill, but good to be sure. Check that the sequence is
    // larger than any loop we will run and larger than the minimum
    // sequence we consider acceptable.
    if (query_sequence.length() < NumberingTools::MINIMUM_SEQUENCE_LENGTH ||
            query_sequence.length() <= this->num_positions ||
            query_sequence.length() <= KMER_SIZE ||
            query_sequence.length() <= KMER_WINDOW_WIDTH)
        return CTERM_FINDER_ERROR;

    std::array<std::vector<size_t>, 3> kmer_locations;
    for (size_t i=0; i < 3; i++) {
        std::vector<size_t> blanks(query_sequence.length(), 0);
        kmer_locations[i] = blanks;
    }
    auto encoded_sequence = std::make_unique<int[]>(query_sequence.length());

    if (!SequenceUtilities::convert_x_sequence_to_array(encoded_sequence.get(),
                query_sequence))
        return CTERM_FINDER_ERROR;

    for (size_t i=0; i < (query_sequence.length() - KMER_SIZE + 1); i++) {
        std::string kmer = query_sequence.substr(i, KMER_SIZE);
        if (this->ig_nterm_kmers.count(kmer) > 0) {
            size_t kmer_type = this->ig_nterm_kmers.at(kmer);
            kmer_locations[kmer_type][i] = 1;
        }
    }


    // We loop up to query_sequence.length() - 5, which is <
    // this->num_positions, in case there is a truncated c-
    // terminal at the end of the sequence. 5 is arbitrary
    // but past a certain point looking for say the first 3
    // letters of the c-terminal would be prone to false
    // positives.
    for (size_t i = 0; i < query_sequence.length() - 5; i++) {
        double cterm_match_score[3];
        size_t score_arr_row = 0;

        for (int k=0; k < 3; k++)
            cterm_match_score[k] = 0;

        for (int j=0; j < this->num_positions; j++) {
            // This is a little clunky and introduces unnecessary
            // branching. TODO: Break into two loops so we don't
            // have to do this.
            if (i + j >= query_sequence.length())
                break;

            for (int k=0; k < 3; k++) {
                size_t position = score_arr_row + encoded_sequence[i+j] *
                    this->score_arr_shape[2] + k;
                cterm_match_score[k] += this->score_array[position];
            }
            score_arr_row += (this->score_arr_shape[1] *
                    this->score_arr_shape[2]);
        }

        for (int k=0; k < 3; k++) {
            if (cterm_match_score[k] > best_cterm_scores[k]) {
                best_cterm_scores[k] = cterm_match_score[k];
                best_cterm_positions[k] = i;
            }
        }
    }

    // Next, run a separate loop for starting positions. TODO:
    // Merge this with the cterminal finder loop for greater
    // efficiency.
    size_t num_matching_kmers[3];

    for (size_t i=0; i < 3; i++)
        num_matching_kmers[i] = 0;

    for (size_t i=0; i < KMER_WINDOW_WIDTH; i++) {
        for (size_t k=0; k < 3; k++)
            num_matching_kmers[k] += kmer_locations[k][i];
    }

    for (size_t i=0; i < 3; i++)
        best_nterm_scores[i] = num_matching_kmers[i];

    for (size_t i=1; i < (query_sequence.length() -
                KMER_WINDOW_WIDTH + 1); i++) {
        for (size_t k=0; k < 3; k++) {
            num_matching_kmers[k] -= kmer_locations[k][i-1];
            num_matching_kmers[k] += kmer_locations[k][i+KMER_WINDOW_WIDTH-1];
            if (num_matching_kmers[k] > best_nterm_scores[k]) {
                best_nterm_scores[k] = num_matching_kmers[k];
                best_nterm_positions[k] = i;
            }
        }
    }

    return CTERM_FINDER_SUCCESS;
}






// A light wrapper on find_c_terminals for use by Python callers.
int PrefilteringTool::pyfind_c_terminals(std::string query_sequence,
        nb::ndarray<double, nb::shape<3>, nb::device::cpu, nb::c_contig> best_cterm_scores,
        nb::ndarray<int32_t, nb::shape<3>, nb::device::cpu, nb::c_contig> best_cterm_positions,
        nb::ndarray<int32_t, nb::shape<3>, nb::device::cpu, nb::c_contig> best_nterm_scores,
        nb::ndarray<int32_t, nb::shape<3>, nb::device::cpu, nb::c_contig> best_nterm_positions) {
    std::array<double, 3> cterm_scores;
    std::array<int, 3> cterm_positions, nterm_scores, nterm_positions;

    int err_code = this->find_start_end_zones(query_sequence,
            cterm_scores, cterm_positions,
            nterm_scores, nterm_positions);

    for (size_t i=0; i < 3; i++) {
        best_cterm_scores(i) = cterm_scores[i];
        best_cterm_positions(i) = cterm_positions[i];
        best_nterm_scores(i) = nterm_scores[i];
        best_nterm_positions(i) = nterm_positions[i];
    }
    return err_code;
}

}  // namespace PrefilteringRoutines
