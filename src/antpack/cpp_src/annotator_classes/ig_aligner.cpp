/* The IGAligner runs a highly efficient profile alignment for mAb sequences.
 * Copyright (C) 2025 Jonathan Parkinson
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
#include <vector>

// Library headers

// Project headers
#include "ig_aligner.h"

namespace NumberingTools {

IGAligner::IGAligner(
             std::string consensus_filepath,
             std::string chain_name,
             std::string scheme):
chain_name(chain_name),
scheme(scheme) {
    // Note that exceptions thrown here are sent back to Python via
    // PyBind as long as this constructor is used within the wrapper.
    if (chain_name != "H" && chain_name != "K" && chain_name != "L") {
        throw std::runtime_error(std::string("The IGAligner class currently "
                    "support only H, K, L chains."));
    }

    std::filesystem::path extension_path = consensus_filepath;
    extension_path = extension_path / "mabs";
    std::string uppercaseScheme = scheme;
    for (auto & c : uppercaseScheme) c = toupper(c);

    std::string npyFName = uppercaseScheme + "_CONSENSUS_" +
        chain_name + ".npy";
    std::string consFName = uppercaseScheme + "_CONSENSUS_" +
        chain_name + ".txt";
    std::filesystem::path npyFPath = extension_path / npyFName;
    std::filesystem::path consFPath = extension_path / consFName;

    std::vector<std::vector<std::string>> position_consensus;
    if (!cnpy::read_consensus_file(consFPath, position_consensus))
        throw std::runtime_error(std::string("Error in library installation (IG)"));

    cnpy::NpyArray raw_score_arr = cnpy::npy_load(npyFPath.string());
    double *raw_score_ptr = raw_score_arr.data<double>();
    if (raw_score_arr.word_size != 8)
        throw std::runtime_error(std::string("Error in library installation (IG)."));

    this->score_arr_shape[0] = raw_score_arr.shape[0];
    this->score_arr_shape[1] = raw_score_arr.shape[1];
    if (this->score_arr_shape[1] != EXPECTED_NUM_AAS_FOR_ALIGNERS)
        throw std::runtime_error(std::string("Error in library installation (IG)."));

    this->score_array = std::make_unique<double[]>(this->score_arr_shape[0] *
        this->score_arr_shape[1]);

    for (size_t k=0; k < raw_score_arr.shape[0] * raw_score_arr.shape[1]; k++)
        this->score_array[k] = raw_score_ptr[k];

    // Check that the number of positions in the scoring and
    // consensus is as expected, and set the list of highly conserved
    // positions according to the selected scheme.
    if (scheme == "imgt") {
        this->highly_conserved_positions = {HIGHLY_CONSERVED_IMGT_1,
            HIGHLY_CONSERVED_IMGT_2, HIGHLY_CONSERVED_IMGT_3,
            HIGHLY_CONSERVED_IMGT_4, HIGHLY_CONSERVED_IMGT_5,
            HIGHLY_CONSERVED_IMGT_6};
        if (this->score_arr_shape[0] != NUM_HEAVY_IMGT_POSITIONS &&
                this->score_arr_shape[0] != NUM_LIGHT_IMGT_POSITIONS) {
            throw std::runtime_error(std::string("The score_array passed to "
                "IGAligner must have the expected number of positions "
                "for the numbering system."));
        }
        if (position_consensus.size() != NUM_HEAVY_IMGT_POSITIONS &&
                position_consensus.size() != NUM_LIGHT_IMGT_POSITIONS) {
            throw std::runtime_error(std::string("The consensus sequence "
                "passed to IGAligner must have the expected number of "
                "positions for the numbering system."));
        }
    } else if (scheme == "aho") {
        this->highly_conserved_positions = {HIGHLY_CONSERVED_AHO_1,
            HIGHLY_CONSERVED_AHO_2, HIGHLY_CONSERVED_AHO_3,
            HIGHLY_CONSERVED_AHO_4, HIGHLY_CONSERVED_AHO_5,
            HIGHLY_CONSERVED_AHO_6};
        if (this->score_arr_shape[0] != NUM_HEAVY_AHO_POSITIONS &&
                this->score_arr_shape[0] != NUM_LIGHT_AHO_POSITIONS) {
            throw std::runtime_error(std::string("The score_array passed to "
                "IGAligner must have the expected number of positions "
                "for the numbering system."));
        }
        if (position_consensus.size() != NUM_HEAVY_AHO_POSITIONS &&
                position_consensus.size() != NUM_LIGHT_AHO_POSITIONS) {
            throw std::runtime_error(std::string("The consensus sequence "
                "passed to IGAligner must have the expected number of "
                "positions for the numbering system."));
        }
    } else if (scheme == "martin" || scheme == "kabat") {
        if (chain_name == "L" || chain_name == "K") {
            this->highly_conserved_positions = {HIGHLY_CONSERVED_KABAT_LIGHT_1,
                HIGHLY_CONSERVED_KABAT_LIGHT_2, HIGHLY_CONSERVED_KABAT_LIGHT_3,
                HIGHLY_CONSERVED_KABAT_LIGHT_4, HIGHLY_CONSERVED_KABAT_LIGHT_5,
                HIGHLY_CONSERVED_KABAT_LIGHT_6};
        } else if (chain_name == "H") {
            this->highly_conserved_positions = {HIGHLY_CONSERVED_KABAT_HEAVY_1,
                HIGHLY_CONSERVED_KABAT_HEAVY_2, HIGHLY_CONSERVED_KABAT_HEAVY_3,
                HIGHLY_CONSERVED_KABAT_HEAVY_4, HIGHLY_CONSERVED_KABAT_HEAVY_5,
                HIGHLY_CONSERVED_KABAT_HEAVY_6};
        }

        if (this->score_arr_shape[0] != NUM_HEAVY_MARTIN_KABAT_POSITIONS &&
                this->score_arr_shape[0] != NUM_LIGHT_MARTIN_KABAT_POSITIONS) {
            throw std::runtime_error(std::string("The score_array passed to "
                "IGAligner must have the expected number of positions "
                "for the numbering system."));
        }
        if (position_consensus.size() != NUM_HEAVY_MARTIN_KABAT_POSITIONS &&
                position_consensus.size() != NUM_LIGHT_MARTIN_KABAT_POSITIONS) {
            throw std::runtime_error(std::string("The consensus sequence "
                        "passed to IGAligner must have the expected number of "
                        "positions for the numbering system."));
        }
    } else {
        throw std::runtime_error(std::string("Currently IGAligner "
                "only recognizes schemes 'martin', 'kabat', 'imgt', 'aho'."));
    }

    this->num_positions = this->score_arr_shape[0];
    // Update the consensus map to indicate which letters are expected at
    // which positions. An empty set at a given position means that
    // ANY letter at that position is considered acceptable. These are
    // excluded from percent identity calculation. num_restricted_positions
    // indicates how many are INCLUDED in percent identity.
    this->num_restricted_positions = 0;

    for (size_t i = 0; i < position_consensus.size(); i++) {
        consensus_map.push_back(std::set<char> {});
        if (position_consensus[i].empty())
            continue;
        num_restricted_positions += 1;
        // This is a new addition in v0.3.5 -- treat X as an allowed character
        // at all positions!
        consensus_map[i].insert('X');

        for (std::string & AA : position_consensus[i]) {
            if (!SequenceUtilities::validate_sequence(AA)) {
                throw std::runtime_error(std::string("Non-default AAs supplied "
                    "in a consensus file."));
            }
            if (AA.length() > 1) {
                throw std::runtime_error(std::string("Non-default AAs supplied "
                    "in a consensus file."));
            }
            consensus_map[i].insert(AA[0]);
        }
    }
}


// Convenience function for retrieving the chain name.
std::string IGAligner::get_chain_name() {
    return this->chain_name;
}





/// @brief Numbers an input sequence
/// @param query_sequence Sequence to number
/// @param encoded_sequence Pointer to array of same length
///        as query sequence with the encoded query seq.
/// @param final_numbering vector that will store the generated
///        numbering
/// @param percent_identity the percent identity to the template.
/// @param error_message the error message (if none, "").
void IGAligner::align(std::string query_sequence,
            int *encoded_sequence, std::vector<std::string> &final_numbering,
            double &percent_identity, std::string &error_message) {

    allowedErrorCodes error_code;
    percent_identity = 0;


    if (query_sequence.length() < MINIMUM_SEQUENCE_LENGTH) {
        error_code = invalidSequence;
        error_message = this->error_code_to_message[error_code];
        return;
    }

    std::vector<int> position_key;
    position_key.reserve(query_sequence.length());

    int row_size = query_sequence.length() + 1;
    int numElements = row_size * (this->num_positions + 1);

    auto path_trace = std::make_unique<uint8_t[]>(numElements);
    auto init_numbering = std::make_unique<unsigned int[]>(this->num_positions);

    // Fill in the scoring table.
    fill_needle_scoring_table(path_trace.get(), query_sequence.length(),
        row_size, encoded_sequence, numElements);

    // Set all members of init_numbering array to 0.
    for (int i=0; i < this->num_positions; i++)
        init_numbering[i] = 0;

    // Backtrace the path, and determine how many insertions there
    // were at each position where there is an insertion. Determine
    // how many gaps exist at the beginning and end of the sequence
    // as well (the "N-terminus" and "C-terminus").
    int i = this->num_positions;
    int j = query_sequence.length();
    int numNTermGaps = 0;
    int numCTermGaps = 0;

    while (i > 0 || j > 0) {
        int grid_pos = i * row_size + j;
        switch (path_trace[grid_pos]) {
            case LEFT_TRANSFER:
                if (i > 0 && i < this->num_positions)
                    init_numbering[i-1] += 1;
                else if (i == 0)
                    numNTermGaps += 1;
                else
                    numCTermGaps += 1;
                j -= 1;
                break;
            case UP_TRANSFER:
                i -= 1;
                break;
            case DIAGONAL_TRANSFER:
                init_numbering[i-1] += 1;
                i -= 1;
                j -= 1;
                break;
            // This should never happen -- if it does, report it as a
            // SERIOUS error.
            default:
                error_code = fatalRuntimeError;
                error_message = this->error_code_to_message[error_code];
                return;
                break;
        }
    }

    // Add gaps at the beginning of the numbering corresponding to the
    // number of N-terminal insertions. This ensures the numbering has
    // the same length as the input sequence.
    // position_key is set to -1 to indicate these positions should not
    // be considered for percent identity calculation.
    for (i=0; i < numNTermGaps; i++) {
        final_numbering.push_back("-");
        position_key.push_back(-1);
    }


    // Build vector of IMGT numbers. Unfortunately the IMGT system adds insertion codes
    // forwards then backwards where there is > 1 insertion at a given position,
    // but ONLY if the insertion is at an expected position in the CDRs!
    // Everywhere else, insertions are recognized by just placing in
    // the usual order. This annoying quirk adds some complications.
    // We also create the position_key here which we can quickly use to determine
    // percent identity by referencing this->consensus_map, which indicates
    // what AAs are expected at each position.
    if (this->scheme == "imgt") {
        for (i=0; i < this->num_positions; i++) {
            if (init_numbering[i] == 0)
                continue;

            // For IMGT, set a maximum of 70 expected insertions anywhere.
            if (init_numbering[i] > this->alphabet.size() ||
                    init_numbering[i] > 70) {
                error_code = tooManyInsertions;
                error_message = this->error_code_to_message[error_code];
                return;
            }
            int ceil_cutpoint, floor_cutpoint;

            // These are the positions for which we need to deploy forward-backward
            // lettering (because the IMGT system is weird, but also the default...)
            // position_key is set to -1 for all insertions to indicate these should
            // not be considered for calculating percent identity.
            switch (i) {
                case CDR1_INSERTION_PT:
                case CDR2_INSERTION_PT:
                    ceil_cutpoint = init_numbering[i] / 2;
                    floor_cutpoint = (init_numbering[i] - 1) / 2;
                    for (j=0; j < floor_cutpoint; j++) {
                        final_numbering.push_back(std::to_string(i) + this->alphabet[j]);
                        position_key.push_back(-1);
                    }
                    for (j=ceil_cutpoint; j > 0; j--) {
                        final_numbering.push_back(std::to_string(i+1) + this->alphabet[j-1]);
                        position_key.push_back(-1);
                    }

                    final_numbering.push_back(std::to_string(i+1));
                    position_key.push_back(i);
                    break;

                case CDR3_INSERTION_PT:
                    final_numbering.push_back(std::to_string(i+1));
                    position_key.push_back(i);
                    ceil_cutpoint = init_numbering[i] / 2;
                    floor_cutpoint = (init_numbering[i] - 1) / 2;
                    for (j=0; j < floor_cutpoint; j++) {
                        final_numbering.push_back(std::to_string(i+1) + this->alphabet[j]);
                        position_key.push_back(-1);
                    }
                    for (j=ceil_cutpoint; j > 0; j--) {
                        final_numbering.push_back(std::to_string(i+2) + this->alphabet[j-1]);
                        position_key.push_back(-1);
                    }
                    break;
                default:
                    final_numbering.push_back(std::to_string(i+1));
                    position_key.push_back(i);
                    for (size_t k=1; k < init_numbering[i]; k++) {
                        final_numbering.push_back(std::to_string(i+1) +
                                this->alphabet[k-1]);
                        position_key.push_back(-1);
                    }
                    break;
            }
        }
    } else {
    // Build vector of numbers for other schemes (much simpler).
        for (i=0; i < this->num_positions; i++) {
            if (init_numbering[i] == 0)
                continue;

            if (init_numbering[i] > this->alphabet.size()) {
                error_code = tooManyInsertions;
                error_message = this->error_code_to_message[error_code];
                return;
            }
            final_numbering.push_back(std::to_string(i+1));
            position_key.push_back(i);
            for (size_t k=1; k < init_numbering[i]; k++) {
                final_numbering.push_back(std::to_string(i+1) +
                        this->alphabet[k-1]);
                position_key.push_back(-1);
            }
        }
    }


    // Add gaps at the end of the numbering corresponding to the
    // number of C-terminal insertions. This ensures the numbering has
    // the same length as the input sequence.
    for (int k=0; k < numCTermGaps; k++) {
        final_numbering.push_back("-");
        position_key.push_back(-1);
    }


    // position_key is now the same length as final_numbering and indicates at
    // each position whether that position maps to a standard
    // position and if so what. Now use this to calculate percent identity,
    // excluding those positions at which the numbering system tolerates any
    // amino acid (i.e. CDRs). At the same time, we can check whether the
    // expected amino acids are present at the highly conserved residues.
    // The consensus map will have only one possible amino acid at
    // those positions.

    int num_required_positions_found = 0;

    for (size_t k=0; k < position_key.size(); k++) {
        int scheme_std_position = position_key[k];
        if (scheme_std_position < 0)
            continue;

        if (this->consensus_map[scheme_std_position].empty())
            continue;

        if (this->consensus_map[scheme_std_position].find(query_sequence[k]) !=
                this->consensus_map[scheme_std_position].end()) {
            percent_identity += 1;
            if (std::binary_search(this->highly_conserved_positions.begin(),
                    this->highly_conserved_positions.end(),
                    scheme_std_position)) {
                num_required_positions_found += 1;
            }
        }
    }

    percent_identity /= this->num_restricted_positions;

    error_code = noError;
    // Check to make sure query sequence and numbering are same length. If not,
    // this is a fatal error. (This should be extremely rare -- we have not
    // encountered it so far).
    if (query_sequence.length() != final_numbering.size()) {
        error_code = alignmentWrongLength;
        error_message = this->error_code_to_message[error_code];
        return;
    }

    if (num_required_positions_found != 6)
        error_code = unacceptableConservedPositions;

    error_message = this->error_code_to_message[error_code];
}





/// @brief Fills in the scoring table to construct an alignment.
/// @param path_trace The array that will store the best path found (aka
/// the scoring table).
/// @param query_seq_len The length of the query sequence.
/// @param encoded_sequence A pointer to the array containing the encoded
/// sequence.
/// @param num_elements The number of elements in the scoring table.
void IGAligner::fill_needle_scoring_table(uint8_t *path_trace,
                int query_seq_len, int row_size, const int *encoded_sequence,
                int &numElements) {
    auto needle_scores = std::make_unique<double[]>( numElements );
    double score_updates[3], best_score;
    int grid_pos, diag_neighbor, upper_neighbor, best_pos;
    double *score_itr = this->score_array.get();

    // Fill in the first row of the tables. We use a default score here
    // to ensure that insertions in the template only are highly tolerated
    // at the beginning of the sequence.
    needle_scores[0] = 0;
    path_trace[0] = LEFT_TRANSFER;
    for (int i=0; i < query_seq_len; i++) {
        double updated_score = needle_scores[i] + NTERMINAL_QUERY_GAP_PENALTY;
        // This seems counterintuitive but realize that the gap penalty is
        // NEGATIVE.
        updated_score = updated_score > MAX_NTERMINAL_GAP_PENALTY ?
            updated_score : MAX_NTERMINAL_GAP_PENALTY;
        needle_scores[i+1] = updated_score;
        path_trace[i+1] = LEFT_TRANSFER;
    }



    // This first loop goes up to the last row only (the last row)
    // is handled separately since it requires special logic).
    // Handling the last row separately is about 10 - 15% faster than
    // using if else statements to determine when the last row is
    // reached inside this loop.
    for (int i=1; i < this->num_positions; i++) {
        grid_pos = i * row_size;
        diag_neighbor = (i - 1) * row_size;
        upper_neighbor = diag_neighbor + 1;
        needle_scores[grid_pos] = needle_scores[diag_neighbor] +
            score_itr[QUERY_GAP_COLUMN];
        path_trace[grid_pos] = UP_TRANSFER;
        grid_pos++;

        for (int j=0; j < (query_seq_len - 1); j++) {
            score_updates[DIAGONAL_TRANSFER] = needle_scores[diag_neighbor] +
                score_itr[encoded_sequence[j]];
            score_updates[LEFT_TRANSFER] = needle_scores[grid_pos - 1] +
                score_itr[TEMPLATE_GAP_COLUMN];
            score_updates[UP_TRANSFER] = needle_scores[upper_neighbor] +
                score_itr[QUERY_GAP_COLUMN];

            // This is mildly naughty -- we don't consider the possibility
            // of a tie, which could lead to a branched alignment. Realistically,
            // ties are going to be rare and low-impact in this situation because
            // of the way the positions are scored. For now we are defaulting to
            // "match" if there is a tie.
            best_score = score_updates[DIAGONAL_TRANSFER];
            best_pos = DIAGONAL_TRANSFER;
            if (score_updates[LEFT_TRANSFER] > best_score) {
                best_pos = LEFT_TRANSFER;
                best_score = score_updates[LEFT_TRANSFER];
            }
            if (score_updates[UP_TRANSFER] > best_score) {
                best_pos = UP_TRANSFER;
                best_score = score_updates[UP_TRANSFER];
            }
            needle_scores[grid_pos] = best_score;
            path_trace[grid_pos] = best_pos;

            grid_pos++;
            diag_neighbor++;
            upper_neighbor++;
        }

        score_updates[DIAGONAL_TRANSFER] = needle_scores[diag_neighbor] +
            score_itr[encoded_sequence[query_seq_len - 1]];
        score_updates[LEFT_TRANSFER] = needle_scores[grid_pos - 1] +
            score_itr[TEMPLATE_GAP_COLUMN];
        // We use a default score for the last column. We exclude however
        // highly conserved positions, for which a large gap penalty should
        // always be assigned.
            if (std::binary_search(this->highly_conserved_positions.begin(),
                    this->highly_conserved_positions.end(), i))
                score_updates[UP_TRANSFER] = needle_scores[upper_neighbor] +
                    score_itr[QUERY_GAP_COLUMN];
            else
                score_updates[UP_TRANSFER] = needle_scores[upper_neighbor] +
                    CTERMINAL_TEMPLATE_GAP_PENALTY;

            best_score = score_updates[DIAGONAL_TRANSFER];
            best_pos = DIAGONAL_TRANSFER;
            if (score_updates[LEFT_TRANSFER] > best_score) {
                best_pos = LEFT_TRANSFER;
                best_score = score_updates[LEFT_TRANSFER];
            }
            if (score_updates[UP_TRANSFER] > best_score) {
                best_pos = UP_TRANSFER;
                best_score = score_updates[UP_TRANSFER];
            }
            needle_scores[grid_pos] = best_score;
            path_trace[grid_pos] = best_pos;

            score_itr += this->score_arr_shape[1];
        }



        // Now handle the last row.
        score_itr = this->score_array.get() + (this->num_positions-1) *
            this->score_arr_shape[1];
        grid_pos = (this->num_positions) * row_size;
        diag_neighbor = (this->num_positions - 1) * row_size;
        upper_neighbor = diag_neighbor + 1;

        needle_scores[grid_pos] = needle_scores[diag_neighbor] +
            score_itr[QUERY_GAP_COLUMN];
        path_trace[grid_pos] = UP_TRANSFER;
        grid_pos++;

        for (int j=0; j < (query_seq_len - 1); j++) {
            score_updates[DIAGONAL_TRANSFER] = needle_scores[diag_neighbor] +
                score_itr[encoded_sequence[j]];
            // If a gap is already open, extending a c-terminal gap should
            // cost nothing.
            if (path_trace[grid_pos - 1] != LEFT_TRANSFER) {
                score_updates[LEFT_TRANSFER] = needle_scores[grid_pos - 1] +
                    CTERMINAL_QUERY_GAP_PENALTY;
            } else {
                score_updates[LEFT_TRANSFER] = needle_scores[grid_pos - 1];
            }
            score_updates[UP_TRANSFER] = needle_scores[upper_neighbor] +
                score_itr[QUERY_GAP_COLUMN];

            best_score = score_updates[DIAGONAL_TRANSFER];
            best_pos = DIAGONAL_TRANSFER;
            if (score_updates[LEFT_TRANSFER] > best_score) {
                best_pos = LEFT_TRANSFER;
                best_score = score_updates[LEFT_TRANSFER];
            }
            if (score_updates[UP_TRANSFER] > best_score) {
                best_pos = UP_TRANSFER;
                best_score = score_updates[UP_TRANSFER];
            }
            needle_scores[grid_pos] = best_score;
            path_trace[grid_pos] = best_pos;

            grid_pos++;
            diag_neighbor++;
            upper_neighbor++;
        }



        // And, finally, the last column of the last row.
        score_updates[DIAGONAL_TRANSFER] = needle_scores[diag_neighbor] +
            score_itr[encoded_sequence[query_seq_len-1]];
        // If a gap is already open, extending a c-terminal gap should
        // cost nothing.
        if (path_trace[grid_pos - 1] != LEFT_TRANSFER) {
            score_updates[LEFT_TRANSFER] = needle_scores[grid_pos - 1] +
                    CTERMINAL_QUERY_GAP_PENALTY;
        } else {
            score_updates[LEFT_TRANSFER] = needle_scores[grid_pos - 1];
        }
        score_updates[UP_TRANSFER] = needle_scores[upper_neighbor] +
            CTERMINAL_TEMPLATE_GAP_PENALTY;

        best_score = score_updates[DIAGONAL_TRANSFER];
        best_pos = DIAGONAL_TRANSFER;
        if (score_updates[LEFT_TRANSFER] > best_score) {
            best_pos = LEFT_TRANSFER;
            best_score = score_updates[LEFT_TRANSFER];
        }
        if (score_updates[UP_TRANSFER] > best_score) {
            best_pos = UP_TRANSFER;
            best_score = score_updates[UP_TRANSFER];
        }
        needle_scores[grid_pos] = best_score;
        path_trace[grid_pos] = best_pos;
}

}  // namespace NumberingTools
