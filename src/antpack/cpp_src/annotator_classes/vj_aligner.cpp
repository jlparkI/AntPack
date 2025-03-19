/* The VJAligner class handles alignment of an input sequence by
 * finding the most appropriate V and J genes from a list then
 * using custom preconstructed scoring matrices. This is different
 * from IGAligner, which aligns to a profile (not individual VJ
 * genes).
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
#include <unordered_map>

// Library headers

// Project headers
#include "vj_aligner.h"


namespace NumberingTools {


// Enforces that the V and J gene sizes we see when loading
// library files are what we expect, and the number of gaps
// that we expect to see at the beginning of a jgene.
static constexpr int EXPECTED_JGENE_SIZE = 16;
static constexpr int EXPECTED_VGENE_SIZE = 112;
static constexpr int NUM_JGENE_GAPS = 5;
static constexpr int EXPECTED_BLOSUM_SHAPE = 22;

// The size of the sliding windows that we want to extract
// from v and jgenes and their locations.
static constexpr int VJ_WINDOW_SIZE = 11;
static constexpr int VGENE_WINDOW1_START = 40;
static constexpr int VGENE_WINDOW2_START = 93;
static constexpr int VGENE_WINDOW3_START = 12;


VJAligner::VJAligner(std::string consensus_filepath,
    std::string chain_name,
    std::string scheme,
    std::string receptor_type):
chain_name(chain_name),
scheme(scheme) {
    // Note that exceptions thrown here are sent back to Python via
    // PyBind as long as this constructor is used within the wrapper.
    // Currently this is set up to be used only with TCRs, but we may
    // modify this if this design turns out to work well.
    if (chain_name != "A" && chain_name != "B" && chain_name != "D" &&
            chain_name != "G") {
        throw std::runtime_error(std::string("VJAligner currently "
                    "support only A, B, D, G chains."));
    }

    std::filesystem::path extension_path = consensus_filepath;
    // First load the BLOSUM matrix.
    extension_path = extension_path / "mabs" / "blosum_matrix.npy";
    cnpy::NpyArray raw_blosum_arr = cnpy::npy_load(extension_path.string());
    int16_t *raw_blosum_ptr = raw_blosum_arr.data<int16_t>();
    if (raw_blosum_arr.word_size != 2 || raw_blosum_arr.shape[0] !=
            EXPECTED_BLOSUM_SHAPE || raw_blosum_arr.shape[1] !=
            EXPECTED_BLOSUM_SHAPE)
        throw std::runtime_error(std::string("Error in library installation "
                    "(VJ)."));

    this->blosum_array = std::make_unique<int32_t[]>(raw_blosum_arr.shape[0] *
            raw_blosum_arr.shape[1]);

    for (size_t i=0; i < EXPECTED_BLOSUM_SHAPE; i++) {
        for (size_t j=0; j < EXPECTED_BLOSUM_SHAPE; j++)
            this->blosum_array[i*EXPECTED_BLOSUM_SHAPE + j] =
                static_cast<int32_t>(raw_blosum_ptr[i*EXPECTED_BLOSUM_SHAPE + j]);
    }


    extension_path = consensus_filepath;
    extension_path = extension_path / "tcrs";
    std::string uppercase_scheme = scheme;
    for (auto & c : uppercase_scheme) c = std::toupper(c);

    // Load jgene and vgene sequence files.
    std::string current_filename = receptor_type + "_" +
        chain_name + "J.txt";
    std::filesystem::path current_filepath = extension_path / current_filename;
    if (!cnpy::read_tcr_vj_gene_file(current_filepath, this->jgenes,
                EXPECTED_JGENE_SIZE)) {
        throw std::runtime_error(std::string("Error in library installation "
                    "(wrong jgene sequence sizes)."));
    }

    current_filename = receptor_type + "_" + chain_name + "V.txt";
    current_filepath = extension_path / current_filename;
    if (!cnpy::read_tcr_vj_gene_file(current_filepath, this->vgenes,
                EXPECTED_VGENE_SIZE)) {
        throw std::runtime_error(std::string("Error in library installation "
                    "(wrong vgene sizes)."));
    }

    // Next, encode selected sliding windows from the j-genes and
    // from the v-genes. We will blosum score the sliding windows on
    // the input data to determine which v- and j-gene most closely
    // resembles a query sequence.
    this->encoded_v_window1 = std::make_unique<int[]>(this->vgenes.size() *
                VJ_WINDOW_SIZE);
    this->encoded_v_window2 = std::make_unique<int[]>(this->vgenes.size() *
                VJ_WINDOW_SIZE);
    this->encoded_v_window3 = std::make_unique<int[]>(this->vgenes.size() *
                VJ_WINDOW_SIZE);
    this->encoded_j_window1 = std::make_unique<int[]>(this->jgenes.size() *
                VJ_WINDOW_SIZE);

    for (size_t i=0; i < this->vgenes.size(); i++) {
        if (this->vgenes.at(i).length() < VGENE_WINDOW1_START +
                VJ_WINDOW_SIZE) {
            throw std::runtime_error(std::string(
                        "Error in library installation, 1A."));
        }
        if (this->vgenes.at(i).length() < VGENE_WINDOW2_START +
                VJ_WINDOW_SIZE || this->vgenes.at(i).length() <
                VGENE_WINDOW3_START) {
            throw std::runtime_error(std::string(
                        "Error in library installation, 1B."));
        }

        std::string window = this->vgenes.at(i).substr(VGENE_WINDOW1_START,
                VJ_WINDOW_SIZE);
        if (!SequenceUtilities::convert_gapped_x_sequence_to_array(
            this->encoded_v_window1.get() + i*VJ_WINDOW_SIZE, window)) {
            throw std::runtime_error(std::string(
                        "Error in library installation, 1C."));
        }

        window = this->vgenes.at(i).substr(VGENE_WINDOW2_START,
                VJ_WINDOW_SIZE);
        if (!SequenceUtilities::convert_gapped_x_sequence_to_array(
            this->encoded_v_window2.get() + i*VJ_WINDOW_SIZE, window)) {
            throw std::runtime_error(std::string(
                        "Error in library installation, 1D."));
        }

        window = this->vgenes.at(i).substr(VGENE_WINDOW3_START,
                VJ_WINDOW_SIZE);
        if (!SequenceUtilities::convert_gapped_x_sequence_to_array(
            this->encoded_v_window3.get() + i*VJ_WINDOW_SIZE, window)) {
            throw std::runtime_error(std::string(
                        "Error in library installation, 1D1."));
        }
    }

    for (size_t i=0; i < this->jgenes.size(); i++) {
        if (this->jgenes.at(i).length() < VJ_WINDOW_SIZE) {
            throw std::runtime_error(std::string(
                        "Error in library installation, 1E."));
        }

        std::string window = this->jgenes.at(i).substr(
                this->jgenes.at(i).length() - VJ_WINDOW_SIZE,
                VJ_WINDOW_SIZE);

        if (!SequenceUtilities::convert_gapped_x_sequence_to_array(
            this->encoded_j_window1.get() +
            i * VJ_WINDOW_SIZE, window)) {
            throw std::runtime_error(std::string(
                        "Error in library installation, 1F."));
        }
    }

    // Load vgene and jgene scoring matrices. These contain scores for all V and
    // j genes for a specific chain type. Contents are transferred to unique_ptr
    // score matrices.
    //
    // VGENES
    //
    current_filename = receptor_type + "_" + chain_name + "V.npy";
    current_filepath = extension_path / current_filename;
    cnpy::NpyArray vraw_score_arr = cnpy::npy_load(current_filepath.string());
    double *raw_score_ptr = vraw_score_arr.data<double>();
    if (vraw_score_arr.word_size != 8)
        throw std::runtime_error(std::string("Error in library installation "
                    "1G."));

    for (size_t i=0; i < 3; i++)
        this->vgene_score_arr_shape[i] = vraw_score_arr.shape[i];

    if (this->vgene_score_arr_shape[1] != EXPECTED_VGENE_SIZE ||
            this->vgene_score_arr_shape[2] != EXPECTED_NUM_AAS_FOR_ALIGNERS ||
            this->vgene_score_arr_shape[0] != this->vgenes.size())
        throw std::runtime_error(std::string("Error in library installation, "
                    "1H."));

    this->vgene_score_array = std::make_unique<double[]>(
            this->vgene_score_arr_shape[0] * this->vgene_score_arr_shape[1] *
            this->vgene_score_arr_shape[2]);

    for (size_t k=0; k < vraw_score_arr.shape[0] *
            vraw_score_arr.shape[1] * vraw_score_arr.shape[2]; k++)
        this->vgene_score_array[k] = raw_score_ptr[k];
    //
    // JGENES
    //
    current_filename = receptor_type + "_" + chain_name + "J.npy";
    current_filepath = extension_path / current_filename;
    cnpy::NpyArray jraw_score_arr = cnpy::npy_load(current_filepath.string());
    raw_score_ptr = jraw_score_arr.data<double>();
    if (jraw_score_arr.word_size != 8)
        throw std::runtime_error(std::string("Error in library installation, "
                    "1I."));

    for (size_t i=0; i < 3; i++)
        this->jgene_score_arr_shape[i] = jraw_score_arr.shape[i];

    if (this->jgene_score_arr_shape[1] != EXPECTED_JGENE_SIZE ||
            this->jgene_score_arr_shape[2] != EXPECTED_NUM_AAS_FOR_ALIGNERS ||
            this->jgene_score_arr_shape[0] != this->jgenes.size())
        throw std::runtime_error(std::string("Error in library installation, "
                    "1J."));

    this->jgene_score_array = std::make_unique<double[]>(
            this->jgene_score_arr_shape[0] * this->jgene_score_arr_shape[1] *
            this->jgene_score_arr_shape[2]);

    for (size_t k=0; k < jraw_score_arr.shape[0] *
            jraw_score_arr.shape[1] * jraw_score_arr.shape[2]; k++)
        this->jgene_score_array[k] = raw_score_ptr[k];


    // Currently only the IMGT scheme is supported. This aligner is primarily
    // used for TCRs, which are not in any case compatible with Martin / Kabat.
    // It may be helpful to add Aho, but interconversion between IMGT / Aho
    // is much more straightforward than interconversion between IMGT and Martin
    // or Kabat, so if Aho is added the Annotator class should make this
    // conversion after VJAligner has set up the alignment.
    if (scheme == "imgt") {
        this->highly_conserved_positions = {HIGHLY_CONSERVED_IMGT_1,
            HIGHLY_CONSERVED_IMGT_2, HIGHLY_CONSERVED_IMGT_3,
            HIGHLY_CONSERVED_IMGT_4, HIGHLY_CONSERVED_IMGT_5,
            HIGHLY_CONSERVED_IMGT_6};
    } else {
        throw std::runtime_error(std::string("Currently VJAligner "
                "only recognizes the IMGT scheme."));
    }
}


// Convenience function for retrieving the chain name.
std::string VJAligner::get_chain_name() {
    return this->chain_name;
}


/// @brief Identifies the vgene which has the most kmers in
/// common with an input sequence.
/// @param query_sequence The input sequence.
/// @param encoded_sequence Pointer to an array containing the
/// input sequence encoded as integers; must be same size as
/// query_sequence.length().
/// @param identity The score for matching to the best vgene;
/// the result is stored in this reference.
/// @param best_vgene_number The id of the best vgene that is
/// found; the result is stored in this reference.
/// @return Returns 1 (VALID_SEQUENCE) or 0 for an error.
int VJAligner::identify_best_vgene(std::string &query_sequence,
        const int *encoded_sequence, int &identity, int &best_vgene_number) {
    if (query_sequence.length() < VJ_WINDOW_SIZE + 1 ||
            query_sequence.length() < VJ_WINDOW_SIZE + 2)
        return INVALID_SEQUENCE;

    std::vector<int32_t> vgene_scores(this->vgenes.size(), 0);
    int *vwindow1 = this->encoded_v_window1.get();
    int *vwindow2 = this->encoded_v_window2.get();
    int *vwindow3 = this->encoded_v_window3.get();
    int32_t *blosum_mat = this->blosum_array.get();

    for (size_t i=0; i < this->vgenes.size(); i++) {
        int32_t best_scores[3] = {0, 0, 0};
        const int *query_ptr = encoded_sequence;
        size_t vletters1[VJ_WINDOW_SIZE],
            vletters2[VJ_WINDOW_SIZE],
            vletters3[VJ_WINDOW_SIZE];

        for (size_t k=0; k < VJ_WINDOW_SIZE; k++) {
            vletters1[k] = *vwindow1 * EXPECTED_BLOSUM_SHAPE;
            vletters2[k] = *vwindow2 * EXPECTED_BLOSUM_SHAPE;
            vletters3[k] = *vwindow3 * EXPECTED_BLOSUM_SHAPE;
            vwindow1++;
            vwindow2++;
            vwindow3++;
        }

        for (size_t j=0; j < query_sequence.length() -
                VJ_WINDOW_SIZE; j++) {
            int32_t scores[3] = {0, 0, 0};

            for (size_t k=0; k < VJ_WINDOW_SIZE; k++) {
                size_t qletter = query_ptr[k];
                scores[0] += blosum_mat[vletters1[k] +
                    qletter];
                scores[1] += blosum_mat[vletters2[k] +
                    qletter];
                scores[2] += blosum_mat[vletters3[k] +
                    qletter];
            }

            query_ptr++;
            for (size_t k=0; k < 3; k++) {
                if (scores[k] > best_scores[k])
                    best_scores[k] = scores[k];
            }
        }
        vgene_scores[i] += best_scores[0] + best_scores[1] + best_scores[2];
    }


    identity = 0;
    best_vgene_number = 0;
    for (size_t i=0; i < vgene_scores.size(); i++) {
        if (vgene_scores[i] > identity) {
            identity = vgene_scores[i];
            best_vgene_number = i;
        }
    }
    return VALID_SEQUENCE;
}


/// @brief Identifies the jgene which is the best match for
/// an input sequence based on a sliding window, and finds
/// the location in the query sequence where the jgene should
/// likely be aligned.
/// @param query_sequence The input sequence.
/// @param encoded_sequence Pointer to an array containing the
/// input sequence encoded as integers; must be same size as
/// query_sequence.length().
/// @param optimal_position The position at which the jgene
/// alignment should most likely occur; the result is stored
/// in this reference.
/// @param identity The score for matching to the best jgene;
/// the result is stored in this reference.
/// @param best_jgene_number The id of the best jgene that is
/// found; the result is stored in this reference.
/// @return Returns 1 (VALID_SEQUENCE) or 0 for an error.
int VJAligner::identify_best_jgene(std::string &query_sequence,
        const int *encoded_sequence, int &optimal_position, int &identity,
        int &best_jgene_number) {
    if (query_sequence.length() < VJ_WINDOW_SIZE + 2)
        return INVALID_SEQUENCE;

    std::vector<int32_t> jgene_scores(this->jgenes.size(), 0);
    std::vector<int> optimal_positions(this->jgenes.size(), 0);
    int *jwindow = this->encoded_j_window1.get();

    for (size_t i=0; i < this->jgenes.size(); i++) {
        int32_t best_score = 0;
        const int *query_ptr = encoded_sequence;
        int jletters[VJ_WINDOW_SIZE];
        for (size_t k=0; k < VJ_WINDOW_SIZE; k++)
            jletters[k] = jwindow[k] * EXPECTED_BLOSUM_SHAPE;

        // Loop up to the length of the sequence minus 7. 7 is arbitrary
        // but past a certain point if we look for matches to say only
        // the first three positions of the jgene we will start to get
        // false positives.
        for (size_t j=0; j < query_sequence.length() - 7; j++) {
            int cap = (j < query_sequence.length() - VJ_WINDOW_SIZE) ?
                VJ_WINDOW_SIZE : (query_sequence.length() - j);
            int32_t score = 0;

            for (size_t k=0; k < cap; k++) {
                score += this->blosum_array[jletters[k] +
                    query_ptr[k]];
            }
            query_ptr++;

            if (score > best_score) {
                optimal_positions[i] = j;
                best_score = score;
            }
        }

        jgene_scores[i] = best_score;
        jwindow += VJ_WINDOW_SIZE;
    }

    identity = 0;
    best_jgene_number = 0;
    for (size_t i=0; i < jgene_scores.size(); i++) {
        if (jgene_scores[i] > identity) {
            identity = jgene_scores[i];
            best_jgene_number = i;
            optimal_position = optimal_positions[i];
        }
    }
    return VALID_SEQUENCE;
}



/// @brief Aligns an input sequence to a specified template V and J gene.
/// @param query_sequence The sequence to be aligned.
/// @param encoded_sequence A pointer to an array of ints containing the
/// sequence encoded as numbers.
/// @param final_number The vector in which the final numbering (as a vector
/// of strings) will be stored.
/// @param percent_identity The variable in which the percent identity to
/// the template vgene will be stored.
/// @param error_message The variable in which the error message for the
/// alignment (if any) will be stored.
/// @param vgene_number The number of the vgene. Determines which vgene in
/// this class's list of vgenes will be used. If the input value is unacceptable
/// this will cause an exception to be thrown.
/// @param jgene_number The number of the jgene. Determines which jgene in
/// this class's list of jgenes will be used. If the input value is unacceptable
/// this will cause an exception to be thrown.
void VJAligner::align(std::string query_sequence,
            int *encoded_sequence, std::vector<std::string> &final_numbering,
            double &percent_identity, std::string &error_message,
            int vgene_number, int jgene_number) {

    allowedErrorCodes error_code;
    percent_identity = 0;


    if (query_sequence.length() < MINIMUM_SEQUENCE_LENGTH) {
        error_code = invalidSequence;
        error_message = this->error_code_to_message[error_code];
        return;
    }
    if (vgene_number >= this->vgenes.size() ||
            jgene_number >= this->jgenes.size() ||
            vgene_number < 0 || jgene_number < 0) {
        throw std::runtime_error("Internal error; invalid v or jgene assigned. "
                "Please report.");
    }

    std::vector<int> position_key;
    position_key.reserve(query_sequence.length());

    int row_size = query_sequence.length() + 1;
    int num_elements = row_size * (this->num_positions + 1);

    auto path_trace = std::make_unique<uint8_t[]>(num_elements);
    auto init_numbering = std::make_unique<unsigned int[]>(this->num_positions);

    // Fill in the scoring table.
    fill_needle_scoring_table(path_trace.get(), query_sequence.length(),
        row_size, encoded_sequence, num_elements,
        vgene_number, jgene_number);

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


    // Build vector of IMGT numbers. Unfortunately the IMGT system adds
    // insertion codes forwards then backwards where there is > 1 insertion
    // at a given position, but ONLY if the insertion is at an expected
    // position in the CDRs! Everywhere else, insertions are recognized by
    // just placing in the usual order. This annoying quirk adds complications.
    // We also create the position_key here which we can quickly use to
    // determine percent identity by matching it to positions in the
    // vgene.
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

        // These are the positions for which we need to deploy
        // forward-backward lettering (because the IMGT system is
        // weird, but also the default...) position_key is set to
        // -1 for all insertions to indicate these should
        // not be considered for calculating percent identity.
        switch (i) {
            case CDR1_INSERTION_PT:
            case CDR2_INSERTION_PT:
                ceil_cutpoint = init_numbering[i] / 2;
                floor_cutpoint = (init_numbering[i] - 1) / 2;
                for (j=0; j < floor_cutpoint; j++) {
                    final_numbering.push_back(std::to_string(i) +
                            this->alphabet[j]);
                    position_key.push_back(-1);
                }
                for (j=ceil_cutpoint; j > 0; j--) {
                    final_numbering.push_back(std::to_string(i+1) +
                            this->alphabet[j-1]);
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
                    final_numbering.push_back(std::to_string(i+1) +
                            this->alphabet[j]);
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


    // Add gaps at the end of the numbering corresponding to the
    // number of C-terminal insertions. This ensures the numbering has
    // the same length as the input sequence.
    for (int k=0; k < numCTermGaps; k++) {
        final_numbering.push_back("-");
        position_key.push_back(-1);
    }


    // position_key is now the same length as final_numbering and indicates at
    // each position whether that position maps to a standard
    // position and if so what. Now use this to calculate percent identity.
    // In an important difference from IGAligner, we calculate percent identity
    // here relative to the assigned Vgene and Jgene rather than to a shared
    // template, and we can ignore any positions which are gaps in the
    // V and Jgene (since these are already gapped).
    // At the same time, we can check whether the expected amino acids are
    // present at the highly conserved residues. The consensus map will have
    // only one possible amino acid at those positions.

    int num_required_positions_found = 0;
    int num_restricted_positions = 0;
    for (const auto & vl : this->vgenes[vgene_number]) {
        if (vl != '-')
            num_restricted_positions += 1;
    }
    for (const auto & jl : this->jgenes[jgene_number]) {
        if (jl != '-')
            num_restricted_positions += 1;
    }

    for (size_t k=0; k < position_key.size(); k++) {
        int scheme_std_position = position_key[k];

        if (scheme_std_position < 0) {
            continue;
        } else if (scheme_std_position < EXPECTED_VGENE_SIZE) {
            if (this->vgenes[vgene_number][scheme_std_position] ==
                    query_sequence[k])
                percent_identity += 1;
        } else if (scheme_std_position < this->num_positions) {
            if (this->jgenes[jgene_number][scheme_std_position -
                    EXPECTED_VGENE_SIZE] == query_sequence[k])
                percent_identity += 1;
        }

        if (std::binary_search(this->highly_conserved_positions.begin(),
            this->highly_conserved_positions.end(), scheme_std_position)) {
            num_required_positions_found += 1;
        }
    }

    percent_identity /= static_cast<double>(num_restricted_positions);

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





/// Fill in the scoring table created by caller, using the position-specific
/// scores for indels and substitutions, and add the appropriate
/// pathway traces.
/// @brief Fills in the scoring table to construct an alignment.
/// @param path_trace The array that will store the best path found (aka
/// the scoring table).
/// @param query_seq_len The length of the query sequence.
/// @param encoded_sequence A pointer to the array containing the encoded
/// sequence.
/// @param num_elements The number of elements in the scoring table.
/// @param vgene_number The vgene we should use.
/// @param jgene_number The jgene we should use.
void VJAligner::fill_needle_scoring_table(uint8_t *path_trace,
                int query_seq_len, int row_size, const int *encoded_sequence,
                const int &num_elements, const int &vgene_number,
                const int &jgene_number) {
    auto needle_scores = std::make_unique<double[]>( num_elements );
    double score_updates[3], best_score;
    int grid_pos, diag_neighbor, upper_neighbor, best_pos;
    // The size of the V and J gene score arrays is checked by the class
    // constructor. We set the pointer to the row corresponding to the
    // selected v and j gene (caller should check this is valid).
    double *vgene_itr = this->vgene_score_array.get() +
        vgene_number * EXPECTED_VGENE_SIZE * EXPECTED_NUM_AAS_FOR_ALIGNERS;
    double *jgene_itr = this->jgene_score_array.get() +
        jgene_number * EXPECTED_JGENE_SIZE * EXPECTED_NUM_AAS_FOR_ALIGNERS;

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
    // It also skips the last EXPECTED_JGENE_SIZE columns since these
    // are filled by a separate loop which will use the jgene scoring
    // table (this->num_positions = EXPECTED_VGENE_SIZE + EXPECTED_JGENE_SIZE).
    for (int i=0; i < EXPECTED_VGENE_SIZE; i++) {
        // Add 1 to i here to skip the first row which is already filled.
        grid_pos = (i+1) * row_size;
        diag_neighbor = i * row_size;
        upper_neighbor = diag_neighbor + 1;
        needle_scores[grid_pos] = needle_scores[diag_neighbor] +
            vgene_itr[QUERY_GAP_COLUMN];
        path_trace[grid_pos] = UP_TRANSFER;
        grid_pos++;

        for (int j=0; j < (query_seq_len - 1); j++) {
            score_updates[DIAGONAL_TRANSFER] = needle_scores[diag_neighbor] +
                vgene_itr[encoded_sequence[j]];
            score_updates[LEFT_TRANSFER] = needle_scores[grid_pos - 1] +
                vgene_itr[TEMPLATE_GAP_COLUMN];
            score_updates[UP_TRANSFER] = needle_scores[upper_neighbor] +
                vgene_itr[QUERY_GAP_COLUMN];

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
            vgene_itr[encoded_sequence[query_seq_len - 1]];
        score_updates[LEFT_TRANSFER] = needle_scores[grid_pos - 1] +
            vgene_itr[TEMPLATE_GAP_COLUMN];
        // We use a default score for the last column. We exclude however
        // highly conserved positions, for which a large gap penalty should
        // always be assigned. Add 1 to i here to skip first row.
        if (std::binary_search(this->highly_conserved_positions.begin(),
                this->highly_conserved_positions.end(), (i+1)))
            score_updates[UP_TRANSFER] = needle_scores[upper_neighbor] +
                vgene_itr[QUERY_GAP_COLUMN];
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

        vgene_itr += this->vgene_score_arr_shape[2];
    }

    //
    // Now handle the next EXPECTED_JGENE_SIZE rows before handling
    // the last row. The last row is handled separately due to
    // special gap conditions.

    // Add 1 and EXPECTED_VGENE_SIZE to grid_row here since the first
    // row and EXPECTED_VGENE_SIZE rows are already filled. We use
    // grid_row to keep track of rows since i will be off by
    // EXPECTED_VGENE_SIZE.
    int grid_row = EXPECTED_VGENE_SIZE + 1;

    for (int i=0; i < (EXPECTED_JGENE_SIZE - 1); i++) {
        grid_pos = grid_row * row_size;
        diag_neighbor = (grid_row - 1) * row_size;
        upper_neighbor = diag_neighbor + 1;
        needle_scores[grid_pos] = needle_scores[diag_neighbor] +
            jgene_itr[QUERY_GAP_COLUMN];
        path_trace[grid_pos] = UP_TRANSFER;
        grid_pos++;

        for (int j=0; j < (query_seq_len - 1); j++) {
            score_updates[DIAGONAL_TRANSFER] = needle_scores[diag_neighbor] +
                jgene_itr[encoded_sequence[j]];
            score_updates[LEFT_TRANSFER] = needle_scores[grid_pos - 1] +
                jgene_itr[TEMPLATE_GAP_COLUMN];
            score_updates[UP_TRANSFER] = needle_scores[upper_neighbor] +
                jgene_itr[QUERY_GAP_COLUMN];

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
            jgene_itr[encoded_sequence[query_seq_len - 1]];
        score_updates[LEFT_TRANSFER] = needle_scores[grid_pos - 1] +
            jgene_itr[TEMPLATE_GAP_COLUMN];
        // We use a default score for the last column. We exclude however
        // highly conserved positions, for which a large gap penalty should
        // always be assigned.
        if (std::binary_search(this->highly_conserved_positions.begin(),
                this->highly_conserved_positions.end(), grid_row))
            score_updates[UP_TRANSFER] = needle_scores[upper_neighbor] +
                jgene_itr[QUERY_GAP_COLUMN];
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

        jgene_itr += this->jgene_score_arr_shape[2];
        grid_row++;
    }


    // Now handle the last row.
    grid_pos = (this->num_positions) * row_size;
    diag_neighbor = (this->num_positions - 1) * row_size;
    upper_neighbor = diag_neighbor + 1;

    needle_scores[grid_pos] = needle_scores[diag_neighbor] +
        jgene_itr[QUERY_GAP_COLUMN];
    path_trace[grid_pos] = UP_TRANSFER;
    grid_pos++;

    for (int j=0; j < (query_seq_len - 1); j++) {
        score_updates[DIAGONAL_TRANSFER] = needle_scores[diag_neighbor] +
            jgene_itr[encoded_sequence[j]];
        // If a gap is already open, extending a c-terminal gap should
        // cost minimally.
        if (path_trace[grid_pos - 1] != LEFT_TRANSFER) {
            score_updates[LEFT_TRANSFER] = needle_scores[grid_pos - 1] +
                CTERMINAL_QUERY_GAP_PENALTY;
        } else {
            score_updates[LEFT_TRANSFER] = needle_scores[grid_pos - 1];
        }
        score_updates[UP_TRANSFER] = needle_scores[upper_neighbor] +
            jgene_itr[QUERY_GAP_COLUMN];

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
        jgene_itr[encoded_sequence[query_seq_len-1]];
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
