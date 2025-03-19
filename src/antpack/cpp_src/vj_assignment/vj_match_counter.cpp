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
#include <memory>
#include <tuple>
#include <vector>
#include <string>
#include <map>
#include <utility>

// Library headers

// Project headers
#include "vj_match_counter.h"


namespace VJAssignment {

VJMatchCounter::VJMatchCounter(std::map<std::string,
        std::vector<std::string>> gene_names,
        std::map<std::string, std::vector<std::string>> gene_seqs,
        nb::ndarray<double, nb::shape<22,22>, nb::device::cpu, nb::c_contig> blosum_matrix,
        std::string scheme,
        std::string consensus_filepath):
gene_seqs(gene_seqs),
gene_names(gene_names),
blosum_matrix(blosum_matrix),
scheme(scheme) {
    // Note that exceptions thrown here go back to Python via
    // PyBind as long as this constructor is used within the wrapper.
    for ( const auto &gene_seq_element : gene_seqs ) {
        std::map<std::string, std::vector<std::string>>::iterator
            gene_name_element = gene_names.find(gene_seq_element.first);

        if (gene_name_element == gene_names.end()) {
            throw std::runtime_error(std::string("One or more keys "
                        "in the gene_seqs dictionary is not in the "
                        "gene_names dictionary."));
        }

        if (gene_name_element->second.size() !=
                gene_seq_element.second.size()) {
            throw std::runtime_error(std::string("The gene name lists "
                        "and gene sequence lists must be of the same sizes."));
        }
        for (const auto &gene_sequence : gene_seq_element.second) {
            if (gene_sequence.length() != REQUIRED_SEQUENCE_LENGTH) {
                throw std::runtime_error(std::string("All sequences passed "
                            "to VJMatchCounter must have the correct length."));
            }
        }
    }

    for ( const auto &gene_name_element : gene_names ) {
        if (gene_seqs.find(gene_name_element.first) == gene_seqs.end()) {
            throw std::runtime_error(std::string("One or more keys in the "
                        "gene_names dictionary is not in the gene_seqs "
                        "dictionary."));
        }
        std::map<std::string, int> name_submap;
        for (size_t i=0; i < gene_name_element.second.size(); i++)
            name_submap[gene_name_element.second[i]] = i;

        this->names_to_positions[gene_name_element.first] =
            std::move(name_submap);
    }

    // Initialize the set of expected imgt positions.
    for (size_t i=1; i < REQUIRED_SEQUENCE_LENGTH + 1; i++)
        this->essential_imgt_map.insert(std::to_string(i));

    // If we are using a scheme OTHER than IMGT, we will
    // map the key positions needed for VJ matching in this
    // other scheme to IMGT using a function from utilities. For now,
    // check to make sure scheme is accepted.
    if (this->scheme != "imgt" && this->scheme != "aho"
            && this->scheme != "kabat" && this->scheme != "martin") {
        throw std::runtime_error(std::string("Unrecognized scheme supplied. "
                    "Scheme must be one of imgt, martin, kabat or aho."));
    }

    // Set up three imgt scheme aligners for internal use (useful for various
    // tasks).
    this->h_aligner = std::make_unique<NumberingTools::IGAligner>
        (consensus_filepath, "H", "imgt");
    this->k_aligner = std::make_unique<NumberingTools::IGAligner>
        (consensus_filepath, "K", "imgt");
    this->l_aligner = std::make_unique<NumberingTools::IGAligner>
        (consensus_filepath, "L", "imgt");
}



/// @brief Assigns the correct VJ genes for a specified species -- or,
/// if unknown is specified, tries to determine the correct species.
/// This function is exposed to Python and to external classes.
/// @param alignment A tuple of numbering, percent identity, chain and
/// error message such as is returned by any of the Annotators.
/// @param sequence The sequence to be checked.
/// @param species A supported species OR "unknown". If "unknown",
/// all available species will be checked.
/// @param mode One of "identity", "evalue". Determines how the closest
/// VJ gene is selected.
/// @return Returns a tuple of v-genes (separated by _), j-genes (separated
/// by _), score for closest v-genes, score for closest j-genes, species.
/// More than one v and j gene is reported when multiple v-genes or j-genes
/// tie for best score. Species is the same as the input species UNLESS
/// the user has supplied unknown, in which case the best identified species
/// is reported.
std::tuple<std::string, std::string, double, double, std::string>
VJMatchCounter::assign_vj_genes(std::tuple<std::vector<std::string>,
        double, std::string, std::string> alignment, std::string sequence,
        std::string species, std::string mode) {

    if (species == "unknown") {
        std::vector<std::string> subspecies_list = {"human", "mouse"};
        std::tuple<std::string, std::string,
                double, double, std::string> best_result =
                {"", "", 0, 0, "unknown"};

        for (const auto & subspecies : subspecies_list) {
            std::tuple<std::string, std::string,
                double, double, std::string> current_result =
                    this->internal_vj_assignment(alignment,
                            sequence, subspecies, mode);
            if (std::get<2>(current_result) > std::get<2>(best_result) &&
                    std::get<3>(current_result) > std::get<3>(best_result))
                best_result = current_result;
        }

        return best_result;
    } else if (species == "alpaca" && std::get<2>(alignment) != "H") {
        throw std::runtime_error(std::string("For alpacas, only heavy "
                    "chain sequences are allowed."));
    } else if (species != "human" && species != "mouse"
            && species != "alpaca" && species != "rabbit") {
        throw std::runtime_error(std::string("Species for VJ gene "
                    "assignment must be one of 'human', 'mouse', "
                    "'alpaca', 'rabbit', 'unknown'."));
    }
    return this->internal_vj_assignment(alignment, sequence, species, mode);
}



/// @brief Assigns the correct VJ genes for a specified species. This
/// function is wrapped by assign_vj_genes and is private so is not
/// exposed to external callers.
/// @param alignment A tuple of numbering, percent identity, chain and
/// error message such as is returned by any of the Annotators.
/// @param sequence The sequence to be checked.
/// @param species A supported species. "unknown" is not allowed.
/// @param mode One of "identity", "evalue". Determines how the closest
/// VJ gene is selected.
/// @return Returns a tuple of v-genes (separated by _), j-genes (separated
/// by _), score for closest v-genes, score for closest j-genes, species.
/// The species reported is the same as input (it is returned however for
/// ease of use by caller, which may need to report species when evaluating
/// multiple possibilities). More than one v and j gene is reported when
/// multiple v-genes or j-genes tie for best score.
std::tuple<std::string, std::string, double, double, std::string>
VJMatchCounter::internal_vj_assignment(std::tuple<std::vector<std::string>,
        double, std::string, std::string> &alignment,
        std::string &sequence,
        const std::string &species,
        const std::string &mode) {
    std::string vkey = species + "_IG" + std::get<2>(alignment) + "V";
    std::string jkey = species + "_IG" + std::get<2>(alignment) + "J";

    if (species != "human" && species != "mouse"
            && species != "alpaca" && species != "rabbit") {
        throw std::runtime_error(std::string("Species for VJ gene "
                    "assignment must be one of 'human', 'mouse'."));
    }

    if (std::get<2>(alignment) == "A" || std::get<2>(alignment) == "B" ||
            std::get<2>(alignment) == "D" || std::get<2>(alignment) == "G") {
        if (species != "human" && species != "mouse") {
            throw std::runtime_error(std::string("For TCRs, species "
                        "must be one of 'human', 'mouse'."));
        }
        if (this->scheme != "imgt") {
            throw std::runtime_error(std::string("For TCRs, the only "
                        "currently-supported scheme is IMGT."));
        }
        vkey = species + "_TR" + std::get<2>(alignment) + "V";
        jkey = species + "_TR" + std::get<2>(alignment) + "J";
    } else if (std::get<2>(alignment) != "H" && std::get<2>(alignment) != "K" &&
            std::get<2>(alignment) != "L") {
        return std::tuple<std::string, std::string,
                double, double, std::string>{"", "", 0, 0, "unknown"};
    }

    if (sequence.length() != std::get<0>(alignment).size()) {
        return std::tuple<std::string, std::string,
                double, double, std::string>{"", "", 0, 0, "unknown"};
    }
    if (!SequenceUtilities::validate_x_sequence(sequence)) {
        return std::tuple<std::string, std::string,
                double, double, std::string>{"", "", 0, 0, "unknown"};
    }


    std::map<std::string, std::vector<std::string>>::iterator
        vname_dict_finder = this->gene_names.find(vkey);
    std::map<std::string, std::vector<std::string>>::iterator
        jname_dict_finder = this->gene_names.find(jkey);
    std::map<std::string, std::vector<std::string>>::iterator
        vgene_dict_finder = this->gene_seqs.find(vkey);
    std::map<std::string, std::vector<std::string>>::iterator
        jgene_dict_finder = this->gene_seqs.find(jkey);

    if (vname_dict_finder == gene_names.end() ||
            jname_dict_finder == gene_names.end() ||
            vgene_dict_finder == gene_seqs.end() ||
            jgene_dict_finder == gene_seqs.end()) {
        throw std::runtime_error(std::string("There was an error "
                    "in the construction of the "
                    "vj gene tool; necessary gene databases "
                    "were not found. Please report."));
    }

    std::vector<std::string> &vgenes = vgene_dict_finder->second;
    std::vector<std::string> &jgenes = jgene_dict_finder->second;
    std::vector<std::string> &vnames = vname_dict_finder->second;
    std::vector<std::string> &jnames = jname_dict_finder->second;

    // This string will store the positions in the
    // sequence that correspond to IMGT 1-128.
    // If the scheme is not IMGT, we can convert the
    // numbering to IMGT temporarily for
    // this purpose.
    std::string prepped_sequence(REQUIRED_SEQUENCE_LENGTH, '-');
    int err_code = this->prep_sequence(prepped_sequence, sequence, alignment);
    if (err_code != VALID_SEQUENCE) {
        return std::tuple<std::string, std::string,
                double, double, std::string>{"", "", 0, 0, "unknown"};
    }

    if (vgenes.size() != vnames.size() || jgenes.size() != jnames.size()) {
        throw std::runtime_error(std::string("There was an error in "
                    "the construction of the vj gene tool; necessary "
                    "gene databases were not found. Please report."));
    }

    std::string vgene_name, jgene_name;
    double videntity, jidentity;

    if (mode == "identity") {
        this->assign_gene_by_identity(vgenes, vnames, prepped_sequence,
                videntity, vgene_name, 'v');
        this->assign_gene_by_identity(jgenes, jnames, prepped_sequence,
                jidentity, jgene_name, 'j');
    } else if (mode == "evalue") {
        auto encoded_query = std::make_unique<int[]>
            (REQUIRED_SEQUENCE_LENGTH);

        if (!SequenceUtilities::
                convert_gapped_x_sequence_to_array(encoded_query.get(),
                    prepped_sequence)) {
            throw std::runtime_error(std::string("The input sequence "
                        "contains invalid amino acids and could "
                        "not be encoded."));
        }

        if (!this->assign_gene_by_evalue(vgenes, vnames, encoded_query.get(),
                videntity, vgene_name, 'v')) {
            throw std::runtime_error(std::string("Invalid sequence supplied."));
        }
        if (!this->assign_gene_by_evalue(jgenes, jnames, encoded_query.get(),
                jidentity, jgene_name, 'j')) {
            throw std::runtime_error(std::string("Invalid sequence supplied."));
        }
    } else {
        throw std::runtime_error(std::string("Unrecognized mode was supplied. "
                    "Mode should be one of 'identity', 'evalue'."));
    }

    return std::tuple<std::string, std::string,
           double, double, std::string>{vgene_name, jgene_name,
                videntity, jidentity, species};
}



/// @brief Finds the closest matching genes to a query sequence that
/// has been converted to IMGT format with standard positions extracted,
/// using a supplied list of gene sequences and names, and determining
/// which is the best match using percent identity. This function is
/// private, not accessible to outside callers.
/// @param gene_seqs A list of V or J genes to which to match. These
/// should be length 128 and appropriately gapped to match IMGT format.
/// @param gene_names A list of gene names of the same length as gene_seqs.
/// @param prepped_sequence An input sequence prepped so that the standard
/// 128 imgt positions have been extracted.
/// @param best_identity The score for the best match is stored here.
/// @param best_gene_name The best gene name is stored here.
/// @param gene_type Indicates whether this is v or j.
void VJMatchCounter::assign_gene_by_identity(
std::vector<std::string> &gene_seqs,
std::vector<std::string> &gene_names,
std::string &prepped_sequence,
double &best_identity,
std::string &best_gene_name,
char gene_type) {
    // Sometimes with e-value or identity there is a tie where multiple genes
    // have the same e-value or percent identity. In this case, we report all
    // of them. To do so, we temporarily keep track of the best matches so far.
    // We reserve 10 which is an arbitrary number much larger than the number
    // of best matches we are likely to see in practice.
    std::vector<int> best_matches;
    best_matches.reserve(10);

    best_identity = 0;
    int start_letter = 0, end_letter = 0;

    if (gene_type == 'j') {
        start_letter = 105;
        end_letter = REQUIRED_SEQUENCE_LENGTH;
    } else {
        start_letter = 0;
        end_letter = 108;
    }

    for (size_t i=0; i < gene_seqs.size(); i++) {
        int matching_positions = 0;
        int nonzero_positions = 0;

        for (int j=start_letter; j < end_letter; j++) {
            if (gene_seqs[i][j] == '-')
                continue;
            nonzero_positions += 1;
            if (gene_seqs[i][j] == prepped_sequence[j] ||
                    prepped_sequence[j] == 'X')
                matching_positions += 1;
        }
        if (nonzero_positions == 0)
            nonzero_positions = 1;
        double current_identity = (static_cast<double>(matching_positions)) /
            (static_cast<double>(nonzero_positions));

        // Since we are comparing doubles here, just check if they are
        // close (for an arbitrary definition of close, here defined as
        // numbers within 1e-10).
        if (current_identity > (best_identity + 1e-10)) {
            best_matches.clear();
            best_matches.push_back(i);
            best_identity = current_identity;
        } else if (current_identity > (best_identity - 1e-10)) {
            best_matches.push_back(i);
        }
    }

    // If only one top match was found, we can use that name...
    // otherwise, use a comma-separated list.
    if (best_matches.size() == 1) {
        best_gene_name = gene_names[best_matches[0]];
    } else if (best_matches.size() > 1) {
        best_gene_name = gene_names[best_matches[0]];
        for (size_t i=1; i < best_matches.size(); i++)
            best_gene_name += "_" + gene_names[best_matches[i]];
    }
}




/// @brief Finds the closest matching genes to a query sequence that
/// has been converted to IMGT format with standard positions extracted,
/// using a supplied list of gene sequences and names, and determining
/// which is the best match using BLOSUM score, which is the same as
/// determining it by evalue. This function is private, not accessible
/// to outside callers.
/// @param gene_seqs A list of V or J genes to which to match. These
/// should be length 128 and appropriately gapped to match IMGT format.
/// @param gene_names A list of gene names of the same length as gene_seqs.
/// @param prepped_sequence An input sequence prepped so that the standard
/// 128 imgt positions have been extracted.
/// @param best_identity The score for the best match is stored here.
/// @param best_gene_name The best gene name is stored here.
/// @param gene_type Indicates whether this is v or j.
/// @return An error code to indicate success or failure.
int VJMatchCounter::assign_gene_by_evalue(
std::vector<std::string> &gene_seqs,
std::vector<std::string> &gene_names,
int *encoded_sequence,
double &best_identity,
std::string &best_gene_name,
char gene_type) {
    // Sometimes with e-value or identity there is a tie where multiple genes
    // have the same e-value or percent identity. In this case, we report all
    // of them. To do so, we temporarily keep track of the best matches so far.
    // We reserve 10 which is an arbitrary number much larger than the number
    // of best matches we are likely to see in practice.
    std::vector<int> best_matches;
    best_matches.reserve(10);

    best_gene_name = "";
    int64_t best_score = 0;

    // For a (very slight) speed gain we loop over only
    // the parts of the numbering that will not contain
    // any gaps.
    int start_letter = 0, end_letter = 0;
    if (gene_type == 'j') {
        start_letter = 105;
        end_letter = REQUIRED_SEQUENCE_LENGTH;
    } else {
        start_letter = 0;
        end_letter = 108;
    }

    for (size_t i=0; i < gene_seqs.size(); i++) {
        int64_t blosum_score = 0;

        for (int j=start_letter; j < end_letter; j++) {
            switch (gene_seqs[i][j]) {
                case 'A':
                    blosum_score += blosum_matrix(0, encoded_sequence[j]);
                    break;
                case 'C':
                    blosum_score += blosum_matrix(1, encoded_sequence[j]);
                    break;
                case 'D':
                    blosum_score += blosum_matrix(2, encoded_sequence[j]);
                    break;
                case 'E':
                    blosum_score += blosum_matrix(3, encoded_sequence[j]);
                    break;
                case 'F':
                    blosum_score += blosum_matrix(4, encoded_sequence[j]);
                    break;
                case 'G':
                    blosum_score += blosum_matrix(5, encoded_sequence[j]);
                    break;
                case 'H':
                    blosum_score += blosum_matrix(6, encoded_sequence[j]);
                    break;
                case 'I':
                    blosum_score += blosum_matrix(7, encoded_sequence[j]);
                    break;
                case 'K':
                    blosum_score += blosum_matrix(8, encoded_sequence[j]);
                    break;
                case 'L':
                    blosum_score += blosum_matrix(9, encoded_sequence[j]);
                    break;
                case 'M':
                    blosum_score += blosum_matrix(10, encoded_sequence[j]);
                    break;
                case 'N':
                    blosum_score += blosum_matrix(11, encoded_sequence[j]);
                    break;
                case 'P':
                    blosum_score += blosum_matrix(12, encoded_sequence[j]);
                    break;
                case 'Q':
                    blosum_score += blosum_matrix(13, encoded_sequence[j]);
                    break;
                case 'R':
                    blosum_score += blosum_matrix(14, encoded_sequence[j]);
                    break;
                case 'S':
                    blosum_score += blosum_matrix(15, encoded_sequence[j]);
                    break;
                case 'T':
                    blosum_score += blosum_matrix(16, encoded_sequence[j]);
                    break;
                case 'V':
                    blosum_score += blosum_matrix(17, encoded_sequence[j]);
                    break;
                case 'W':
                    blosum_score += blosum_matrix(18, encoded_sequence[j]);
                    break;
                case 'Y':
                    blosum_score += blosum_matrix(19, encoded_sequence[j]);
                    break;
                case '-':
                    break;
                default:
                    return INVALID_SEQUENCE;
                    break;
            }
        }

        if (blosum_score > best_score) {
            best_matches.clear();
            best_matches.push_back(i);
            best_score = blosum_score;
        } else if (blosum_score == best_score) {
            best_matches.push_back(i);
        }
    }

    // If only one top match was found, we can use that name...
    // otherwise, use a comma-separated list.
    if (best_matches.size() == 1) {
        best_gene_name = gene_names[best_matches[0]];
    } else if (best_matches.size() > 1) {
        best_gene_name = gene_names[best_matches[0]];
        for (size_t i=1; i < best_matches.size(); i++)
            best_gene_name += "_" + gene_names[best_matches[i]];
    }

    best_identity = best_score;

    return VALID_SEQUENCE;
}






/// @brief Converts a sequence to a standard format by extracting the
/// 128 standard IMGT positions for easy comparison to pre-gapped
/// V and J genes.
/// @param prepped_sequence The output is stored here.
/// @param sequence The input sequence.
/// @param alignment A tuple of numbering, percent identity, chain
/// and error message such as is returned by the Annotator classes.
/// @return An error code indicating success or failure.
int VJMatchCounter::prep_sequence(std::string &prepped_sequence,
std::string &sequence,
std::tuple<std::vector<std::string>, double,
std::string, std::string> &alignment) {
    if (this->scheme == "imgt") {
        std::vector<std::string> &numbering = std::get<0>(alignment);
        for (size_t i=0; i < sequence.length(); i++) {
            if (this->essential_imgt_map.find(numbering[i]) !=
                    this->essential_imgt_map.end())
                prepped_sequence.at(std::stoi(numbering[i]) - 1) = sequence[i];
        }
    } else {
        // If the scheme is not IMGT, we have several options.
        // We can convert numbering to IMGT, and indeed some programs
        // (ANARCI) rely heavily on interconversion between schemes,
        // although this involves some heuristics so for now we are
        // avoiding this. Alternatively we can realign the trimmed
        // sequence. This is actually quite fast because 1) we trim
        // the sequence and 2) we only need to do one chain alignment
        // (we already know the chain type). So for now we are doing
        // this (although it is not perfect either). The third option is
        // to do alignments of VJ genes against the original sequence without
        // using numbering, and while some tools use this this is slow,
        // so we are not doing that here.

        // Finally, it should be noted that for TCRs we currently only
        // support IMGT and thus if the chain type indicates TCR should
        // abort.
        if (std::get<2>(alignment) == "A" || std::get<2>(alignment) == "B" ||
                std::get<2>(alignment) == "D" ||
                std::get<2>(alignment) == "G") {
            throw std::runtime_error("For TCRs, only the IMGT scheme is "
                    "currently supported.");
        }
        std::vector<std::string> trimmed_numbering;
        std::vector<char> trimmed_seq_vector;
        int exstart, exend;

        // The only circumstance that could cause trim_alignment_utility
        // to return an error is already checked for by caller.
        int err_code = SequenceUtilities::trim_alignment_utility(sequence,
                alignment, trimmed_numbering,
                exstart, exend, trimmed_seq_vector);
        std::string trimmed_seq(trimmed_seq_vector.begin(),
                trimmed_seq_vector.end());

        auto queryAsIdx = std::make_unique<int[]>(trimmed_seq.length());
        std::vector<std::string> imgt_numbering;
        std::string errorMessage = "";
        double percent_identity = -1;
        imgt_numbering.reserve(trimmed_seq.length() * 2);

        if (!SequenceUtilities::convert_x_sequence_to_array(queryAsIdx.get(),
                    trimmed_seq))
            return INVALID_SEQUENCE;

        // Now that the alignment is trimmed, we can renumber the
        // trimmed sequence. Caller already checks that chain is one
        // of H, K, L so no need to recheck here. The order H = 0, K = 1,
        // L = 2 is set in the class constructor.
        if (std::get<2>(alignment) == "H") {
            this->h_aligner->align(trimmed_seq,
                    queryAsIdx.get(), imgt_numbering, percent_identity,
                    errorMessage);
        } else if (std::get<2>(alignment) == "K") {
            this->k_aligner->align(trimmed_seq,
                    queryAsIdx.get(), imgt_numbering, percent_identity,
                    errorMessage);
        } else if (std::get<2>(alignment) == "L") {
            this->l_aligner->align(trimmed_seq,
                    queryAsIdx.get(), imgt_numbering, percent_identity,
                    errorMessage);
        }

        if (imgt_numbering.size() != trimmed_seq.length())
            return INVALID_SEQUENCE;

        for (size_t i=0; i < trimmed_seq.length(); i++) {
            if (this->essential_imgt_map.find(imgt_numbering[i]) !=
                    this->essential_imgt_map.end())
                prepped_sequence.at(std::stoi(imgt_numbering[i]) - 1) =
                    sequence[i];
        }
    }
    return VALID_SEQUENCE;
}



/// @brief Gets the sequence of a named V or J gene for an indicated species.
/// @param query_name The name of the gene.
/// @param species The name of the species.
/// @return The sequence in gapped IMGT format. If the gene is not found,
/// this string is empty.
std::string VJMatchCounter::get_vj_gene_sequence(std::string query_name,
    std::string species) {
    size_t matchID;
    std::string output = "";

    if (query_name.length() < 4)
        return output;

    std::string dict_key = species + "_" + query_name.substr(0, 4);

    // We assume here that if the key is in names to positions,
    // it is also in gene_seqs, and that the vectors stored in
    // both are the same length, because of the way we initialized
    // names to positions (see class constructor).
    if (this->names_to_positions.find(dict_key) !=
            this->names_to_positions.end()) {
        std::map<std::string, int> &name_to_position_map =
            this->names_to_positions[dict_key];

        if (name_to_position_map.find(query_name) !=
                name_to_position_map.end()) {
            matchID = name_to_position_map[query_name];
            output = this->gene_seqs[dict_key][matchID];
        }
    }
    return output;
}



// This function just returns the maps stored by the object,
// which are converted to python dicts if used through the
// wrapper. This is used only for testing.
std::tuple<std::map<std::string, std::vector<std::string>>,
  std::map<std::string, std::vector<std::string>>>  VJMatchCounter::get_seq_lists() {
    return std::tuple<std::map<std::string, std::vector<std::string>>,
        std::map<std::string, std::vector<std::string>>>
        {this->gene_seqs, this->gene_names};
}

}  // namespace VJAssignment
