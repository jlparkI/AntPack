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
#ifndef SRC_ANTPACK_CPP_SRC_VJ_ASSIGNMENT_VJ_MATCH_COUNTER_H_
#define SRC_ANTPACK_CPP_SRC_VJ_ASSIGNMENT_VJ_MATCH_COUNTER_H_

// C++ headers
#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <vector>
#include <string>
#include <tuple>
#include <map>
#include <unordered_set>
#include <memory>

// Library headers


// Project headers
#include "../annotator_classes/ig_aligner.h"
#include "../utilities/utilities.h"


namespace nb = nanobind;



namespace VJAssignment {

// The input sequence must have this length in order to be
// considered. Caller should already have extracted the IMGT
// standard positions which are all that is used for VJ
// assignment.
static const int REQUIRED_SEQUENCE_LENGTH = 128;
static const int VALID_SEQUENCE = 1;
static const int INVALID_SEQUENCE = 0;

class VJMatchCounter {
 public:
    VJMatchCounter(
            std::map<std::string, std::vector<std::string>> gene_names,
            std::map<std::string, std::vector<std::string>> gene_seqs,
            nb::ndarray<double, nb::shape<22,22>, nb::device::cpu, nb::c_contig> blosum_matrix,
            std::string scheme,
            std::string consensus_filepath);

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
        assign_vj_genes(std::tuple<std::vector<std::string>,
            double, std::string, std::string> alignment,
            std::string sequence,
            std::string species, std::string mode);

    /// @brief Gets the sequence of a named V or J gene for an indicated
    /// species.
    /// @param query_name The name of the gene.
    /// @param species The name of the species.
    /// @return The sequence in gapped IMGT format. If the gene is not found,
    /// this string is empty.
    std::string get_vj_gene_sequence(std::string query_name,
            std::string species);

    std::tuple<std::map<std::string, std::vector<std::string>>,
        std::map<std::string, std::vector<std::string>>> get_seq_lists();

 protected:
    std::map<std::string, std::vector<std::string>> gene_seqs;
    std::map<std::string, std::vector<std::string>> gene_names;
    nb::ndarray<double, nb::shape<22, 22>, nb::device::cpu,
        nb::c_contig> blosum_matrix;
    std::string scheme;

    std::map<std::string, std::map<std::string, int>> names_to_positions;
    std::unordered_set<std::string> essential_imgt_map;

    std::unique_ptr<NumberingTools::IGAligner> h_aligner;
    std::unique_ptr<NumberingTools::IGAligner> k_aligner;
    std::unique_ptr<NumberingTools::IGAligner> l_aligner;

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
    internal_vj_assignment(std::tuple<std::vector<std::string>,
            double, std::string, std::string> &alignment,
            std::string &sequence,
            const std::string &species,
            const std::string &mode);

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
    void assign_gene_by_identity(std::vector<std::string> &gene_seqs,
            std::vector<std::string> &gene_names,
            std::string &prepped_sequence,
            double &best_identity,
            std::string &best_gene_name,
            char gene_type);
    
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
    int assign_gene_by_evalue(std::vector<std::string> &gene_seqs,
            std::vector<std::string> &gene_names,
            int *encoded_sequence,
            double &best_identity,
            std::string &best_gene_name,
            char gene_type);

    /// @brief Converts a sequence to a standard format by extracting the
    /// 128 standard IMGT positions for easy comparison to pre-gapped
    /// V and J genes.
    /// @param prepped_sequence The output is stored here.
    /// @param sequence The input sequence.
    /// @param alignment A tuple of numbering, percent identity, chain
    /// and error message such as is returned by the Annotator classes.
    /// @return An error code indicating success or failure.
    int prep_sequence(std::string &prepped_sequence, std::string &sequence,
           std::tuple<std::vector<std::string>,
           double, std::string, std::string> &alignment);
};

}  // namespace VJAssignment

#endif  // SRC_ANTPACK_CPP_SRC_VJ_ASSIGNMENT_VJ_MATCH_COUNTER_H_
