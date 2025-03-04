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

    std::tuple<std::string, std::string, double, double> assign_vj_genes(
            std::tuple<std::vector<std::string>,
            double, std::string, std::string> alignment,
            std::string sequence,
            std::string species, std::string mode);

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

    void assign_gene_by_identity(std::vector<std::string> &gene_seqs,
            std::vector<std::string> &gene_names,
            std::string &prepped_sequence,
            double &best_identity,
            std::string &best_gene_name,
            char gene_type);
    
    int assign_gene_by_evalue(std::vector<std::string> &gene_seqs,
            std::vector<std::string> &gene_names,
            int *encoded_sequence,
            double &best_identity,
            std::string &best_gene_name,
            char gene_type);

    int prep_sequence(std::string &prepped_sequence, std::string &sequence,
           std::tuple<std::vector<std::string>,
           double, std::string, std::string> &alignment);
};

}  // namespace VJAssignment

#endif  // SRC_ANTPACK_CPP_SRC_VJ_ASSIGNMENT_VJ_MATCH_COUNTER_H_
