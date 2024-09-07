#ifndef VJ_MATCH_COUNTER_HEADER_H
#define VJ_MATCH_COUNTER_HEADER_H

#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <vector>
#include <string>
#include <tuple>
#include <map>
#include <iostream>
#include "../utilities/utilities.h"
#include "../utilities/numbering_utilities.h"


namespace py = pybind11;

// The input sequence must have this length in order to be
// considered. Caller should already have extracted the IMGT
// standard positions which are all that is used for VJ
// assignment.
#define REQUIRED_SEQUENCE_LENGTH 128
#define VALID_SEQUENCE 1
#define INVALID_SEQUENCE 0




class VJMatchCounter {
    public:
        VJMatchCounter(std::map<std::string, std::vector<std::string>> gene_names,
                std::map<std::string, std::vector<std::string>> gene_seqs,
                py::array_t<int16_t, py::array::c_style> blosum_matrix,
                std::string scheme);

        std::tuple<std::string, std::string, double, double> assign_vj_genes(std::tuple<std::vector<std::string>,
                double, std::string, std::string> alignment, std::string sequence,
                std::string species, std::string mode);

        std::string get_vj_gene_sequence(std::string query_name, std::string species);

        std::tuple<std::map<std::string, std::vector<std::string>>,
                std::map<std::string, std::vector<std::string>>>  get_seq_lists();

    protected:
        std::map<std::string, std::vector<std::string>> gene_seqs;
        std::map<std::string, std::vector<std::string>> gene_names;
        py::array_t<int16_t, py::array::c_style> blosum_matrix;
        std::string scheme;

        std::map<std::string, std::map<std::string, int>> names_to_positions;

        std::unordered_set<std::string> essential_imgt_map;

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

        void prep_sequence(std::string &prepped_sequence, std::string &sequence,
                std::vector<std::string> &numbering);

};

#endif
