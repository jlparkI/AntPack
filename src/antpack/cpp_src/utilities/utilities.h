#ifndef UTILITIES_HEADER_H
#define UTILITIES_HEADER_H

#include <string>
#include <vector>
#include <tuple>
#include <map>
#include <set>
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <unordered_map>
#include "../numbering_constants.h"

int validate_sequence(std::string &query_sequence);
int validate_x_sequence(std::string &query_sequence);
int validate_gapped_sequence(std::string &query_sequence);

int convert_sequence_to_array(int *queryAsIdx, std::string &query_sequence);
int convert_x_sequence_to_array(int *queryAsIdx, std::string &query_sequence);

int sort_position_codes_utility(std::vector<std::string> &position_codes,
        std::string scheme, std::vector<std::string> &ordered_translated_codes);

int build_msa_utility(std::vector<std::string> &sequences,
        std::vector<std::tuple<std::vector<std::string>, double, std::string, std::string>> &annotations,
        std::vector<std::string> &position_codes,
        std::vector<std::string> &aligned_seqs,
        const std::string &scheme);

int trim_alignment_utility(const std::string &sequence,
        std::tuple<std::vector<std::string>, double, std::string, std::string> &alignment,
        std::vector<std::string> &trimmed_alignment, int &exstart, int &exend,
        std::vector<char> &trimmed_seq);

#endif
