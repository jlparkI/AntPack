#ifndef UTILITIES_HEADER_H
#define UTILITIES_HEADER_H

#include <string>
#include <vector>
#include <tuple>
#include <map>
#include <set>
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <unordered_map>
#include "numbering_constants.h"

int validate_sequence(std::string query_sequence);

int convert_sequence_to_array(int *queryAsIdx, std::string query_sequence);

int sort_position_codes_utility(std::vector<std::string> position_codes,
        std::string scheme, std::vector<std::string> &orderedTranslatedCodes);

int read_consensus_file(std::filesystem::path consFPath,
        std::vector<std::vector<std::string>> &allowedAAs);

int build_msa_utility(std::vector<std::string> sequences,
        std::vector<std::tuple<std::vector<std::string>, double, std::string, std::string>> annotations,
        std::vector<std::string> &positionCodes,
        std::vector<std::string> &alignedSeqs,
        std::string scheme);

#endif
