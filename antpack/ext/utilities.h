#ifndef UTILITIES_HEADER_H
#define UTILITIES_HEADER_H

#include <string>
#include <vector>
#include <tuple>
#include <map>
#include <set>
#include <algorithm>
#include <cmath>
#include "numbering_constants.h"

int validate_sequence(std::string query_sequence);
int convert_sequence_to_array(int *queryAsIdx, std::string query_sequence);
std::tuple<std::vector<std::string>, int> sort_position_codes_cpp(std::vector<std::string> position_codes,
        std::string scheme);

#endif
