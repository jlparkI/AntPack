#ifndef UTILITIES_HEADER_H
#define UTILITIES_HEADER_H

#include <string>

int validate_sequence(std::string query_sequence);
int convert_sequence_to_array(int *queryAsIdx, std::string query_sequence);

#endif