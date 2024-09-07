#ifndef NUMBERING_UTILITIES_HEADER_H
#define NUMBERING_UTILITIES_HEADER_H

#include <string>
#include <vector>
#include "../numbering_constants.h"


void convert_numbering_to_imgt(std::vector<std::string> &og_numbering,
        std::vector<std::string> &con_numbering,
        std::string &scheme);


#endif
