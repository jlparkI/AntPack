/* Tools for handling consensus files, including numpy consensus files.
 * Much of the code in this file is adapted from Carl Rogers'
 * https://github.com/rogersce/cnpy library, which is very nice but
 * also contains functionality we do not need (writing npy files,
 * reading npy files with complex types etc.), so this has been
 * simplified somewhat.
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

// Library headers

// Project headers
#include "consensus_file_utilities.h"



/// @brief Reads a specially formatted text file into a
/// vector of vectors of strings. Each entry in the outer
/// vector is a vector of amino acids allowed at that
/// position. Note that an empty vector at a given position
/// indicates any AA is tolerated at that position. Therefore,
/// if a '-' is found at a given postiion, assume any AA is
/// tolerated there.
/// @param consFPath The location of the consensus file.
/// @param consensusAAs The list in which the results will
/// be stored.
/// @return An integer error code.
int cnpy::read_consensus_file(std::filesystem::path consFPath,
        std::vector<std::vector<std::string>> &consensusAAs) {
    if (!std::filesystem::exists(consFPath))
        return INVALID_CONSENSUS_FILE;

    std::ifstream file(consFPath.string());

    int lastPosition = 0;
    bool readNow = false;
    std::string current_line;

    while (std::getline(file, current_line)) {
        if (current_line.at(0) == '#') {
            readNow = true;
            continue;
        }
        if (current_line.at(0) == '/')
            break;
        if (!readNow)
            continue;

        std::stringstream split_string(current_line);
        std::string segment;
        bool firstSegment = true, allAAsAllowed = false;
        std::vector<std::string> allowedAAs;

        while (std::getline(split_string, segment, ',')) {
            if (firstSegment) {
                int position = std::stoi(segment);
                if (position - 1 != lastPosition)
                    return INVALID_CONSENSUS_FILE;
                firstSegment = false;
                lastPosition += 1;
            } else {
                if (segment == "-") {
                    consensusAAs.push_back(std::vector<std::string>{});
                    allAAsAllowed = true;
                    break;
                }
                allowedAAs.push_back(segment);
            }
        }

        if (!allAAsAllowed)
            consensusAAs.push_back(allowedAAs);
    }

    return VALID_CONSENSUS_FILE;
}





/// @brief Reads a specially formatted text file into a
/// vector of strings. Each entry is a V-gene or J-gene,
/// which must be of a prespecified length in order to
/// be accepted.
/// @param filepath The location of the file.
/// @param gene_list The vector in which the results will
/// be stored.
/// @param expected_length Enforces that all loaded genes will
/// be of an expected length.
/// @return An integer error code.
int cnpy::read_tcr_vj_gene_file(std::filesystem::path filepath,
        std::vector<std::string> &gene_list,
        size_t expected_length) {
    if (!std::filesystem::exists(filepath))
        return INVALID_CONSENSUS_FILE;

    std::ifstream file(filepath.string());

    int lastPosition = 0;
    std::string current_line;

    while (std::getline(file, current_line)) {
        if (current_line.length() < 2)
            break;

        if (current_line.length() != expected_length)
            return INVALID_CONSENSUS_FILE;

        gene_list.push_back(current_line);
    }

    return VALID_CONSENSUS_FILE;
}




// Test for big-endianness.
char cnpy::BigEndianTest() {
    int x = 1;
    return (((char *)&x)[0]) ? '<' : '>';
}


char cnpy::map_type(const std::type_info& t)
{
    if(t == typeid(float) ) return 'f';
    if(t == typeid(double) ) return 'f';
    if(t == typeid(long double) ) return 'f';

    if(t == typeid(int) ) return 'i';
    if(t == typeid(char) ) return 'i';
    if(t == typeid(short) ) return 'i';
    if(t == typeid(long) ) return 'i';
    if(t == typeid(long long) ) return 'i';

    if(t == typeid(unsigned char) ) return 'u';
    if(t == typeid(unsigned short) ) return 'u';
    if(t == typeid(unsigned long) ) return 'u';
    if(t == typeid(unsigned long long) ) return 'u';
    if(t == typeid(unsigned int) ) return 'u';

    if(t == typeid(bool) ) return 'b';

    // Note that we do not check for complex types here. However,
    // we are loading npy files we have generated, which will
    // NEVER be complex-valued. If we have saved a complex-valued
    // npy file there has been an error and '?' is the appropriate
    // thing to return.

    else return '?';
}


void cnpy::parse_npy_header(unsigned char* buffer,size_t& word_size, std::vector<size_t>& shape, bool& fortran_order) {
    //std::string magic_string(buffer,6);
    uint16_t header_len = *reinterpret_cast<uint16_t*>(buffer+8);
    std::string header(reinterpret_cast<char*>(buffer+9),header_len);

    size_t loc1, loc2;

    //fortran order
    loc1 = header.find("fortran_order")+16;
    fortran_order = (header.substr(loc1,4) == "True" ? true : false);

    //shape
    loc1 = header.find("(");
    loc2 = header.find(")");

    std::regex num_regex("[0-9][0-9]*");
    std::smatch sm;
    shape.clear();

    std::string str_shape = header.substr(loc1+1,loc2-loc1-1);
    while(std::regex_search(str_shape, sm, num_regex)) {
        shape.push_back(std::stoi(sm[0].str()));
        str_shape = sm.suffix().str();
    }

    //endian, word size, data type
    //byte order code | stands for not applicable. 
    //not sure when this applies except for byte array
    loc1 = header.find("descr")+9;
    bool littleEndian = (header[loc1] == '<' || header[loc1] == '|' ? true : false);
    if (!littleEndian){
        throw std::runtime_error(std::string("There is a serious error with library installation; "
                    "please report."));
    }

    std::string str_ws = header.substr(loc1+2);
    loc2 = str_ws.find("'");
    word_size = atoi(str_ws.substr(0,loc2).c_str());
}


void cnpy::parse_npy_header(FILE* fp, size_t& word_size, std::vector<size_t>& shape, bool& fortran_order) {  
    char buffer[256];
    size_t res = fread(buffer,sizeof(char),11,fp);       
    if(res != 11)
        throw std::runtime_error("parse_npy_header: failed fread");
    std::string header = fgets(buffer,256,fp);
    if (!(header[header.size()-1] == '\n')){
        throw std::runtime_error(std::string("There is a serious error with library installation; "
                    "please report."));
    }

    size_t loc1, loc2;

    //fortran order
    loc1 = header.find("fortran_order");
    if (loc1 == std::string::npos)
        throw std::runtime_error("parse_npy_header: failed to find header keyword: 'fortran_order'");
    loc1 += 16;
    fortran_order = (header.substr(loc1,4) == "True" ? true : false);

    //shape
    loc1 = header.find("(");
    loc2 = header.find(")");
    if (loc1 == std::string::npos || loc2 == std::string::npos)
        throw std::runtime_error("parse_npy_header: failed to find header keyword: '(' or ')'");

    std::regex num_regex("[0-9][0-9]*");
    std::smatch sm;
    shape.clear();

    std::string str_shape = header.substr(loc1+1,loc2-loc1-1);
    while(std::regex_search(str_shape, sm, num_regex)) {
        shape.push_back(std::stoi(sm[0].str()));
        str_shape = sm.suffix().str();
    }

    //endian, word size, data type
    //byte order code | stands for not applicable. 
    //not sure when this applies except for byte array
    loc1 = header.find("descr");
    if (loc1 == std::string::npos)
        throw std::runtime_error("parse_npy_header: failed to find header keyword: 'descr'");
    loc1 += 9;
    bool littleEndian = (header[loc1] == '<' || header[loc1] == '|' ? true : false);
    if (!littleEndian){
        throw std::runtime_error(std::string("There is a serious error with library installation; "
                    "please report."));
    }

    std::string str_ws = header.substr(loc1+2);
    loc2 = str_ws.find("'");
    word_size = atoi(str_ws.substr(0,loc2).c_str());
}

cnpy::NpyArray load_the_npy_file(FILE* fp) {
    std::vector<size_t> shape;
    size_t word_size;
    bool fortran_order;
    cnpy::parse_npy_header(fp,word_size,shape,fortran_order);

    cnpy::NpyArray arr(shape, word_size, fortran_order);
    size_t nread = fread(arr.data<char>(),1,arr.num_bytes(),fp);
    if(nread != arr.num_bytes())
        throw std::runtime_error("load_the_npy_file: failed fread");
    return arr;
}


cnpy::NpyArray cnpy::npy_load(std::string fname) {

    FILE* fp = fopen(fname.c_str(), "rb");

    if(!fp) throw std::runtime_error("npy_load: Unable to open file "+fname);

    NpyArray arr = load_the_npy_file(fp);

    fclose(fp);
    return arr;
}
