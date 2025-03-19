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
#ifndef CONSENSUS_FILE_UTILITIES_HEADER_H
#define CONSENSUS_FILE_UTILITIES_HEADER_H

// C++ headers
#include <vector>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <regex>

// Library headers


// Project headers




namespace cnpy {


static constexpr int VALID_CONSENSUS_FILE = 1;
static constexpr int INVALID_CONSENSUS_FILE = 0;


struct NpyArray {
    NpyArray(const std::vector<size_t>& _shape, size_t _word_size, bool _fortran_order) :
        shape(_shape), word_size(_word_size), fortran_order(_fortran_order)
    {
        num_vals = 1;
        for(size_t i = 0;i < shape.size();i++) num_vals *= shape[i];
        data_holder = std::shared_ptr<std::vector<char>>(
            new std::vector<char>(num_vals * word_size));
    }

    NpyArray() : shape(0), word_size(0), fortran_order(0), num_vals(0) { }

    template<typename T>
    T* data() {
        return reinterpret_cast<T*>(&(*data_holder)[0]);
    }

    template<typename T>
    const T* data() const {
        return reinterpret_cast<T*>(&(*data_holder)[0]);
    }

    template<typename T>
    std::vector<T> as_vec() const {
        const T* p = data<T>();
        return std::vector<T>(p, p+num_vals);
    }

    size_t num_bytes() const {
        return data_holder->size();
    }

    std::shared_ptr<std::vector<char>> data_holder;
    std::vector<size_t> shape;
    size_t word_size;
    bool fortran_order;
    size_t num_vals;
};

int read_consensus_file(std::filesystem::path consFPath,
        std::vector<std::vector<std::string>> &allowedAAs);

int read_tcr_vj_gene_file(std::filesystem::path filepath,
        std::vector<std::string> &gene_list,
        size_t expected_length);

char BigEndianTest();

char map_type(const std::type_info& t);

void parse_npy_header(FILE* fp, size_t& word_size,
        std::vector<size_t>& shape, bool& fortran_order);

void parse_npy_header(unsigned char* buffer, size_t& word_size,
        std::vector<size_t>& shape, bool& fortran_order);

NpyArray npy_load(std::string fname);



}  // namespace cnpy

#endif
