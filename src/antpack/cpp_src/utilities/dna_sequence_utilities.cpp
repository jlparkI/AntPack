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
// C++ headers
#include <string>
#include <stdexcept>

// Library headers

// Project headers
#include "dna_sequence_utilities.h"



namespace DNASequenceUtilities {


/// @brief Gets the reverse complement letter for an input DNA letter.
/// @param letter The input letter.
/// @return The output reverse complement letter.
char reverse_complement_letter(const char &letter) {
    switch (letter) {
        case 'A':
            return 'T';
            break;
        case 'G':
            return 'C';
            break;
        case 'T':
            return 'A';
            break;
        case 'C':
            return 'G';
            break;
        case 'N':
            return 'N';
            break;
        default:
            // Exceptions thrown here are handled by the Python
            // wrapper.
            throw std::runtime_error(std::string("Invalid character "
                   "found in DNA sequence. Only accepted characters "
                        "are A, G, T, C and (possibly) N."));
            break;
    }
}



/// @brief A hash function for DNA codons so that each is mapped to a unique
/// value.
/// @param codon The DNA codon.
/// @return An integer. All possible DNA codons will have a unique integer
/// assigned.
constexpr int64_t codon_hash(const char *codon) {
    const int64_t p = 131;
    int64_t total = 0;
    int64_t multiplier = 1;
    for (size_t i=0; codon[i] != '\0'; i++) {
        total += multiplier * codon[i];
        multiplier *= p;
    }
    return total;
}



/// @brief Converts an input DNA codon to an amino acid. Throws an exception
/// which can be handled by Python caller if unrecognized codon is supplied.
/// @param codon The DNA codon.
/// @return An amino acid letter.
char codon_to_aa(const char *codon) {
    switch (codon_hash(codon)) {
        case codon_hash("TTT"):
        case codon_hash("TTC"):
            return 'F';
            break;
        case codon_hash("TTA"):
        case codon_hash("TTG"):
        case codon_hash("CTT"):
        case codon_hash("CTC"):
        case codon_hash("CTA"):
        case codon_hash("CTG"):
            return 'L';
            break;
        case codon_hash("ATT"):
        case codon_hash("ATC"):
        case codon_hash("ATA"):
            return 'I';
            break;
        case codon_hash("ATG"):
            return 'M';
            break;
        case codon_hash("GTT"):
        case codon_hash("GTC"):
        case codon_hash("GTA"):
        case codon_hash("GTG"):
            return 'V';
            break;
        case codon_hash("TCT"):
        case codon_hash("TCC"):
        case codon_hash("TCA"):
        case codon_hash("TCG"):
            return 'S';
            break;
        case codon_hash("CCT"):
        case codon_hash("CCC"):
        case codon_hash("CCA"):
        case codon_hash("CCG"):
            return 'P';
            break;
        case codon_hash("ACT"):
        case codon_hash("ACC"):
        case codon_hash("ACA"):
        case codon_hash("ACG"):
            return 'T';
            break;
        case codon_hash("GCT"):
        case codon_hash("GCC"):
        case codon_hash("GCA"):
        case codon_hash("GCG"):
            return 'A';
            break;
        case codon_hash("TAT"):
        case codon_hash("TAC"):
            return 'Y';
            break;
        case codon_hash("TAA"):
        case codon_hash("TAG"):
        case codon_hash("TGA"):
            return '*';
            break;
        case codon_hash("CAT"):
        case codon_hash("CAC"):
            return 'H';
            break;
        case codon_hash("CAA"):
        case codon_hash("CAG"):
            return 'Q';
            break;
        case codon_hash("AAT"):
        case codon_hash("AAC"):
            return 'N';
            break;
        case codon_hash("AAA"):
        case codon_hash("AAG"):
            return 'K';
            break;
        case codon_hash("GAT"):
        case codon_hash("GAC"):
            return 'D';
            break;
        case codon_hash("GAA"):
        case codon_hash("GAG"):
            return 'E';
            break;
        case codon_hash("TGT"):
        case codon_hash("TGC"):
            return 'C';
            break;
        case codon_hash("TGG"):
            return 'W';
            break;
        case codon_hash("CGT"):
        case codon_hash("CGC"):
        case codon_hash("CGA"):
        case codon_hash("CGG"):
            return 'R';
            break;
        case codon_hash("AGT"):
        case codon_hash("AGC"):
            return 'S';
            break;
        case codon_hash("AGA"):
        case codon_hash("AGG"):
            return 'R';
            break;
        case codon_hash("GGT"):
        case codon_hash("GGC"):
        case codon_hash("GGA"):
        case codon_hash("GGG"):
            return 'G';
            break;
        default:
            // Exceptions thrown here are handled by the Python
            // wrapper.
            throw std::runtime_error(std::string("Invalid character "
                   "found in DNA sequence. Only accepted characters "
                        "are A, G, T, C and (possibly) N."));
            break;
    }
}


}  // namespace DNASequenceUtilities
