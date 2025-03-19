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
#ifndef DNA_SEQUENCE_UTILITIES_HEADER_H
#define DNA_SEQUENCE_UTILITIES_HEADER_H

// C++ headers
#include <cstdint>


namespace DNASequenceUtilities {

/// @brief Gets the reverse complement letter for an input DNA letter.
/// @param letter The input letter.
/// @return The output reverse complement letter.
char reverse_complement_letter(const char &letter);


/// @brief A hash function for DNA codons so that each is mapped to a unique
/// value.
/// @param codon The DNA codon.
/// @return An integer. All possible DNA codons will have a unique integer
/// assigned.
constexpr int64_t codon_hash(const char *codon);

/// @brief Converts an input DNA codon to an amino acid. Throws an exception
/// which can be handled by Python caller if unrecognized codon is supplied.
/// @param codon The DNA codon.
/// @return An amino acid letter.
char codon_to_aa(const char *codon);


}  // namespace DNASequenceUtilities

#endif
