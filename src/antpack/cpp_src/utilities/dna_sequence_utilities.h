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
