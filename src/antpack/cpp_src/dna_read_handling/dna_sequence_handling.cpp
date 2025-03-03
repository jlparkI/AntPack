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
#include <unordered_set>
#include <string>
#include <stdexcept>
#include <vector>

// Library headers

// Project headers
#include "dna_sequence_handling.h"
#include "../utilities/dna_sequence_utilities.h"


namespace DNASequenceTools {


DNASeqTranslatorCpp::DNASeqTranslatorCpp(
        std::unordered_set<std::string> known_mab_kmers):
        known_mab_kmers(known_mab_kmers) {
}




std::string DNASeqTranslatorCpp::translate_dna_unknown_rf(std::string sequence,
              bool check_reverse_complement) {
    std::string output_string = "";
    int best_score = 0;

    for (size_t i=0; i < 3; i++) {
        std::string translation = this->translate_dna_known_rf(sequence, i,
                false);
        int sequence_length = translation.length();
        int translation_score = 0;

        for (size_t j=0; j < sequence_length - 8; j++) {
            std::string kmer = translation.substr(j, 9);
            if (this->known_mab_kmers.count(kmer) > 0)
                translation_score += 1;
        }

        if (translation_score > best_score) {
            best_score = translation_score;
            output_string = translation;
        }
    }
    if (check_reverse_complement) {
        for (size_t i=0; i < 3; i++) {
            std::string translation = this->translate_dna_known_rf(sequence, i,
                    true);
            int sequence_length = translation.length();
            int translation_score = 0;

            for (size_t j=0; j < sequence_length - 8; j++) {
                std::string kmer = translation.substr(j, 9);
                if (this->known_mab_kmers.count(kmer) > 0)
                    translation_score += 1;
            }

            if (translation_score > best_score) {
                best_score = translation_score;
                output_string = translation;
            }
        }
    }
    return output_string;
}



/// @brief Translates from DNA to protein when the reading frame
///        is known.
/// @param sequence The input DNA sequence. Should contain A, G, C, T
///        and (possibly) N only.
/// @param reading_frame How many positions from start to begin translating.
///        Should be one of 0, 1, 2.
/// @param reverse_complement If True, use the reverse complement of the
///        sequence rather than the sequence itself.
/// @return The DNA sequence translated to amino acids.
std::string DNASeqTranslatorCpp::translate_dna_known_rf(std::string sequence,
              int reading_frame,
              bool reverse_complement) {
    std::string rev_complement = "";
    std::string *target_string;
    std::vector<char> output_letters;
    std::string output_sequence = "";
    int sequence_length = sequence.length();

    output_letters.reserve(sequence.length());

    if (reading_frame < 0 || reading_frame > 2)
        throw std::runtime_error("Reading frame should be one of 0, 1 or 2.");

    if (reverse_complement) {
        std::vector<char> revcomp_letters;
        revcomp_letters.reserve(sequence.length());

        for (int i=sequence_length - 1; i >= 0; i--) {
            char revcomp_letter = DNASequenceUtilities::reverse_complement_letter(
                    sequence.at(i));
            revcomp_letters.push_back(revcomp_letter);
        }

        rev_complement = std::string(revcomp_letters.begin(),
                revcomp_letters.end());
        target_string = &rev_complement;
    } else {
        target_string = &sequence;
    }

    for (size_t i=reading_frame; i < sequence_length - 2; i+=3) {
        std::string codon = target_string->substr(i, 3);

        // Any codon containing an N is converted to X as
        // undefined.
        if (codon.at(0) == 'N' || codon.at(1) == 'N' ||
                codon.at(2) == 'N') {
            output_letters.push_back('X');
            continue;
        }

        char next_letter = DNASequenceUtilities::codon_to_aa(codon.c_str());
        output_letters.push_back(next_letter);
    }

    if (output_letters.size() > 0) {
        output_sequence = std::string(output_letters.begin(),
                output_letters.end());
    }
    return output_sequence;
}





}  // namespace DNASequenceTools
