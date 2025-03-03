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
#ifndef DNA_SEQUENCE_TRANSLATOR_HEADER_H
#define DNA_SEQUENCE_TRANSLATOR_HEADER_H

// C++ headers
#include <string>
#include <unordered_set>

// Library headers


// Project headers



namespace DNASequenceTools {


class DNASeqTranslatorCpp {
 public:
    explicit DNASeqTranslatorCpp(std::unordered_set<std::string>
            known_mab_kmers);

      /// @brief Translate from DNA to protein when the reading frame
      ///        is unknown and must be determined from the data.
      /// @param sequence The input DNA sequence. Should contain A, G, C,
      ///        T and (possibly) N only.
      /// @param check_reverse_complement If True, also check the reverse
      ///        complement.
      /// @return The DNA sequence translated into amino acids.
      std::string translate_dna_unknown_rf(std::string sequence,
              bool check_reverse_complement = false);

      /// @brief Translates from DNA to protein when the reading frame
      ///        is known.
      /// @param sequence The input DNA sequence. Should contain A, G, C, T
      ///        and (possibly) N only.
      /// @param reading_frame How many positions from start to begin translating.
      /// @param reverse_complement If True, use the reverse complement of the
      ///        sequence rather than the sequence itself.
      /// @return The DNA sequence translated to amino acids.
      std::string translate_dna_known_rf(std::string sequence,
              int reading_frame = 0,
              bool reverse_complement = false);

 private:
      std::unordered_set<std::string> known_mab_kmers;

};


}  // namespace DNASequenceTools



#endif
