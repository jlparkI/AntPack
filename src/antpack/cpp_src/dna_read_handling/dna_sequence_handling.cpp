// C++ headers
#include <unordered_set>
#include <string>
#include <stdexcept>
#include <vector>

// Library headers

// Project headers
#include "dna_sequence_handling.h"


namespace DNASequenceTools {

DNASeqTranslatorCpp::DNASeqTranslatorCpp(
        std::unordered_set<std::string> known_mab_kmers):
        known_mab_kmers(known_mab_kmers) {
}


std::string DNASeqTranslatorCpp::translate_dna_unknown_rf(std::string sequence,
              bool check_reverse_complement) {
}



std::string DNASeqTranslatorCpp::translate_dna_known_rf(std::string sequence,
              int reading_frame,
              bool reverse_complement) {
    std::string rev_complement = "";
    std::string *target_string;
    std::vector<char> output_letters;
    std::string output_sequence = "";
    int sequence_length = sequence.length();

    output_letters.reserve(sequence.length());

    if (reverse_complement) {
        std::vector<char> revcomp_letters;
        revcomp_letters.reserve(sequence.length());

        for (size_t i=sequence_length - 1; i >= 0; i--) {
            char letter = sequence.at(i);

            switch (letter) {
                case 'A':
                    revcomp_letters.push_back('T');
                    break;
                case 'G':
                    revcomp_letters.push_back('C');
                    break;
                case 'T':
                    revcomp_letters.push_back('A');
                    break;
                case 'C':
                    revcomp_letters.push_back('G');
                    break;
                case 'N':
                    revcomp_letters.push_back('N');
                    break;
                default:
                    // Exceptions thrown here are handled by the Python
                    // wrapper.
                    throw std::runtime_error(std::string("Invalid character "
                           "found in DNA sequence. Only accepted characters "
                                "are A, G, T, C and (possibly) N."));
                    break;
            }
            rev_complement = std::string(revcomp_letters.begin(),
                    revcomp_letters.end());
            target_string = &rev_complement;
        }
    } else {
        target_string = &sequence;
    }

    for (size_t i=0; i < sequence_length; i++) {
        char letter = target_string->at(i);

        switch (letter) {
            case 'A':
                output_letters.push_back('T');
                break;
            case 'G':
                output_letters.push_back('C');
                break;
            case 'T':
                output_letters.push_back('A');
                break;
            case 'C':
                output_letters.push_back('G');
                break;
            case 'N':
                output_letters.push_back('N');
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

    output_sequence = std::string(output_letters.begin(),
            output_letters.end());
    return output_sequence;
}


}  // namespace DNASequenceTools
