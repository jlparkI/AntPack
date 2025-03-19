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
#include <set>
#include <string>
#include <map>

// Project headers
#include "utilities.h"


namespace SequenceUtilities {


// Checks a query sequence to ensure it contains only recognized
// amino acids.
int validate_sequence(std::string &query_sequence) {
    bool valid_query = true;

    if (query_sequence.length() == 0)
        return INVALID_SEQUENCE;

    for (char & c : query_sequence) {
        // Switch will provide a (small) speed gain over map
        // since compiler should convert it to a lookup table
        switch (c) {
            case 'A':
            case 'C':
            case 'D':
            case 'E':
            case 'F':
            case 'G':
            case 'H':
            case 'I':
            case 'K':
            case 'L':
            case 'M':
            case 'N':
            case 'P':
            case 'Q':
            case 'R':
            case 'S':
            case 'T':
            case 'V':
            case 'W':
            case 'Y':
                break;
            default:
                valid_query = false;
                break;
        }
        if (!valid_query)
            break;
    }
    if (!valid_query)
        return INVALID_SEQUENCE;

    return VALID_SEQUENCE;
}


// Checks a query sequence to ensure it contains only recognized
// amino acids OR the letter x. Use this in preference to
// validate_sequence when 'X' is allowed or ok.
int validate_x_sequence(std::string &query_sequence) {
    bool valid_query = true;

    if (query_sequence.length() == 0)
        return INVALID_SEQUENCE;

    for (char & c : query_sequence) {
        // Switch will provide a (small) speed gain over map
        // since compiler should convert it to a lookup table
        switch (c) {
            case 'A':
            case 'C':
            case 'D':
            case 'E':
            case 'F':
            case 'G':
            case 'H':
            case 'I':
            case 'K':
            case 'L':
            case 'M':
            case 'N':
            case 'P':
            case 'Q':
            case 'R':
            case 'S':
            case 'T':
            case 'V':
            case 'W':
            case 'Y':
            case 'X':
                break;
            default:
                valid_query = false;
                break;
        }
        if (!valid_query)
            break;
    }
    if (!valid_query)
        return INVALID_SEQUENCE;

    return VALID_SEQUENCE;
}


// Checks a query sequence to ensure it contains only recognized
// amino acids OR GAPS. Use this in preference to validate_sequence
// when gaps are allowed / ok.
int validate_gapped_sequence(std::string &query_sequence) {
    bool valid_query = true;

    if (query_sequence.length() == 0)
        return INVALID_SEQUENCE;

    for (char & c : query_sequence) {
        // Switch will provide a (small) speed gain over map
        // since compiler should convert it to a lookup table
        switch (c) {
            case 'A':
            case 'C':
            case 'D':
            case 'E':
            case 'F':
            case 'G':
            case 'H':
            case 'I':
            case 'K':
            case 'L':
            case 'M':
            case 'N':
            case 'P':
            case 'Q':
            case 'R':
            case 'S':
            case 'T':
            case 'V':
            case 'W':
            case 'Y':
            case '-':
                break;
            default:
                valid_query = false;
                break;
        }
        if (!valid_query)
            break;
    }
    if (!valid_query)
        return INVALID_SEQUENCE;

    return VALID_SEQUENCE;
}




// Converts an input sequence to a numeric representation as an
// array of integers. This function does not do any safety checks
// on the input so it should never be accessed by Python functions --
// only ever by C++ functions that DO perform safety checks on
// input.
int convert_sequence_to_array(int *queryAsIdx, std::string &query_sequence) {    
    // Translate the query sequence into an integer encoding. This
    // is more verbose than using std::map but should be slightly faster
    // since compiler will convert to lookup table. If unexpected characters
    // are encountered, abort and return an error code.
    for (size_t i=0; i < query_sequence.length(); i++) {
        switch (query_sequence[i]) {
            case 'A':
                queryAsIdx[i] = 0;
                break;
            case 'C':
                queryAsIdx[i] = 1;
                break;
            case 'D':
                queryAsIdx[i] = 2;
                break;
            case 'E':
                queryAsIdx[i] = 3;
                break;
            case 'F':
                queryAsIdx[i] = 4;
                break;
            case 'G':
                queryAsIdx[i] = 5;
                break;
            case 'H':
                queryAsIdx[i] = 6;
                break;
            case 'I':
                queryAsIdx[i] = 7;
                break;
            case 'K':
                queryAsIdx[i] = 8;
                break;
            case 'L':
                queryAsIdx[i] = 9;
                break;
            case 'M':
                queryAsIdx[i] = 10;
                break;
            case 'N':
                queryAsIdx[i] = 11;
                break;
            case 'P':
                queryAsIdx[i] = 12;
                break;
            case 'Q':
                queryAsIdx[i] = 13;
                break;
            case 'R':
                queryAsIdx[i] = 14;
                break;
            case 'S':
                queryAsIdx[i] = 15;
                break;
            case 'T':
                queryAsIdx[i] = 16;
                break;
            case 'V':
                queryAsIdx[i] = 17;
                break;
            case 'W':
                queryAsIdx[i] = 18;
                break;
            case 'Y':
                queryAsIdx[i] = 19;
                break;

            default:
                return INVALID_SEQUENCE;
                break;
        }
    }
    return VALID_SEQUENCE;
}


// Converts an input sequence to a numeric representation as an
// array of integers. This function does not do any safety checks
// on the input so it should never be accessed by Python functions --
// only ever by C++ functions that DO perform safety checks on
// input. This function unlike convert_sequence_to_array
// allows the letter X to be present in the input.
int convert_x_sequence_to_array(int *queryAsIdx, std::string &query_sequence) {  
    // Translate the query sequence into an integer encoding. This
    // is more verbose than using std::map but should be slightly faster
    // since compiler will convert to lookup table. If unexpected characters
    // are encountered, abort and return an error code.
    for (size_t i=0; i < query_sequence.length(); i++) {
        switch (query_sequence[i]) {
            case 'A':
                queryAsIdx[i] = 0;
                break;
            case 'C':
                queryAsIdx[i] = 1;
                break;
            case 'D':
                queryAsIdx[i] = 2;
                break;
            case 'E':
                queryAsIdx[i] = 3;
                break;
            case 'F':
                queryAsIdx[i] = 4;
                break;
            case 'G':
                queryAsIdx[i] = 5;
                break;
            case 'H':
                queryAsIdx[i] = 6;
                break;
            case 'I':
                queryAsIdx[i] = 7;
                break;
            case 'K':
                queryAsIdx[i] = 8;
                break;
            case 'L':
                queryAsIdx[i] = 9;
                break;
            case 'M':
                queryAsIdx[i] = 10;
                break;
            case 'N':
                queryAsIdx[i] = 11;
                break;
            case 'P':
                queryAsIdx[i] = 12;
                break;
            case 'Q':
                queryAsIdx[i] = 13;
                break;
            case 'R':
                queryAsIdx[i] = 14;
                break;
            case 'S':
                queryAsIdx[i] = 15;
                break;
            case 'T':
                queryAsIdx[i] = 16;
                break;
            case 'V':
                queryAsIdx[i] = 17;
                break;
            case 'W':
                queryAsIdx[i] = 18;
                break;
            case 'Y':
                queryAsIdx[i] = 19;
                break;
            case 'X':
                queryAsIdx[i] = 20;
                break;

            default:
                return INVALID_SEQUENCE;
                break;
        }
    }
    return VALID_SEQUENCE;
}




// Converts an input sequence to a numeric representation as an
// array of integers. This function does not do any safety checks
// on the input so it should never be accessed by Python functions --
// only ever by C++ functions that DO perform safety checks on
// input. This function unlike convert_sequence_to_array
// allows EITHER the letter X OR gaps to be present in the input.
int convert_gapped_x_sequence_to_array(int *queryAsIdx, std::string &query_sequence) { 
    // Translate the query sequence into an integer encoding. This
    // is more verbose than using std::map but should be slightly faster
    // since compiler will convert to lookup table. If unexpected characters
    // are encountered, abort and return an error code.
    for (size_t i=0; i < query_sequence.length(); i++) {
        switch (query_sequence[i]) {
            case 'A':
                queryAsIdx[i] = 0;
                break;
            case 'C':
                queryAsIdx[i] = 1;
                break;
            case 'D':
                queryAsIdx[i] = 2;
                break;
            case 'E':
                queryAsIdx[i] = 3;
                break;
            case 'F':
                queryAsIdx[i] = 4;
                break;
            case 'G':
                queryAsIdx[i] = 5;
                break;
            case 'H':
                queryAsIdx[i] = 6;
                break;
            case 'I':
                queryAsIdx[i] = 7;
                break;
            case 'K':
                queryAsIdx[i] = 8;
                break;
            case 'L':
                queryAsIdx[i] = 9;
                break;
            case 'M':
                queryAsIdx[i] = 10;
                break;
            case 'N':
                queryAsIdx[i] = 11;
                break;
            case 'P':
                queryAsIdx[i] = 12;
                break;
            case 'Q':
                queryAsIdx[i] = 13;
                break;
            case 'R':
                queryAsIdx[i] = 14;
                break;
            case 'S':
                queryAsIdx[i] = 15;
                break;
            case 'T':
                queryAsIdx[i] = 16;
                break;
            case 'V':
                queryAsIdx[i] = 17;
                break;
            case 'W':
                queryAsIdx[i] = 18;
                break;
            case 'Y':
                queryAsIdx[i] = 19;
                break;
            case 'X':
                queryAsIdx[i] = 20;
                break;
            case '-':
                queryAsIdx[i] = 21;
                break;

            default:
                return INVALID_SEQUENCE;
                break;
        }
    }
    return VALID_SEQUENCE;
}








// Takes as input a list of position codes (e.g. IMGT) and sorts it in a way
// that respects the properties of the scheme. If any invalid position
// codes are encountered, the function will stop and return an invalid
// error code. Note that since this function constructs a map on
// each call, calling it 1000s of times is not efficient or advisable.
// Indeed, since it should only be necessary to call it once per
// dataset, that should never actually happen.
int sort_position_codes_utility(std::vector<std::string> &position_codes,
    std::string scheme, std::vector<std::string> &ordered_translated_codes) {
    std::vector<float> ordered_float_codes;
    float arbitraryDivisor = 1000.0;
    int max_position = 128;
    if (scheme == "aho")
        max_position = 149;

    // Do not attempt to sort if no codes were supplied; simply
    // return what we were given.
    if (position_codes.size() == 0)
        return VALID_SEQUENCE;

    // Map from a standard alphabet Our preference would be to
    // number insertions as _1, _2 etc, but most numbering programs
    // use letters, so we do the same here for consistency and ease
    // of comparison. This does limit the number of possible insertions
    // BUT there should never be a case where the user actually needs more than
    // 78 insertion codes (if they really do there's something wrong with the
    // alignment and/or it's not really an antibody...)
    const std::map<std::string, int> alphabet = {{"A", 1}, {"B", 2}, {"C", 3}, {"D", 4}, {"E", 5},
                                {"F", 6}, {"G", 7}, {"H", 8}, {"I", 9}, {"J", 10}, {"K", 11}, {"L", 12},
                                {"M", 13}, {"N", 14}, {"O", 15}, {"P", 16}, {"Q", 17}, {"R", 18}, {"S", 19},
                                {"T", 20}, {"U", 21}, {"V", 22}, {"W", 23}, {"X", 24}, {"Y", 25}, {"Z", 26},
                                {"AA", 27}, {"BB", 28}, {"CC", 29}, {"DD", 30}, {"EE", 31},
                                {"FF", 32}, {"GG", 33}, {"HH", 34}, {"II", 35}, {"JJ", 36}, {"KK", 37}, {"LL", 38},
                                {"MM", 39}, {"NN", 40}, {"OO", 41}, {"PP", 42}, {"QQ", 43}, {"RR", 44}, {"SS", 45},
                                {"TT", 46}, {"UU", 47}, {"VV", 48}, {"WW", 49}, {"XX", 50}, {"YY", 51}, {"ZZ", 52},
                                {"AAA", 53}, {"BBB", 54}, {"CCC", 55}, {"DDD", 56}, {"EEE", 57},
                                {"FFF", 58}, {"GGG", 59}, {"HHH", 60}, {"III", 61}, {"JJJ", 62}, {"KKK", 63}, {"LLL", 64},
                                {"MMM", 65}, {"NNN", 66}, {"OOO", 67}, {"PPP", 68}, {"QQQ", 69}, {"RRR", 70}, {"SSS", 71},
                                {"TTT", 72}, {"UUU", 73}, {"VVV", 74}, {"WWW", 75}, {"XXX", 76}, {"YYY", 77}, {"ZZZ", 78}
                                };
    std::map<int, std::string> reverse_alphabet;
    for ( const auto & p : alphabet )
        reverse_alphabet[p.second] = p.first;


    // The IMGT scheme has several positions at which we have to
    // "count backwards", which is an annoying quirk of
    // that scheme. We have to hard-code these unfortunately.
    std::set<int> backwards_count_positions;
    if (scheme == "imgt")
        backwards_count_positions = {33, 61, 112};


    for (size_t i=0; i < position_codes.size(); i++) {
        std::string input_code = position_codes[i];
        std::string non_numeric_portion;
        int numeric_portion;
        float converted_code;

        // This will throw if the string starts with a letter but
        // will otherwise extract the integer piece,
        // which rules out certain kinds of invalid codes the user
        // might pass in error (e.g. '-', 'A128').
        try {
            numeric_portion = std::stoi(input_code);
        }
        catch (...) {
            return INVALID_SEQUENCE;
        }
        if (numeric_portion <= 0 || numeric_portion > max_position)
            return INVALID_SEQUENCE;

        converted_code = static_cast<float>(numeric_portion);

        for (char & c : input_code) {
            if (!std::isdigit(c))
                non_numeric_portion += c;
        }

        if (non_numeric_portion.length() == 0) {
            ordered_float_codes.push_back(converted_code);
            continue;
        }

        if (alphabet.find(non_numeric_portion) != alphabet.end()) {
            if (backwards_count_positions.find(numeric_portion) !=
                    backwards_count_positions.end())
                converted_code -= (static_cast<float>
                        (alphabet.at(non_numeric_portion)) /
                        arbitraryDivisor);
            else
                converted_code += (static_cast<float>
                        (alphabet.at(non_numeric_portion)) /
                        arbitraryDivisor);
        } else {
            return INVALID_SEQUENCE;
        }
        ordered_float_codes.push_back(converted_code);
    }

    std::sort(ordered_float_codes.begin(), ordered_float_codes.end());


    for (size_t i=0; i < ordered_float_codes.size(); i++) {
        int numeric_portion = std::round(ordered_float_codes[i]);
        int non_numeric_portion = std::round((ordered_float_codes[i] -
                    static_cast<float>(numeric_portion)) * arbitraryDivisor);
        std::string ostring = std::to_string(numeric_portion);

        if (non_numeric_portion == 0) {
            ordered_translated_codes.push_back(ostring);
            continue;
        }

        if (backwards_count_positions.find(numeric_portion) !=
                backwards_count_positions.end())
            non_numeric_portion = std::abs(non_numeric_portion);

        if (reverse_alphabet.find(non_numeric_portion) !=
                reverse_alphabet.end())
            ostring += reverse_alphabet.at(non_numeric_portion);
        else
            return INVALID_SEQUENCE;

        ordered_translated_codes.push_back(ostring);
    }

    return VALID_SEQUENCE;
}




/// Converts a list of sequences and a corresponding list
/// of annotations into an MSA.
/// Convenient for a smaller number of sequences that can fit in memory.
int build_msa_utility(std::vector<std::string> &sequences,
        std::vector<std::tuple<std::vector<std::string>, double, std::string, std::string>> &annotations,
        std::vector<std::string> &position_codes,
        std::vector<std::string> &aligned_seqs,
        const std::string &scheme,
        bool add_unobserved_positions) {
    if (sequences.size() != annotations.size() || sequences.size() == 0)
        throw std::runtime_error(std::string("The number of sequences and "
                    "annotations must match."));

    std::set<std::string> all_position_codes;
    int errCode;
    int chain_type = -1;

    // Add all codes to the set of observed codes, and check that only
    // one chain type is present.

    for (auto & annotation : annotations) {
        for (auto & pos_code : std::get<0>(annotation)) {
            if (pos_code != "-")
                all_position_codes.insert(pos_code);
        }
        if (std::get<2>(annotation) == "H") {
            if (chain_type > 0 && chain_type != MSA_HEAVY_CHAIN_ONLY) {
                throw std::runtime_error(std::string("An MSA can only be built "
                        "from either heavy chains or light chains. TCRs and "
                        "mAbs cannot be combined into an MSA."));
            } else {
                chain_type = MSA_HEAVY_CHAIN_ONLY;
            }
        } else if (std::get<2>(annotation) == "K" ||
                std::get<2>(annotation) == "L") {
            if (chain_type > 0 && chain_type != MSA_LIGHT_CHAIN_ONLY) {
                throw std::runtime_error(std::string("An MSA can only be built "
                        "from either heavy chains or light chains. TCRs and "
                        "mAbs cannot be combined into an MSA."));
            } else {
                chain_type = MSA_LIGHT_CHAIN_ONLY;
            }
        } else if (std::get<2>(annotation) == "B" ||
                std::get<2>(annotation) == "D") {
            if (chain_type > 0 && chain_type != MSA_TCR_LIGHT_CHAIN_ONLY) {
                throw std::runtime_error(std::string("An MSA can only be built "
                        "from either heavy chains or light chains. TCRs and "
                        "mAbs cannot be combined into an MSA."));
            } else {
                chain_type = MSA_TCR_LIGHT_CHAIN_ONLY;
            }
        } else if (std::get<2>(annotation) == "A" ||
                std::get<2>(annotation) == "G") {
            if (chain_type > 0 && chain_type != MSA_TCR_HEAVY_CHAIN_ONLY) {
                throw std::runtime_error(std::string("An MSA can only be built "
                        "from either heavy chains or light chains. TCRs and "
                        "mAbs cannot be combined into an MSA."));
            } else {
                chain_type = MSA_TCR_HEAVY_CHAIN_ONLY;
            }
        } else {
                throw std::runtime_error(std::string("An unrecognized chain "
                            "code was supplied."));
        }
    }

    // Check that there ARE position codes (other than '-', which is
    // ignored). We used to raise an exception if there were no
    // position codes but for the GUI and some user applications this
    // may be undesirable, so if this occurs, we now just return empty
    // lists.
    if (all_position_codes.size() == 0) {
        for (size_t i=0; i < sequences.size(); i++) {
            std::string empty_string = "";
            aligned_seqs.push_back(empty_string);
        }
        return VALID_SEQUENCE;
    }

    // If the user has requested addition of unobserved position codes,
    // add any as appropriate for the scheme. We have hard-coded the
    // expected positions for each scheme here, which is...not great...
    // but on the other hand, the numbering schemes are relatively ancient
    // and will not change at any time in the forseeable future,
    // so this may be ok for now.
    if (add_unobserved_positions) {
        std::vector<std::string> codes_to_add;
        // IMGT and AHO are chain invariant.
        if (scheme == "imgt") {
            for (int i=1; i <= 128; i++)
                all_position_codes.insert(std::to_string(i));
        } else if (scheme == "aho") {
            for (int i=1; i <= 149; i++)
                all_position_codes.insert(std::to_string(i));
        } else if (chain_type == MSA_LIGHT_CHAIN_ONLY) {
            // The other schemes sadly are not.
            if (scheme == "martin") {
                for (int i=1; i <= 107; i++)
                    all_position_codes.insert(std::to_string(i));
            } else if (scheme == "kabat") {
                for (int i=1; i <= 107; i++)
                    all_position_codes.insert(std::to_string(i));
            }
        } else if (chain_type == MSA_HEAVY_CHAIN_ONLY) {
            if (scheme == "martin") {
                for (int i=1; i <= 113; i++)
                    all_position_codes.insert(std::to_string(i));
            } else if (scheme == "kabat") {
                for (int i=1; i <= 113; i++)
                    all_position_codes.insert(std::to_string(i));
            }
        }
    }


    // Sort the position codes.

    std::vector<std::string> all_position_codes_vec(all_position_codes.begin(),
            all_position_codes.end());
    if (!sort_position_codes_utility(all_position_codes_vec,
            scheme, position_codes)) {
        throw std::runtime_error(std::string("Invalid position codes "
                    "were supplied."));
    }


    std::unordered_map<std::string, int> code_to_location;

    for (size_t i=0; i < position_codes.size(); i++)
        code_to_location[position_codes[i]] = i;

    for (size_t i=0; i < sequences.size(); i++) {
        std::string aligned_seq(position_codes.size(), '-');

        // The only situations where the annotation could be a different
        // length from the sequence are if the user made changes to an
        // annotation OR there was an alignment error (which is reported as
        // an error code). We could skip the sequence, but this creates more
        // problems than it solves -- same for throwing an exception. For now,
        // we report a blank alignment if this has occurred.
        if (std::get<0>(annotations[i]).size() != sequences[i].length()) {
            aligned_seqs.push_back(aligned_seq);
            continue;
        }

        for (size_t j=0; j < sequences[i].length(); j++) {
            std::string pos_code = std::get<0>(annotations[i])[j];
            if (pos_code == "-")
                continue;
            if (code_to_location.find(pos_code) == code_to_location.end())
                throw std::runtime_error(std::string("Invalid position codes "
                            "were supplied."));

            int location = code_to_location[pos_code];
            aligned_seq[location] = sequences[i][j];
        }
        aligned_seqs.push_back(aligned_seq);
    }
    return VALID_SEQUENCE;
}


// Trims an annotation / alignment to remove gaps at either end that correspond to
// non-numbered AAs.
int trim_alignment_utility(const std::string &sequence,
        std::tuple<std::vector<std::string>, double, std::string, std::string> &alignment,
        std::vector<std::string> &trimmed_alignment, int &exstart, int &exend,
        std::vector<char> &trimmed_seq) {
    if (std::get<0>(alignment).size() != sequence.length())
        return INVALID_SEQUENCE;

    exstart = -1;
    exend = -1;

    for (size_t i=0; i < std::get<0>(alignment).size(); i++) {
        if (std::get<0>(alignment)[i] == "-") {
            if (exstart >= 0 && exend < 0)
                exend = i;
            continue;
        }
        if (exstart < 0)
            exstart = i;
        trimmed_seq.push_back(sequence[i]);
        trimmed_alignment.push_back(std::get<0>(alignment)[i]);
    }

    if (exend < 0)
        exend = sequence.length();

    return VALID_SEQUENCE;
}

}  // namespace SequenceUtilities
