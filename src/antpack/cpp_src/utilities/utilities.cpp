/* Utilities, tools and other extensions for sequence processing.
*/
#include "utilities.h"


// Codes for sequence validation.
#define VALID_SEQUENCE 1
#define INVALID_SEQUENCE 0

// Codes for chain type identification when building MSAs.
#define MSA_HEAVY_CHAIN_ONLY 1
#define MSA_LIGHT_CHAIN_ONLY 2


//Checks a query sequence to ensure it contains only recognized
//amino acids.
int validate_sequence(std::string query_sequence){
    bool validQuery = true;

    if (query_sequence.length() == 0){
        return INVALID_SEQUENCE;
    }

    for (char & c : query_sequence)
    {
        //Switch will provide a (small) speed gain over map
        //since compiler should convert it to a lookup table
        switch (c){
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
                validQuery = false;
                break;

        }
        if (!validQuery){
            break;
        }
    }
    if (!validQuery){
        return INVALID_SEQUENCE;
    }
    return VALID_SEQUENCE;
}


// Converts an input sequence to a numeric representation as an
// array of integers. This function does not do any safety checks
// on the input so it should never be accessed by Python functions --
// only ever by C++ functions that DO perform safety checks on
// input. In particular, caller MUST verify that seqAsArray is
// the same length as sequence.
int convert_sequence_to_array(int *queryAsIdx, std::string query_sequence){
    
    // Translate the query sequence into an integer 0-21 encoding. This
    // is more verbose than using std::map but should be slightly faster
    // since compiler will convert to lookup table. If unexpected characters
    // are encountered, abort and return an error code.
    for (size_t i=0; i < query_sequence.length(); i++){
        switch (query_sequence[i]){
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
            case '-':
                queryAsIdx[i] = 20;
                break;

            default:
                return INVALID_SEQUENCE;
                break;

        }
    }
    return VALID_SEQUENCE;
}




//Takes as input a list of position codes (e.g. IMGT) and sorts it in a way
//that respects the properties of the scheme. If any invalid position
//codes are encountered, the function will stop and return an invalid
//error code. Note that since this function constructs a map on
//each call, calling it 1000s of times is not efficient or advisable.
//Indeed, since it should only be necessary to call it once per
//dataset, that should never actually happen.
int sort_position_codes_utility(std::vector<std::string> &position_codes,
        std::string scheme, std::vector<std::string> &orderedTranslatedCodes){
    std::vector<float> orderedFloatCodes;
    float arbitraryDivisor = 1000.0;
    int max_position = 128;
    if (scheme == "aho")
        max_position = 149;

    if (position_codes.size() == 0)
        return INVALID_SEQUENCE;

    // Map from a standard alphabet Our preference would be to number insertions
    // as _1, _2 etc, but most numbering programs use letters, so we do the same here
    // for consistency and ease of comparison. This does limit the number of possible
    // insertions BUT there should never be a case where the user actually needs more than
    // 78 insertion codes (if they really do there's something wrong with their
    // alignments and/or it's not really an antibody...)
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
    std::map<int, std::string> reverseAlphabet;
    for ( const auto &p : alphabet )
        reverseAlphabet[p.second] = p.first;


    // The IMGT scheme has several positions at which we have to "count backwards", which is an annoying quirk of
    // that scheme. We have to hard-code these unfortunately.
    std::set<int> backwardsCountPositions;
    if (scheme == "imgt")
        backwardsCountPositions = {33, 61, 112};


    for (size_t i=0; i < position_codes.size(); i++){
        std::string inputCode = position_codes[i];
        std::string nonNumericPortion;
        int numericPortion;
        float convertedCode;

        // This will throw if the string starts with a letter but will otherwise extract the integer piece,
        // which rules out certain kinds of invalid codes the user might pass in error (e.g. '-', 'A128').
        try{
            numericPortion = std::stoi(inputCode);
        }
        catch (...){
            return INVALID_SEQUENCE;
        }
        if (numericPortion <= 0 || numericPortion > max_position)
            return INVALID_SEQUENCE;

        convertedCode = static_cast<float>(numericPortion);

        for (char & c : inputCode) {
	        if (!std::isdigit(c))
                nonNumericPortion += c;
        }

        if (nonNumericPortion.length() == 0){
            orderedFloatCodes.push_back(convertedCode);
            continue;
        }
        
        if (alphabet.find(nonNumericPortion)!=alphabet.end()){
            if (backwardsCountPositions.find(numericPortion)!=backwardsCountPositions.end())
                convertedCode -= (static_cast<float>(alphabet.at(nonNumericPortion)) / arbitraryDivisor);
            else
                convertedCode += (static_cast<float>(alphabet.at(nonNumericPortion)) / arbitraryDivisor);
        }
        else
            return INVALID_SEQUENCE;

        orderedFloatCodes.push_back(convertedCode);
    }

    std::sort(orderedFloatCodes.begin(), orderedFloatCodes.end());


    for (size_t i=0; i < orderedFloatCodes.size(); i++){
        int numericPortion = std::round(orderedFloatCodes[i]);
        int nonNumericPortion = std::round((orderedFloatCodes[i] - numericPortion) * arbitraryDivisor);
        std::string ostring = std::to_string(numericPortion);

        if (nonNumericPortion == 0){
            orderedTranslatedCodes.push_back(ostring);
            continue;
        }

        if (backwardsCountPositions.find(numericPortion)!=backwardsCountPositions.end())
            nonNumericPortion = std::abs(nonNumericPortion);
        
        if (reverseAlphabet.find(nonNumericPortion)!=reverseAlphabet.end())
            ostring += reverseAlphabet.at(nonNumericPortion);
        else
            return INVALID_SEQUENCE;

        orderedTranslatedCodes.push_back(ostring);
    }

    return VALID_SEQUENCE;
}




// Converts a list of sequences and a corresponding list of annotations into an MSA.
// Convenient for a smaller number of sequences that can fit in memory.
int build_msa_utility(std::vector<std::string> &sequences,
        std::vector<std::tuple<std::vector<std::string>, double, std::string, std::string>> &annotations,
        std::vector<std::string> &positionCodes,
        std::vector<std::string> &alignedSeqs,
        const std::string &scheme){
    if (sequences.size() != annotations.size() || sequences.size() == 0)
        throw std::runtime_error(std::string("The number of sequences and annotations must match."));

    std::set<std::string> allPositionCodes;
    int errCode;
    int chainType = -1;

    for (auto &annotation : annotations){
        for (auto &posCode : std::get<0>(annotation)){
            if (posCode != "-")
                allPositionCodes.insert(posCode);
        }
        if (std::get<2>(annotation) == "H"){
            if (chainType > 0 && chainType != MSA_HEAVY_CHAIN_ONLY){
                throw std::runtime_error(std::string("An MSA can only be built "
                            "from one chain type."));
            }
            else
                chainType = MSA_HEAVY_CHAIN_ONLY;
        }
        else if (std::get<2>(annotation) == "K" || std::get<2>(annotation) == "L"){
            if (chainType > 0 && chainType != MSA_LIGHT_CHAIN_ONLY){
                throw std::runtime_error(std::string("An MSA can only be built "
                            "from one chain type."));
            }
            else
                chainType = MSA_LIGHT_CHAIN_ONLY;
        }
    }

    std::vector<std::string> allPositionCodesVec(allPositionCodes.begin(), allPositionCodes.end()); 
    errCode = sort_position_codes_utility(allPositionCodesVec, scheme, positionCodes);
    if (errCode != VALID_SEQUENCE)
        throw std::runtime_error(std::string("Invalid position codes were supplied."));


    std::unordered_map<std::string, int> codeToLocation;

    for (size_t i=0; i < positionCodes.size(); i++)
        codeToLocation[positionCodes[i]] = i;

    for (size_t i=0; i < sequences.size(); i++){
        std::string alignedSeq(positionCodes.size(), '-');

        for (size_t j=0; j < sequences[i].length(); j++){
            std::string posCode = std::get<0>(annotations[i])[j];
            if (codeToLocation.find(posCode) == codeToLocation.end())
                throw std::runtime_error(std::string("Invalid position codes were supplied."));
            
            int location = codeToLocation[posCode];
            alignedSeq[location] = sequences[i][j];
        }
        alignedSeqs.push_back(alignedSeq);
    }
    return VALID_SEQUENCE;
}


// Trims an annotation / alignment to remove gaps at either end that correspond to
// non-numbered AAs.
int trim_alignment_utility(const std::string &sequence,
        std::tuple<std::vector<std::string>, double, std::string, std::string> &alignment,
        std::vector<std::string> &trimmedAlignment, int &exstart, int &exend,
        std::vector<char> &trimmedSeq){
    if (std::get<0>(alignment).size() != sequence.length())
        return INVALID_SEQUENCE;

    exstart = -1;
    exend = -1;

    for (size_t i=0; i < std::get<0>(alignment).size(); i++){
        if (std::get<0>(alignment)[i] == "-"){
            if (exstart >= 0 && exend < 0)
                exend = i;
            continue;
        }
        if (exstart < 0)
            exstart = i;
        trimmedSeq.push_back(sequence.at(i));
        trimmedAlignment.push_back(std::get<0>(alignment)[i]);
    }

    if (exend < 0)
        exend = sequence.length();

    return VALID_SEQUENCE;
}