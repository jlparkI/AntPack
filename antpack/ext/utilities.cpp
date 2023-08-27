/* Utilities, tools and other extensions for sequence processing.
*/
#include "utilities.h"


// Codes for sequence validation.
#define VALID_SEQUENCE 1
#define INVALID_SEQUENCE 0


//Checks a query sequence to ensure it contains only recognized
//amino acids. This function is directly available to Python
// callers through the wrapper.
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