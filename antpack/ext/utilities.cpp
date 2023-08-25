/* Utilities, tools and other extensions for sequence processing.
*/
#include "utilities.h"

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