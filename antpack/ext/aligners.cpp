#include "aligners.h"

// Codes for the pathways that can link a score
// to the best-scoring parent.
#define LEFT_TRANSFER 0
#define DIAGONAL_TRANSFER 1
#define UP_TRANSFER 2


IMGTAligner::IMGTAligner(
                 py::array_t<double> scoreArray
):
    scoreArray(scoreArray)
{
    py::buffer_info info = scoreArray.request();
    if (info.shape.size() != 2){
        throw std::runtime_error(std::string("The scoreArray passed to "
                "IMGTAligner must be a 2d array"));
    }
    if (info.shape[0] < 1){
        throw std::runtime_error(std::string("The scoreArray passed to "
                "IMGTAligner must have >= 1 row"));
    }
    if (info.shape[1] != this->numAAs){
        throw std::runtime_error(std::string("The scoreArray passed to "
                "IMGTAligner must have 21 columns (1 per AA)"));
    }
    numPositions = info.shape[0];
}



std::vector<std::string> IMGTAligner::align(std::string query_sequence){

    if (query_sequence.length() == 0){
        throw std::runtime_error(std::string("An empty query "
                "sequence was passed to IMGTAligner."));
    }

    int rowSize = query_sequence.length() + 1;

    std::vector<std::string> numbering;
    int numElements = rowSize * (this->numPositions + 1);
    double *needleScores = new double[ numElements ];
    unsigned int *queryAsIdx = new unsigned int[query_sequence.length()];
    unsigned int *pathTrace = new unsigned int[ numElements ];
    auto scoreItr = this->scoreArray.unchecked<2>();

    // Translate the query sequence into an integer 0-21 encoding. This
    // is more verbose than using std::map but should be slightly faster
    // since compiler will convert to lookup table. If unexpected characters
    // are encountered, abort (PyBind11 will ensure the exception is passed
    // back to Python).
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
                delete[] queryAsIdx;
                delete[] needleScores;
                delete[] pathTrace;
                throw std::runtime_error(std::string("The query string "
                "passed to IMGTAligner.align contains nonstandard "
                        "amino acids or lowercase letters."));
                break;

        }
    }

    // Fill in the first row of the tables.
    needleScores[0] = 0;
    pathTrace[0] = 0;
    for (size_t i=0; i < query_sequence.length(); i++){
        needleScores[i+1] = needleScores[i] + scoreItr(0,20);
        pathTrace[i+1] = LEFT_TRANSFER;
    }

    // Now fill in the scoring grid, using the position-specific
    // scores for indels and substitutions, and add the appropriate
    // pathway traces.
    for (int i=1; i < this->numPositions + 1; i++){
        int gridPos = i * rowSize;
        int diagNeighbor = (i - 1) * rowSize;
        int upperNeighbor = diagNeighbor + 1;
        needleScores[gridPos] = diagNeighbor + scoreIter(i,20);
        pathTrace[gridPos] = UP_TRANSFER;
        gridPos++;

        for (size_t j=0; j < query_sequence.length(); j++){
            double lscore = needleScores[gridPos - 1] +
                    scoreItr(i,20);
            double uscore = needleScores[upperNeighbor] + 
                    scoreItr(i,20);
            double dscore = needleScores[diagNeighbor] + 
                    scoreItr(i,queryAsIdx[j]);

            // This is mildly naughty -- we don't consider the possibility
            // of a tie. Realistically, ties are going to be both extremely
            // rare and low-impact in this situation because of the way
            // the positions are scored. For now we are defaulting to "match"
            // if there is a tie.
            if (lscore > uscore && lscore > dscore){
                needleScores[gridPos] = lscore;
                pathTrace[gridPos] = LEFT_TRANSFER;
            }
            else if (uscore > dscore){
                needleScores[gridPos] = uscore;
                pathTrace[gridPos] = UP_TRANSFER;
            }
            else{
                needleScores[gridPos] = dscore;
                pathTrace[gridPos] = DIAGONAL_TRANSFER;
            }

            gridPos++;
            diagNeighbor++;
            upperNeighbor++;
        }
    }

    // Backtrace the path, and convert the backtrace to position numbering.
    int i = this->numPositions + 1;
    int j = query_sequence.length() + 1;
    for (int i=1; i < this->numPositions + 1; i++){
        int gridPos = i * query_sequence.length();
        int diagNeighbor = (i - 1) * query_sequence.length();
        int upperNeighbor = diagNeighbor + 1;
        needleScores[gridPos] = diagNeighbor + scoreIter(i,20);
        pathTrace[gridPos] = UP_TRANSFER;
        gridPos++;

        for (size_t j=0; j < query_sequence.length(); j++){
            double lscore = needleScores[gridPos - 1] +
                    scoreItr(i,20);
            double uscore = needleScores(upperNeighbor) + 
                    scoreItr(i,20);
            double dscore = needleScores[diagNeighbor] + 
                    scoreItr(i,queryAsIdx[j]);

            // This is mildly naughty -- we don't consider the possibility
            // of a tie. Realistically, ties are going to be both extremely
            // rare and low-impact in this situation because of the way
            // the positions are scored.
            if (lscore > uscore && lscore > dscore){
                needleScores[gridPos] = lscore;
                pathTrace[gridPos] = LEFT_TRANSFER;
            }
            else if (uscore > dscore){
                needleScores[gridPos] = uscore;
                pathTrace[gridPos] = UP_TRANSFER;
            }
            else{
                needleScores[gridPos] = dscore;
                pathTrace[gridPos] = DIAGONAL_TRANSFER;
            }

            gridPos++;
            diagNeighbor++;
            upperNeighbor++;
        }
    }
    

    delete[] queryAsIdx;
    delete[] needleScores;
    delete[] pathTrace;
    return numbering;
}