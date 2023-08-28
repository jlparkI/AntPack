#include "aligners.h"
#include "utilities.h"

// Codes for the pathways that can link a score
// to the best-scoring parent.
#define LEFT_TRANSFER 0
#define DIAGONAL_TRANSFER 1
#define UP_TRANSFER 2

// Codes for errors that may be encountered.
#define INVALID_SEQUENCE 0
#define NO_ERROR 1
#define FATAL_RUNTIME_ERROR 2

// These are "magic number" positions in the IMGT framework at
// which "forwards-backwards" insertion numbering must be applied.
// This is a nuisance, but is out of our control -- the IMGT #ing
// system has this quirk built-in... Note that because IMGT numbers
// from 1, these positions are the actual IMGT position - 1.
#define CDR1_INSERTION_PT 32
#define CDR2_INSERTION_PT 60
#define CDR3_INSERTION_PT 110




IMGTAligner::IMGTAligner(
                 py::array_t<double> scoreArray,
                 std::vector<std::vector<std::string>> consensus
):
    scoreArray(scoreArray)
{
    py::buffer_info info = scoreArray.request();
    // Note that exceptions thrown here are go back to Python via
    // PyBind as long as this constructor is used within the wrapper.
    if (info.shape.size() != 2){
        throw std::runtime_error(std::string("The scoreArray passed to "
                "IMGTAligner must be a 2d array"));
    }
    if (info.shape[1] != this->numAAs){
        throw std::runtime_error(std::string("The scoreArray passed to "
                "IMGTAligner must have 21 columns (1 per AA plus 1 for gap penalties)"));
    }
    if (info.shape[0] != NUM_IMGT_POSITIONS){
        throw std::runtime_error(std::string("The scoreArray passed to "
                "IMGTAligner must have the expected number of positions "
                "for the IMGT numbering system."));
    }
    if (consensus.size() != NUM_IMGT_POSITIONS){
        throw std::runtime_error(std::string("The consensus sequence passed to "
                "IMGTAligner must have the expected number of positions "
                "for the IMGT numbering system."));
    }

    // Update the consensus map to indicate which letters are expected at
    // which IMGT positions. An empty set at a given position means that
    // ANY letter at that position is considered acceptable. These are
    // excluded from percent identity calculation. numRestrictedPositions
    // indicates how many are INCLUDED in percent identity.
    numRestrictedPositions = 0;

    for (size_t i = 0; i < consensus.size(); i++){
        if (consensus[i].empty()){
            continue;
        }
        numRestrictedPositions += 1;
        for (std::string & AA : consensus[i]){
            if (!validate_sequence(AA)){
                throw std::runtime_error(std::string("Non-default AAs supplied "
                    "in a consensus file."));
            }
            if (AA.length() > 1){
                throw std::runtime_error(std::string("Non-default AAs supplied "
                    "in a consensus file."));
            }
            consensusMap[i].insert(AA[0]);
        }
    }

    numPositions = info.shape[0];
}


std::tuple<std::vector<std::string>, double, int> IMGTAligner::align(std::string query_sequence){

    std::vector<std::string> finalNumbering;
    double percentIdentity = 0;
    if (query_sequence.length() == 0){
        return std::tuple<std::vector<std::string>, double, int>{finalNumbering,
                percentIdentity, INVALID_SEQUENCE};
    }

    std::vector<int> positionKey;

    int rowSize = query_sequence.length() + 1;
    int numElements = rowSize * (this->numPositions + 1);

    double *needleScores = new double[ numElements ];
    int *queryAsIdx = new int[query_sequence.length()];
    int *pathTrace = new int[ numElements ];
    int unsigned *initNumbering = new unsigned int[ this->numPositions ];
    auto scoreItr = this->scoreArray.unchecked<2>();


    if (!convert_sequence_to_array(queryAsIdx, query_sequence)){
        delete[] queryAsIdx;
        delete[] needleScores;
        delete[] pathTrace;
        delete[] initNumbering;
        return std::tuple<std::vector<std::string>, double, int>{finalNumbering,
                        percentIdentity, INVALID_SEQUENCE};
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
        needleScores[gridPos] = needleScores[diagNeighbor] + scoreItr(i-1,20);
        pathTrace[gridPos] = UP_TRANSFER;
        gridPos++;

        for (size_t j=0; j < query_sequence.length(); j++){
            double lscore = needleScores[gridPos - 1] +
                    scoreItr(i-1,20);
            double uscore = needleScores[upperNeighbor] + 
                    scoreItr(i-1,20);
            double dscore = needleScores[diagNeighbor] + 
                    scoreItr(i-1,queryAsIdx[j]);

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

    // Set all members of initNumbering array to 0.
    for (int i=0; i < this->numPositions; i++){
        initNumbering[i] = 0;
    }

    // Backtrace the path, and determine how many insertions there
    // were at each position where there is an insertion. Determine
    // how many gaps exist at the beginning and end of the sequence
    // as well (the "N-terminus" and "C-terminus").
    int i = this->numPositions;
    int j = query_sequence.length();
    int numNTermGaps = 0;
    int numCTermGaps = 0;

    while (i > 0 || j > 0){
        int gridPos = i * rowSize + j;
        switch (pathTrace[gridPos]){
            case LEFT_TRANSFER:
                if (i > 0 && i < this->numPositions){
                    initNumbering[i-1] += 1;
                }
                else if (i == 0){
                    numNTermGaps += 1;
                }
                else{
                    numCTermGaps += 1;
                }
                j -= 1;
                break;
            case UP_TRANSFER:
                i -= 1;
                break;
            case DIAGONAL_TRANSFER:
                initNumbering[i-1] += 1;
                i -= 1;
                j -= 1;
                break;
            // This should never happen -- if it does, report it as a
            // SERIOUS error.
            default:
                delete[] queryAsIdx;
                delete[] needleScores;
                delete[] pathTrace;
                delete[] initNumbering;
                return std::tuple<std::vector<std::string>, double, int>{finalNumbering,
                        percentIdentity, FATAL_RUNTIME_ERROR};
                break;
        }
        // This should never happen, but just in case...
        if (i < 0){
            i = 0;
        }
        if (j < 0){
            j = 0;
        }
    }

    // Add gaps at the beginning of the numbering corresponding to the
    // number of N-terminal insertions. This ensures the numbering has
    // the same length as the input sequence.
    // positionKey is set to -1 to indicate these positions should not
    // be considered for percent identity calculation.
    for (i=0; i < numNTermGaps; i++){
        finalNumbering.push_back("-");
        positionKey.push_back(-1);
    }

    // Build vector of IMGT numbers. Unfortunately the IMGT system adds insertion codes
    // forwards then backwards where there is > 1 insertion at a given position,
    // but ONLY if the insertion is at an expected position in the CDRs!
    // Everywhere else, insertions are recognized by just placing in
    // the usual order (?!!!?) This annoying quirk adds some complications.
    // We also create the positionKey here which we can quickly use to determine
    // percent identity by referencing this->consensusMap, which indicates
    // what AAs are expected at each position.
    for (i=0; i < this->numPositions; i++){
        if (initNumbering[i] == 0){
            continue;
        }
        // The +1 here and hereafter is because IMGT numbers from 1.
        // Position key maps the IMGT numbering to integers that access
        // this->consensusMap.
        finalNumbering.push_back(std::to_string(i+1));
        positionKey.push_back(i);
        if (initNumbering[i] == 1){
            continue;
        }
        // These are the positions for which we need to deploy forward-backward
        // lettering (because the IMGT system is weird, but also the default...)
        // positionKey is set to -1 for all insertions to indicate these should
        // not be considered for calculating percent identity.
        switch (i){
            case CDR1_INSERTION_PT:
            case CDR2_INSERTION_PT:
            case CDR3_INSERTION_PT:
                if (initNumbering[i] > this->alphabet.size()){
                    delete[] queryAsIdx;
                    delete[] needleScores;
                    delete[] pathTrace;
                    delete[] initNumbering;
                    return std::tuple<std::vector<std::string>, double, int>{finalNumbering,
                        percentIdentity, FATAL_RUNTIME_ERROR};
                }
                int ceil_cutpoint, floor_cutpoint;
                ceil_cutpoint = initNumbering[i] / 2;
                floor_cutpoint = (initNumbering[i] - 1) / 2;
                for (j=0; j < floor_cutpoint; j++){
                    finalNumbering.push_back(std::to_string(i+1) + this->alphabet[j]);
                    positionKey.push_back(-1);
                }
                for (j=ceil_cutpoint; j > 0; j--){
                    finalNumbering.push_back(std::to_string(i+2) + this->alphabet[j-1]);
                    positionKey.push_back(-1);
                }
                break;
            default:
                if (initNumbering[i] > this->alphabet.size()){
                    delete[] queryAsIdx;
                    delete[] needleScores;
                    delete[] pathTrace;
                    delete[] initNumbering;
                    return std::tuple<std::vector<std::string>, double, int>{finalNumbering,
                        percentIdentity, FATAL_RUNTIME_ERROR};
                }
                for (j=1; j < initNumbering[i]; j++){
                    finalNumbering.push_back(std::to_string(i+1) + this->alphabet[j-1]);
                    positionKey.push_back(-1);
                }
                break;
        }
    }


    // Add gaps at the end of the numbering corresponding to the
    // number of C-terminal insertions. This ensures the numbering has
    // the same length as the input sequence.
    for (i=0; i < numCTermGaps; i++){
        finalNumbering.push_back("-");
        positionKey.push_back(-1);
    }

    // positionKey is now the same length as finalNumbering and indicates at
    // each position whether that position maps to a standard 1 - 128 IMGT
    // position and if so what. Now use this to calculate percent identity,
    // excluding those positions at which IMGT tolerates any amino acid.
    for (size_t k=0; k < positionKey.size(); k++){
        int imgtStdPosition = positionKey[k];
        if (imgtStdPosition < 0){
            continue;
        }
        if (this->consensusMap[imgtStdPosition].empty()){
            continue;
        }
        if (this->consensusMap[imgtStdPosition].find(query_sequence[k]) !=
                this->consensusMap[imgtStdPosition].end()){
            percentIdentity += 1;
        }
    }
    
    percentIdentity /= this->numRestrictedPositions;

    delete[] queryAsIdx;
    delete[] needleScores;
    delete[] pathTrace;
    delete[] initNumbering;
    return std::tuple<std::vector<std::string>, double, int>{finalNumbering,
                    percentIdentity, NO_ERROR};
}