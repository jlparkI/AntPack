#include "cterm_finder.h"
#include "utilities.h"


CTermFinder::CTermFinder(
                 py::array_t<double, py::array::c_style> scoreArray
):
    scoreArray(scoreArray)
{
    py::buffer_info info = scoreArray.request();
    // Note that exceptions thrown here are go back to Python via
    // PyBind as long as this constructor is used within the wrapper.
    if (info.shape.size() != 3){
        throw std::runtime_error(std::string("The scoreArray passed to "
                "CTermFinder must be a 3d array"));
    }
    if (info.shape[1] != 20){
        throw std::runtime_error(std::string("The scoreArray passed to "
                "CTermFinder must have shape[1] == 20 (1 per AA)"));
    }
    if (info.shape[2] != 3){
        throw std::runtime_error(std::string("The scoreArray passed to "
                "CTermFinder must have shape[2] == 3 (1 per chain type)"));
    }


    numPositions = info.shape[0];
}









// Does a fast (relative to typical alignments) search for templates
// which indicate the c-terminal of a kappa, lambda or heavy chain.
std::string CTermFinder::find_c_terminals(std::string query_sequence,
        py::array_t<double, py::array::c_style> bestScores,
        py::array_t<int32_t, py::array::c_style> bestPositions){

    cTermAllowedErrorCodes errorCode;
    std::string errorMessage;
    size_t numQueryAlignments = query_sequence.length() - this->numPositions;
    auto scoreMatItr = this->scoreArray.unchecked<3>();
    auto bestScoresItr = bestScores.mutable_unchecked<1>();
    auto bestPositionsItr = bestPositions.mutable_unchecked<1>();

    if (bestScores.shape(0) != 3 || bestPositions.shape(0) != 3){
        throw std::runtime_error(std::string("Arguments passed to CTermFinder do "
                    "not have the correct shape.")); 
    }

    for (int k=0; k < 3; k++){
        bestScoresItr[k] = 0;
        bestPositionsItr[k] = 0;
    }

    if (query_sequence.length() < MINIMUM_SEQUENCE_LENGTH ||
            numQueryAlignments <= 0){
        errorCode = invalidSequence;
        errorMessage = this->errorCodeToMessage[errorCode];
        return errorMessage;
    }

    auto queryAsIdx = std::make_unique<int[]>( query_sequence.length() );

    if (!convert_sequence_to_array(queryAsIdx.get(), query_sequence)){
        errorCode = invalidSequence;
        errorMessage = this->errorCodeToMessage[errorCode];
        return errorMessage;
    }


    // Look for the first two cysteines in the sequence,
    // then start searching for template matches which occur
    // after this.
    size_t startPosition = 0, ncysteines = 0;

    for (size_t i=0; i < query_sequence.length(); i++){
        if (queryAsIdx[i] == 1){
            ncysteines += 1;
            if (ncysteines >= 2){
                startPosition = i + 1;
                break;
            }
        }
    }

    // If we could not find two cysteines, or if the resulting starting
    // position is too close to the end of the sequence, there is
    // something wrong; abort.
    if ( ncysteines < 2 || startPosition >= numQueryAlignments ){
        errorCode = invalidSequence;
        errorMessage = this->errorCodeToMessage[errorCode];
        return errorMessage;
    }


    for (size_t i=startPosition; i < numQueryAlignments; i++){
       double matchScore[3];
       for (int k=0; k < 3; k++)
           matchScore[k] = 0;

       for (int j=0; j < this->numPositions; j++){
            for (int k=0; k < 3; k++)
                matchScore[k] += scoreMatItr(j, queryAsIdx[i+j], k);
       }

       for (int k=0; k < 3; k++){
           if (matchScore[k] > bestScoresItr(k)){
               bestScoresItr[k] = matchScore[k];
               bestPositionsItr[k] = i;
           }
       }
    }

    errorCode = noError;
    errorMessage = this->errorCodeToMessage[errorCode];
    return errorMessage;
}
