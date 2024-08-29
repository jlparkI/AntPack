#include "cterm_finder.h"



CTermFinder::CTermFinder(
                 py::array_t<double, py::array::c_style> score_array
):
    score_array(score_array)
{
    py::buffer_info info = score_array.request();
    // Note that exceptions thrown here are go back to Python via
    // PyBind as long as this constructor is used within the wrapper.
    if (info.shape.size() != 3){
        throw std::runtime_error(std::string("The score_array passed to "
                "CTermFinder must be a 3d array"));
    }
    if (info.shape[1] != 20){
        throw std::runtime_error(std::string("The score_array passed to "
                "CTermFinder must have shape[1] == 20 (1 per AA)"));
    }
    if (info.shape[2] != 3){
        throw std::runtime_error(std::string("The score_array passed to "
                "CTermFinder must have shape[2] == 3 (1 per chain type)"));
    }


    this->num_positions = info.shape[0];
}








// Does a fast (relative to typical alignments) search for templates
// which indicate the c-terminal of a kappa, lambda or heavy chain.
int CTermFinder::find_c_terminals(std::string query_sequence,
        std::array<double, 3> &best_scores, std::array<int, 3> &best_positions){

    std::string errorMessage;
    size_t numQueryAlignments = query_sequence.length() - this->num_positions;
    auto scoreMatItr = this->score_array.unchecked<3>();

    for (int k=0; k < 3; k++){
        best_scores[k] = 0;
        best_positions[k] = 0;
    }

    if (query_sequence.length() < MINIMUM_SEQUENCE_LENGTH ||
            numQueryAlignments <= 0)
        return INVALID_SEQUENCE;

    auto encoded_sequence = std::make_unique<int[]>( query_sequence.length() );

    if (!convert_sequence_to_array(encoded_sequence.get(), query_sequence))
        return INVALID_SEQUENCE;


    // Look for the first two cysteines in the sequence,
    // then start searching for template matches which occur
    // after this.
    size_t startPosition = 0, ncysteines = 0;

    for (size_t i=0; i < query_sequence.length(); i++){
        if (encoded_sequence[i] == 1){
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
    if ( ncysteines < 2 || startPosition >= numQueryAlignments )
        return INVALID_SEQUENCE;


    for (size_t i=startPosition; i < numQueryAlignments; i++){
       double match_score[3];
       for (int k=0; k < 3; k++)
           match_score[k] = 0;

       for (int j=0; j < this->num_positions; j++){
            for (int k=0; k < 3; k++)
                match_score[k] += scoreMatItr(j, encoded_sequence[i+j], k);
       }

       for (int k=0; k < 3; k++){
           if (match_score[k] > best_scores[k]){
               best_scores[k] = match_score[k];
               best_positions[k] = i;
           }
       }
    }

    return VALID_SEQUENCE;
}
