#include "cterm_finder.h"



CTermFinder::CTermFinder(std::string consensus_filepath){
    std::filesystem::path extensionPath = consensus_filepath;
    std::string npyFName = "CTERMFINDER_CONSENSUS_H.npy";
    std::filesystem::path npyFPath = extensionPath / npyFName;
    cnpy::NpyArray raw_score_arr;

    try{
        raw_score_arr = cnpy::npy_load(npyFPath.string());
    }
    catch (...){
        throw std::runtime_error(std::string("The consensus file / library installation "
                            "has an issue."));
    }

    std::vector<std::string> boundary_chains = {"H", "K", "L"};
    this->score_arr_shape[0] = raw_score_arr.shape[0];
    this->score_arr_shape[1] = raw_score_arr.shape[1];
    this->score_arr_shape[2] = boundary_chains.size();
    if (this->score_arr_shape[1] != 20){
        throw std::runtime_error(std::string("The consensus file / library installation "
                            "has an issue."));
    }
    this->score_array = std::make_unique<double[]>( this->score_arr_shape[0] *
           this->score_arr_shape[1] * this->score_arr_shape[2] );

    for (size_t i=0; i < boundary_chains.size(); i++){
        std::string npyFName = "CTERMFINDER_CONSENSUS_" + boundary_chains[i] + ".npy";
        std::filesystem::path npyFPath = extensionPath / npyFName;

        try{
            raw_score_arr = cnpy::npy_load(npyFPath.string());
        }
        catch (...){
            throw std::runtime_error(std::string("The consensus file / library installation "
                            "has an issue."));
        }
        int ld_shape0 = raw_score_arr.shape[0], ld_shape1 = raw_score_arr.shape[1];
        if (raw_score_arr.word_size != 8 || ld_shape0 != this->score_arr_shape[0] ||
                ld_shape1 != this->score_arr_shape[1]){
            throw std::runtime_error(std::string("The consensus file / library installation "
                            "has an issue."));
        }

        double *raw_score_ptr = raw_score_arr.data<double>();

        for (int j=0; j < this->score_arr_shape[0]; j++){
            for (int k=0; k < this->score_arr_shape[1]; k++){
                score_array[j*this->score_arr_shape[1]*this->score_arr_shape[2] +
                    k*this->score_arr_shape[2] + i] = *raw_score_ptr;
                raw_score_ptr++;
            }
        }
    }
    this->num_positions = this->score_arr_shape[0];
}








// Does a fast (relative to typical alignments) search for templates
// which indicate the c-terminal of a kappa, lambda or heavy chain.
int CTermFinder::find_c_terminals(std::string query_sequence,
        std::array<double, 3> &best_scores, std::array<int, 3> &best_positions){

    std::string errorMessage;
    size_t numQueryAlignments = query_sequence.length() - this->num_positions;

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
       size_t score_arr_row = 0;

       for (int k=0; k < 3; k++)
           match_score[k] = 0;

       for (int j=0; j < this->num_positions; j++){
            for (int k=0; k < 3; k++){
                size_t position = score_arr_row + encoded_sequence[i+j] * this->score_arr_shape[2] + k;
                match_score[k] += this->score_array[position];
            }
            score_arr_row += (this->score_arr_shape[1] * this->score_arr_shape[2]);
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
