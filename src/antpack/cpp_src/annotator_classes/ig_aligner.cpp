#include "ig_aligner.h"


IGAligner::IGAligner(
                 py::array_t<double, py::array::c_style> score_array,
                 std::vector<std::vector<std::string>> consensus,
                 std::string chain_name,
                 std::string scheme,
                 double terminalTemplateGapPenalty,
                 double CterminalQueryGapPenalty,
                 bool compress_initial_gaps
):
    score_array(score_array),
    chain_name(chain_name),
    scheme(scheme),
    terminalTemplateGapPenalty(terminalTemplateGapPenalty),
    CterminalQueryGapPenalty(CterminalQueryGapPenalty),
    compress_initial_gaps(compress_initial_gaps)
{
    py::buffer_info info = score_array.request();
    // Note that exceptions thrown here are go back to Python via
    // PyBind as long as this constructor is used within the wrapper.
    if (info.shape.size() != 2){
        throw std::runtime_error(std::string("The score_array passed to "
                "IGAligner must be a 2d array"));
    }
    if (info.shape[1] != NUM_AAS){
        throw std::runtime_error(std::string("The score_array passed to "
                "IGAligner must have 22 columns (1 per AA plus 2 for gap penalties)"));
    }

    // Check that the number of positions in the scoring and consensus is as expected,
    // and set the list of highly conserved positions according to the selected scheme.
    if (scheme == "imgt"){
        this->highly_conserved_positions = {HIGHLY_CONSERVED_IMGT_1, HIGHLY_CONSERVED_IMGT_2,
                                        HIGHLY_CONSERVED_IMGT_3, HIGHLY_CONSERVED_IMGT_4,
                                        HIGHLY_CONSERVED_IMGT_5, HIGHLY_CONSERVED_IMGT_6};
        if (info.shape[0] != NUM_HEAVY_IMGT_POSITIONS && info.shape[0] !=
                NUM_LIGHT_IMGT_POSITIONS){
            throw std::runtime_error(std::string("The score_array passed to "
                "IGAligner must have the expected number of positions "
                "for the numbering system."));
        }
        if (consensus.size() != NUM_HEAVY_IMGT_POSITIONS && consensus.size() !=
                NUM_LIGHT_IMGT_POSITIONS){
            throw std::runtime_error(std::string("The consensus sequence passed to "
                "IGAligner must have the expected number of positions "
                "for the numbering system."));
        }
    }
    else if (scheme == "aho"){
        this->highly_conserved_positions = {HIGHLY_CONSERVED_AHO_1, HIGHLY_CONSERVED_AHO_2,
                                        HIGHLY_CONSERVED_AHO_3, HIGHLY_CONSERVED_AHO_4,
                                        HIGHLY_CONSERVED_AHO_5, HIGHLY_CONSERVED_AHO_6};
        if (info.shape[0] != NUM_HEAVY_AHO_POSITIONS && info.shape[0] !=
                NUM_LIGHT_AHO_POSITIONS){
            throw std::runtime_error(std::string("The score_array passed to "
                "IGAligner must have the expected number of positions "
                "for the numbering system."));
        }
        if (consensus.size() != NUM_HEAVY_AHO_POSITIONS && consensus.size() !=
                NUM_LIGHT_AHO_POSITIONS){
            throw std::runtime_error(std::string("The consensus sequence passed to "
                "IGAligner must have the expected number of positions "
                "for the numbering system."));
        }
    }
    else if (scheme == "martin" || scheme == "kabat"){
        if (chain_name == "L" || chain_name == "K"){
            this->highly_conserved_positions = {HIGHLY_CONSERVED_KABAT_LIGHT_1,
                    HIGHLY_CONSERVED_KABAT_LIGHT_2, HIGHLY_CONSERVED_KABAT_LIGHT_3,
                    HIGHLY_CONSERVED_KABAT_LIGHT_4, HIGHLY_CONSERVED_KABAT_LIGHT_5,
                    HIGHLY_CONSERVED_KABAT_LIGHT_6};
        }
        else if (chain_name == "H"){
            this->highly_conserved_positions = {HIGHLY_CONSERVED_KABAT_HEAVY_1,
                    HIGHLY_CONSERVED_KABAT_HEAVY_2, HIGHLY_CONSERVED_KABAT_HEAVY_3,
                    HIGHLY_CONSERVED_KABAT_HEAVY_4, HIGHLY_CONSERVED_KABAT_HEAVY_5,
                    HIGHLY_CONSERVED_KABAT_HEAVY_6};
        }
        else{
            throw std::runtime_error(std::string("Martin, Kabat currently "
                        "support only H, K, L chains."));
        }
        if (info.shape[0] != NUM_HEAVY_MARTIN_KABAT_POSITIONS && info.shape[0] !=
                NUM_LIGHT_MARTIN_KABAT_POSITIONS){
            throw std::runtime_error(std::string("The score_array passed to "
                "IGAligner must have the expected number of positions "
                "for the numbering system."));
        }
        if (consensus.size() != NUM_HEAVY_MARTIN_KABAT_POSITIONS && consensus.size() !=
                NUM_LIGHT_MARTIN_KABAT_POSITIONS){
            throw std::runtime_error(std::string("The consensus sequence passed to "
                "IGAligner must have the expected number of positions "
                "for the numbering system."));
        }
    }
    else{
        throw std::runtime_error(std::string("Currently IGAligner only recognizes "
                    "schemes 'martin', 'kabat', 'imgt', 'aho'."));
    }

    num_positions = info.shape[0];
    // Update the consensus map to indicate which letters are expected at
    // which positions. An empty set at a given position means that
    // ANY letter at that position is considered acceptable. These are
    // excluded from percent identity calculation. num_restricted_positions
    // indicates how many are INCLUDED in percent identity.
    num_restricted_positions = 0;

    for (size_t i = 0; i < consensus.size(); i++){
        consensus_map.push_back(std::set<char> {});
        if (consensus[i].empty()){
            continue;
        }
        num_restricted_positions += 1;
        for (std::string & AA : consensus[i]){
            if (!validate_sequence(AA)){
                throw std::runtime_error(std::string("Non-default AAs supplied "
                    "in a consensus file."));
            }
            if (AA.length() > 1){
                throw std::runtime_error(std::string("Non-default AAs supplied "
                    "in a consensus file."));
            }
            consensus_map[i].insert(AA[0]);
        }
    }
}


// Convenience function for retrieving the chain name.
std::string IGAligner::get_chain_name(){
    return this->chain_name;
}





// This core alignment function is wrapped by other alignment functions.
// Previously this function was available to Python but is no longer made
// available; indeed, it now expects queryAsidx, which is an int array that
// must be of the same length as query_sequence -- caller must guarantee this.
// Therefore this function is now accessed only by the SingleChainAnnotator /
// PairedChainAnnotator classes (for an align function that can be accessed
// directly by Python wrappers, see align_test_only).
void IGAligner::align(std::string query_sequence,
                int *encoded_sequence, std::vector<std::string> &final_numbering,
                double &percent_identity, std::string &error_message){

    allowedErrorCodes error_code;
    percent_identity = 0;


    if (query_sequence.length() < MINIMUM_SEQUENCE_LENGTH){
        error_code = invalidSequence;
        error_message = this->error_code_to_message[error_code];
        return;
    }

    std::vector<int> position_key;

    int row_size = query_sequence.length() + 1;
    int numElements = row_size * (this->num_positions + 1);

    auto needle_scores = std::make_unique<double[]>( numElements );
    auto path_trace = std::make_unique<uint8_t[]>( numElements );
    auto init_numbering = std::make_unique<unsigned int[]>( this->num_positions );

    // Fill in the scoring table.
    fill_needle_scoring_table(needle_scores.get(), path_trace.get(),
                    query_sequence.length(), row_size, encoded_sequence);

    // Set all members of init_numbering array to 0.
    for (int i=0; i < this->num_positions; i++)
        init_numbering[i] = 0;

    // Backtrace the path, and determine how many insertions there
    // were at each position where there is an insertion. Determine
    // how many gaps exist at the beginning and end of the sequence
    // as well (the "N-terminus" and "C-terminus").
    int i = this->num_positions;
    int j = query_sequence.length();
    int numNTermGaps = 0;
    int numCTermGaps = 0;

    while (i > 0 || j > 0){
        int grid_pos = i * row_size + j;
        switch (path_trace[grid_pos]){
            case LEFT_TRANSFER:
                if (i > 0 && i < this->num_positions)
                    init_numbering[i-1] += 1;
                else if (i == 0)
                    numNTermGaps += 1;
                else
                    numCTermGaps += 1;
                j -= 1;
                break;
            case UP_TRANSFER:
                i -= 1;
                break;
            case DIAGONAL_TRANSFER:
                init_numbering[i-1] += 1;
                i -= 1;
                j -= 1;
                break;
            // This should never happen -- if it does, report it as a
            // SERIOUS error.
            default:
                error_code = fatalRuntimeError;
                error_message = this->error_code_to_message[error_code];
                return;
                break;
        }
    }

    // Add gaps at the beginning of the numbering corresponding to the
    // number of N-terminal insertions. This ensures the numbering has
    // the same length as the input sequence.
    // position_key is set to -1 to indicate these positions should not
    // be considered for percent identity calculation.
    for (i=0; i < numNTermGaps; i++){
        final_numbering.push_back("-");
        position_key.push_back(-1);
    }

    
    // Next, let's check to see if there are gaps in positions 1 - 5. The convention
    // in most numbering tools is to push those gaps to the beginning of the
    // sequence. I'm personally not sure I agree with this, but it is what both ANARCI
    // and AbNum do. As a result, we only do this if the user so requests (hence only
    // if compress_initial_gaps is true). At this point, if the sequence has gaps in the
    // first 5 numbered positions, we shuffle them around so the gaps are at the beginning --
    // but again only if user so requested.
    if (this->compress_initial_gaps && query_sequence.length() > 5){
        bool gaps_filled = false;
        while (!gaps_filled){
            gaps_filled = true;
            for (int i=0; i < 5; i++){
                if (init_numbering[i] ==1 && init_numbering[i+1] == 0){
                        init_numbering[i+1] = init_numbering[i];
                        init_numbering[i] = 0;
                        gaps_filled = false;
                }
            }
        }
    }


    // Build vector of IMGT numbers. Unfortunately the IMGT system adds insertion codes
    // forwards then backwards where there is > 1 insertion at a given position,
    // but ONLY if the insertion is at an expected position in the CDRs!
    // Everywhere else, insertions are recognized by just placing in
    // the usual order (?!!!?) This annoying quirk adds some complications.
    // We also create the position_key here which we can quickly use to determine
    // percent identity by referencing this->consensus_map, which indicates
    // what AAs are expected at each position.
    if (this->scheme == "imgt"){
        for (i=0; i < this->num_positions; i++){
            if (init_numbering[i] == 0)
                continue;

            // Highly unlikely given the size of the alphabet we've used that
            // we will ever run into a problem where there are too many insertions,
            // but at least possible.
            if (init_numbering[i] > this->alphabet.size()){
                    error_code = tooManyInsertions;
                    error_message = this->error_code_to_message[error_code];
                    return;
            }
            int ceil_cutpoint, floor_cutpoint;

            // These are the positions for which we need to deploy forward-backward
            // lettering (because the IMGT system is weird, but also the default...)
            // position_key is set to -1 for all insertions to indicate these should
            // not be considered for calculating percent identity.
            switch (i){
                case CDR1_INSERTION_PT:
                case CDR2_INSERTION_PT:
                    ceil_cutpoint = init_numbering[i] / 2;
                    floor_cutpoint = (init_numbering[i] - 1) / 2;
                    for (j=0; j < floor_cutpoint; j++){
                        final_numbering.push_back(std::to_string(i) + this->alphabet[j]);
                        position_key.push_back(-1);
                    }
                    for (j=ceil_cutpoint; j > 0; j--){
                        final_numbering.push_back(std::to_string(i+1) + this->alphabet[j-1]);
                        position_key.push_back(-1);
                    }
    
                    final_numbering.push_back(std::to_string(i+1));
                    position_key.push_back(i);
                    break;

                case CDR3_INSERTION_PT:
                    final_numbering.push_back(std::to_string(i+1));
                    position_key.push_back(i);
                    ceil_cutpoint = init_numbering[i] / 2;
                    floor_cutpoint = (init_numbering[i] - 1) / 2;
                    for (j=0; j < floor_cutpoint; j++){
                        final_numbering.push_back(std::to_string(i+1) + this->alphabet[j]);
                        position_key.push_back(-1);
                    }
                    for (j=ceil_cutpoint; j > 0; j--){
                        final_numbering.push_back(std::to_string(i+2) + this->alphabet[j-1]);
                        position_key.push_back(-1);
                    }
                    break;
                default:
                    final_numbering.push_back(std::to_string(i+1));
                    position_key.push_back(i);
                    for (size_t k=1; k < init_numbering[i]; k++){
                        final_numbering.push_back(std::to_string(i+1) + this->alphabet[k-1]);
                        position_key.push_back(-1);
                    }
                    break;
            }
        }
    }
    // Build vector of numbers for other schemes (much simpler).
    else{
        for (i=0; i < this->num_positions; i++){
            if (init_numbering[i] == 0)
                continue;
            
            if (init_numbering[i] > this->alphabet.size()){
                    error_code = tooManyInsertions;
                    error_message = this->error_code_to_message[error_code];
                    return;
            }
            final_numbering.push_back(std::to_string(i+1));
            position_key.push_back(i);
            for (size_t k=1; k < init_numbering[i]; k++){
                final_numbering.push_back(std::to_string(i+1) + this->alphabet[k-1]);
                position_key.push_back(-1);
            }
        }
    }


    // Add gaps at the end of the numbering corresponding to the
    // number of C-terminal insertions. This ensures the numbering has
    // the same length as the input sequence.
    for (int k=0; k < numCTermGaps; k++){
        final_numbering.push_back("-");
        position_key.push_back(-1);
    }


    // position_key is now the same length as final_numbering and indicates at
    // each position whether that position maps to a standard 
    // position and if so what. Now use this to calculate percent identity,
    // excluding those positions at which the numbering system tolerates any
    // amino acid (i.e. CDRs). At the same time, we can check whether the
    // expected amino acids are present at the highly conserved residues.
    // The consensus map will have only one possible amino acid at those positions.

    int num_required_positions_found = 0;

    for (size_t k=0; k < position_key.size(); k++){
        int scheme_std_position = position_key[k];
        if (scheme_std_position < 0)
            continue;
        
        if (this->consensus_map[scheme_std_position].empty())
            continue;

        if (this->consensus_map[scheme_std_position].find(query_sequence[k]) !=
                this->consensus_map[scheme_std_position].end()){
            percent_identity += 1;
            if(std::binary_search(this->highly_conserved_positions.begin(),
                this->highly_conserved_positions.end(), scheme_std_position)){
                    num_required_positions_found += 1;
            }
        }
    }

    percent_identity /= this->num_restricted_positions;

    error_code = noError;
    // Check to make sure query sequence and numbering are same length. If not,
    // this is a fatal error. (This should be extremely rare -- we have not
    // encountered it so far).
    if (query_sequence.length() != final_numbering.size()){
        error_code = alignmentWrongLength;
        error_message = this->error_code_to_message[error_code];
        return;
    }

    if (num_required_positions_found != 6)
        error_code = unacceptableConservedPositions;

    error_message = this->error_code_to_message[error_code];
}





// Fill in the scoring table created by caller, using the position-specific
// scores for indels and substitutions, and add the appropriate
// pathway traces.
void IGAligner::fill_needle_scoring_table(double *needle_scores, uint8_t *path_trace,
                    int query_seq_len, int row_size, int *encoded_sequence){
    double lscore, uscore, dscore;
    int i, j, grid_pos, diagNeighbor, upperNeighbor;
    auto scoreItr = this->score_array.unchecked<2>();

    // Fill in the first row of the tables. We use a default score here
    // to ensure that insertions in the template only are highly tolerated
    // at the beginning of the sequence.
    needle_scores[0] = 0;
    path_trace[0] = 0;
    for (i=0; i < query_seq_len; i++){
        needle_scores[i+1] = needle_scores[i] + this->terminalTemplateGapPenalty;
        path_trace[i+1] = LEFT_TRANSFER;
    }




    // This first loop goes up to the last row only (the last row)
    // is handled separately since it requires special logic).
    // Handling the last row separately is about 10 - 15% faster than
    // using if else statements to determine when the last row is
    // reached inside this loop.
    for (i=1; i < this->num_positions; i++){

        grid_pos = i * row_size;
        diagNeighbor = (i - 1) * row_size;
        upperNeighbor = diagNeighbor + 1;
        // The first column assigns a modified penalty for gaps so that n-terminal
        // deletions if encountered are accepted.
        needle_scores[grid_pos] = needle_scores[diagNeighbor] + scoreItr(i-1,QUERY_GAP_COLUMN);
        path_trace[grid_pos] = UP_TRANSFER;
        grid_pos++;

        for (j=0; j < (query_seq_len - 1); j++){
            dscore = needle_scores[diagNeighbor] + scoreItr(i-1,encoded_sequence[j]);
            lscore = needle_scores[grid_pos - 1] + scoreItr(i-1,TEMPLATE_GAP_COLUMN);
            uscore = needle_scores[upperNeighbor] + scoreItr(i-1,QUERY_GAP_COLUMN);

            // This is mildly naughty -- we don't consider the possibility
            // of a tie, which could lead to a branched alignment. Realistically,
            // ties are going to be rare and low-impact in this situation because
            // of the way the positions are scored. For now we are defaulting to
            // "match" if there is a tie.
            if (lscore > uscore && lscore > dscore){
                needle_scores[grid_pos] = lscore;
                path_trace[grid_pos] = LEFT_TRANSFER;
            }
            else if (uscore > dscore){
                needle_scores[grid_pos] = uscore;
                path_trace[grid_pos] = UP_TRANSFER;
            }
            else{
                needle_scores[grid_pos] = dscore;
                path_trace[grid_pos] = DIAGONAL_TRANSFER;
            }

            grid_pos++;
            diagNeighbor++;
            upperNeighbor++;
        }

        j = query_seq_len - 1;
        dscore = needle_scores[diagNeighbor] + scoreItr(i-1,encoded_sequence[j]);
        lscore = needle_scores[grid_pos - 1] + scoreItr(i-1,TEMPLATE_GAP_COLUMN);
        // We use a default score for the last column, so that c-terminal
        // deletions if encountered are well-tolerated. We exclude however
        // highly conserved positions, for which a large gap penalty should
        // always be assigned.
        if (std::binary_search(this->highly_conserved_positions.begin(),
                    this->highly_conserved_positions.end(), i))
            uscore = needle_scores[upperNeighbor] + scoreItr(i-1,QUERY_GAP_COLUMN);
        else
            uscore = needle_scores[upperNeighbor] + this->CterminalQueryGapPenalty;

        if (lscore > uscore && lscore > dscore){
            needle_scores[grid_pos] = lscore;
            path_trace[grid_pos] = LEFT_TRANSFER;
        }
        else if (uscore > dscore){
            needle_scores[grid_pos] = uscore;
            path_trace[grid_pos] = UP_TRANSFER;
        }
        else{
            needle_scores[grid_pos] = dscore;
            path_trace[grid_pos] = DIAGONAL_TRANSFER;
        }
    }



    // Now handle the last row.
    i = this->num_positions;
    grid_pos = i * row_size;
    diagNeighbor = (i - 1) * row_size;
    upperNeighbor = diagNeighbor + 1;
    // The first column assigns a low penalty for gaps so that n-terminal
    // deletions if encountered are accepted, UNLESS we are at a highly
    // conserved position.
    needle_scores[grid_pos] = needle_scores[diagNeighbor] + scoreItr(i-1,QUERY_GAP_COLUMN);
    path_trace[grid_pos] = UP_TRANSFER;
    grid_pos++;

    for (j=0; j < (query_seq_len - 1); j++){
        double lscore, uscore, dscore;
        dscore = needle_scores[diagNeighbor] + scoreItr(i-1,encoded_sequence[j]);
        lscore = needle_scores[grid_pos - 1] + this->terminalTemplateGapPenalty;
        uscore = needle_scores[upperNeighbor] + scoreItr(i-1,QUERY_GAP_COLUMN);

        if (lscore > uscore && lscore > dscore){
            needle_scores[grid_pos] = lscore;
            path_trace[grid_pos] = LEFT_TRANSFER;
        }
        else if (uscore > dscore){
            needle_scores[grid_pos] = uscore;
            path_trace[grid_pos] = UP_TRANSFER;
        }
        else{
            needle_scores[grid_pos] = dscore;
            path_trace[grid_pos] = DIAGONAL_TRANSFER;
        }

        grid_pos++;
        diagNeighbor++;
        upperNeighbor++;
    }



    // And, finally, the last column of the last row.
    j = query_seq_len - 1;
    dscore = needle_scores[diagNeighbor] + scoreItr(i-1,encoded_sequence[j]);
    lscore = needle_scores[grid_pos - 1] + scoreItr(i-1,TEMPLATE_GAP_COLUMN);
    uscore = needle_scores[upperNeighbor] + this->CterminalQueryGapPenalty;
    if (lscore > uscore && lscore > dscore){
        needle_scores[grid_pos] = lscore;
        path_trace[grid_pos] = LEFT_TRANSFER;
    }
    else if (uscore > dscore){
        needle_scores[grid_pos] = uscore;
        path_trace[grid_pos] = UP_TRANSFER;
    }
    else{
        needle_scores[grid_pos] = dscore;
        path_trace[grid_pos] = DIAGONAL_TRANSFER;
    }
}





// Wraps IGAligner::core_align_test_only. Unlike align, this function can be
// (and should be) accessed by external Python callers.
std::tuple<std::vector<std::string>, double, std::string,
    std::string> IGAligner::align_test_only(std::string query_sequence,
            py::array_t<double> score_matrix, py::array_t<uint8_t> path_trace){

    std::vector<std::string> final_numbering;
    allowedErrorCodes error_code = noError;

    double percent_identity = this->core_align_test_only(query_sequence,\
                        final_numbering, error_code, score_matrix, path_trace);

    std::string error_message = this->error_code_to_message[error_code];
    return std::tuple<std::vector<std::string>, double, std::string,
                            std::string>{final_numbering,
                            percent_identity, this->chain_name, error_message};
}





// This alignment function is a test-only version of core_align. It ensures that caller
// can access the filled-out scoring matrix, which is useful for testing and diagnostics
// but not really essential for a typical run.
double IGAligner::core_align_test_only(std::string const &query_sequence,
                std::vector<std::string> &final_numbering,
                allowedErrorCodes &error_code,
                py::array_t<double> score_matrix,
                py::array_t<uint8_t> path_traceMat){


    double percent_identity = 0;


    if (query_sequence.length() < MINIMUM_SEQUENCE_LENGTH){
        error_code = invalidSequence;
        return percent_identity;
    }

    std::vector<int> position_key;

    int row_size = query_sequence.length() + 1;

    if (score_matrix.shape(1) != row_size || score_matrix.shape(0) != this->num_positions + 1){
        throw std::runtime_error(std::string("The score_matrix passed to a test function "
                    "does not have the correct shape.")); 
    }
    if (path_traceMat.shape(1) != row_size || path_traceMat.shape(0) != this->num_positions + 1){
        throw std::runtime_error(std::string("The path_trace passed to a test function "
                    "does not have the correct shape.")); 
    }

    auto buffer_info = score_matrix.request();
    double *needle_scores = static_cast<double*>(buffer_info.ptr);
    auto ptrace_buffer_info = path_traceMat.request();
    uint8_t *path_trace = static_cast<uint8_t*>(ptrace_buffer_info.ptr);

    auto encoded_sequence = std::make_unique<int[]>( query_sequence.length() );
    auto init_numbering = std::make_unique<unsigned int[]>( this->num_positions );


    if (!convert_sequence_to_array(encoded_sequence.get(), query_sequence)){
        error_code = invalidSequence;
        return percent_identity;
    }

    // Fill in the scoring table.
    fill_needle_scoring_table(needle_scores, path_trace,
                    query_sequence.length(), row_size, encoded_sequence.get());

    // Set all members of init_numbering array to 0.
    for (int i=0; i < this->num_positions; i++)
        init_numbering[i] = 0;

    // Backtrace the path, and determine how many insertions there
    // were at each position where there is an insertion. Determine
    // how many gaps exist at the beginning and end of the sequence
    // as well (the "N-terminus" and "C-terminus").
    int i = this->num_positions;
    int j = query_sequence.length();
    int numNTermGaps = 0;
    int numCTermGaps = 0;

    while (i > 0 || j > 0){
        int grid_pos = i * row_size + j;
        switch (path_trace[grid_pos]){
            case LEFT_TRANSFER:
                if (i > 0 && i < this->num_positions)
                    init_numbering[i-1] += 1;
                else if (i == 0)
                    numNTermGaps += 1;
                else
                    numCTermGaps += 1;
                j -= 1;
                break;
            case UP_TRANSFER:
                i -= 1;
                break;
            case DIAGONAL_TRANSFER:
                init_numbering[i-1] += 1;
                i -= 1;
                j -= 1;
                break;
            // This should never happen -- if it does, report it as a
            // SERIOUS error.
            default:
                error_code = fatalRuntimeError;
                return percent_identity;
                break;
        }
    }

    // Add gaps at the beginning of the numbering corresponding to the
    // number of N-terminal insertions. This ensures the numbering has
    // the same length as the input sequence.
    // position_key is set to -1 to indicate these positions should not
    // be considered for percent identity calculation.
    for (i=0; i < numNTermGaps; i++){
        final_numbering.push_back("-");
        position_key.push_back(-1);
    }

    
    // Next, let's check to see if there are gaps in positions 1 - 5. The convention
    // in most numbering tools is to push those gaps to the beginning of the
    // sequence. I'm personally not sure I agree with this, but it is what both ANARCI
    // and AbNum do. As a result, we only do this if the user so requests (hence only
    // if compress_initial_gaps is true). At this point, if the sequence has gaps in the
    // first 5 numbered positions, we shuffle them around so the gaps are at the beginning --
    // but again only if user so requested.
    if (this->compress_initial_gaps && query_sequence.length() > 5){
        bool gaps_filled = false;
        while (!gaps_filled){
            gaps_filled = true;
            for (int i=0; i < 5; i++){
                if (init_numbering[i] ==1 && init_numbering[i+1] == 0){
                        init_numbering[i+1] = init_numbering[i];
                        init_numbering[i] = 0;
                        gaps_filled = false;
                }
            }
        }
    }


    // Build vector of IMGT numbers. Unfortunately the IMGT system adds insertion codes
    // forwards then backwards where there is > 1 insertion at a given position,
    // but ONLY if the insertion is at an expected position in the CDRs!
    // Everywhere else, insertions are recognized by just placing in
    // the usual order (?!!!?) This annoying quirk adds some complications.
    // We also create the position_key here which we can quickly use to determine
    // percent identity by referencing this->consensus_map, which indicates
    // what AAs are expected at each position.
    if (this->scheme == "imgt"){
        for (i=0; i < this->num_positions; i++){
            if (init_numbering[i] == 0)
                continue;

            // Highly unlikely given the size of the alphabet we've used that
            // we will ever run into a problem where there are too many insertions,
            // but at least possible.
            if (init_numbering[i] > this->alphabet.size()){
                    error_code = tooManyInsertions;
                    return percent_identity;
            }
            int ceil_cutpoint, floor_cutpoint;

            // These are the positions for which we need to deploy forward-backward
            // lettering (because the IMGT system is weird, but also the default...)
            // position_key is set to -1 for all insertions to indicate these should
            // not be considered for calculating percent identity.
            switch (i){
                case CDR1_INSERTION_PT:
                case CDR2_INSERTION_PT:
                    ceil_cutpoint = init_numbering[i] / 2;
                    floor_cutpoint = (init_numbering[i] - 1) / 2;
                    for (j=0; j < floor_cutpoint; j++){
                        final_numbering.push_back(std::to_string(i) + this->alphabet[j]);
                        position_key.push_back(-1);
                    }
                    for (j=ceil_cutpoint; j > 0; j--){
                        final_numbering.push_back(std::to_string(i+1) + this->alphabet[j-1]);
                        position_key.push_back(-1);
                    }
    
                    final_numbering.push_back(std::to_string(i+1));
                    position_key.push_back(i);
                    break;

                case CDR3_INSERTION_PT:
                    final_numbering.push_back(std::to_string(i+1));
                    position_key.push_back(i);
                    ceil_cutpoint = init_numbering[i] / 2;
                    floor_cutpoint = (init_numbering[i] - 1) / 2;
                    for (j=0; j < floor_cutpoint; j++){
                        final_numbering.push_back(std::to_string(i+1) + this->alphabet[j]);
                        position_key.push_back(-1);
                    }
                    for (j=ceil_cutpoint; j > 0; j--){
                        final_numbering.push_back(std::to_string(i+2) + this->alphabet[j-1]);
                        position_key.push_back(-1);
                    }
                    break;
                default:
                    final_numbering.push_back(std::to_string(i+1));
                    position_key.push_back(i);
                    for (size_t k=1; k < init_numbering[i]; k++){
                        final_numbering.push_back(std::to_string(i+1) + this->alphabet[k-1]);
                        position_key.push_back(-1);
                    }
                    break;
            }
        }
    }
    // Build vector of numbers for other schemes (much simpler).
    else{
        for (i=0; i < this->num_positions; i++){
            if (init_numbering[i] == 0)
                continue;
            
            if (init_numbering[i] > this->alphabet.size()){
                    error_code = tooManyInsertions;
                    return percent_identity;
            }
            final_numbering.push_back(std::to_string(i+1));
            position_key.push_back(i);
            for (size_t k=1; k < init_numbering[i]; k++){
                final_numbering.push_back(std::to_string(i+1) + this->alphabet[k-1]);
                position_key.push_back(-1);
            }
        }
    }


    // Add gaps at the end of the numbering corresponding to the
    // number of C-terminal insertions. This ensures the numbering has
    // the same length as the input sequence.
    for (int k=0; k < numCTermGaps; k++){
        final_numbering.push_back("-");
        position_key.push_back(-1);
    }


    // position_key is now the same length as final_numbering and indicates at
    // each position whether that position maps to a standard 
    // position and if so what. Now use this to calculate percent identity,
    // excluding those positions at which the numbering system tolerates any
    // amino acid (i.e. CDRs). At the same time, we can check whether the
    // expected amino acids are present at the highly conserved residues.
    // The consensus map will have only one possible amino acid at those positions.

    int num_required_positions_found = 0;

    for (size_t k=0; k < position_key.size(); k++){
        int scheme_std_position = position_key[k];
        if (scheme_std_position < 0)
            continue;
        
        if (this->consensus_map[scheme_std_position].empty())
            continue;

        if (this->consensus_map[scheme_std_position].find(query_sequence[k]) !=
                this->consensus_map[scheme_std_position].end()){
            percent_identity += 1;
            if(std::binary_search(this->highly_conserved_positions.begin(),
                this->highly_conserved_positions.end(), scheme_std_position)){
                    num_required_positions_found += 1;
            }
        }
    }

    percent_identity /= this->num_restricted_positions;

    error_code = noError;
    // Check to make sure query sequence and numbering are same length. If not,
    // this is a fatal error. (This should be extremely rare -- we have not
    // encountered it so far).
    if (query_sequence.length() != final_numbering.size()){
        error_code = alignmentWrongLength;
        return percent_identity;
    }

    if (num_required_positions_found != 6)
        error_code = unacceptableConservedPositions;

    return percent_identity;
}
