#include "ig_aligner.h"
#include "utilities.h"


IGAligner::IGAligner(
                 py::array_t<double> scoreArray,
                 std::vector<std::vector<std::string>> consensus,
                 std::string chainName,
                 std::string scheme,
                 double terminalTemplateGapPenalty,
                 double CterminalQueryGapPenalty,
                 bool compressInitialGaps
):
    scoreArray(scoreArray),
    chainName(chainName),
    scheme(scheme),
    terminalTemplateGapPenalty(terminalTemplateGapPenalty),
    CterminalQueryGapPenalty(CterminalQueryGapPenalty),
    compressInitialGaps(compressInitialGaps)
{
    py::buffer_info info = scoreArray.request();
    // Note that exceptions thrown here are go back to Python via
    // PyBind as long as this constructor is used within the wrapper.
    if (info.shape.size() != 2){
        throw std::runtime_error(std::string("The scoreArray passed to "
                "IGAligner must be a 2d array"));
    }
    if (info.shape[1] != NUM_AAS){
        throw std::runtime_error(std::string("The scoreArray passed to "
                "IGAligner must have 22 columns (1 per AA plus 2 for gap penalties)"));
    }

    // Check that the number of positions in the scoring and consensus is as expected,
    // and set the list of highly conserved positions according to the selected scheme.
    if (scheme == "imgt"){
        this->highlyConservedPositions = {HIGHLY_CONSERVED_IMGT_1, HIGHLY_CONSERVED_IMGT_2,
                                        HIGHLY_CONSERVED_IMGT_3, HIGHLY_CONSERVED_IMGT_4,
                                        HIGHLY_CONSERVED_IMGT_5, HIGHLY_CONSERVED_IMGT_6};
        if (info.shape[0] != NUM_HEAVY_IMGT_POSITIONS && info.shape[0] !=
                NUM_LIGHT_IMGT_POSITIONS){
            throw std::runtime_error(std::string("The scoreArray passed to "
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
    else if (scheme == "martin" || scheme == "kabat"){
        if (chainName == "L" || chainName == "K"){
            this->highlyConservedPositions = {HIGHLY_CONSERVED_KABAT_LIGHT_1,
                    HIGHLY_CONSERVED_KABAT_LIGHT_2, HIGHLY_CONSERVED_KABAT_LIGHT_3,
                    HIGHLY_CONSERVED_KABAT_LIGHT_4, HIGHLY_CONSERVED_KABAT_LIGHT_5,
                    HIGHLY_CONSERVED_KABAT_LIGHT_6};
        }
        else if (chainName == "H"){
            this->highlyConservedPositions = {HIGHLY_CONSERVED_KABAT_HEAVY_1,
                    HIGHLY_CONSERVED_KABAT_HEAVY_2, HIGHLY_CONSERVED_KABAT_HEAVY_3,
                    HIGHLY_CONSERVED_KABAT_HEAVY_4, HIGHLY_CONSERVED_KABAT_HEAVY_5,
                    HIGHLY_CONSERVED_KABAT_HEAVY_6};
        }
        else{
            throw std::runtime_error(std::string("TMartin, Kabat currently "
                        "support only H, K, L chains."));
        }
        if (info.shape[0] != NUM_HEAVY_MARTIN_KABAT_POSITIONS && info.shape[0] !=
                NUM_LIGHT_MARTIN_KABAT_POSITIONS){
            throw std::runtime_error(std::string("The scoreArray passed to "
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
                    "schemes 'martin', 'kabat', 'imgt'."));
    }

    numPositions = info.shape[0];
    // Update the consensus map to indicate which letters are expected at
    // which positions. An empty set at a given position means that
    // ANY letter at that position is considered acceptable. These are
    // excluded from percent identity calculation. numRestrictedPositions
    // indicates how many are INCLUDED in percent identity.
    numRestrictedPositions = 0;

    for (size_t i = 0; i < consensus.size(); i++){
        consensusMap.push_back(std::set<char> {});
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




    // Set up a list of breakpoints that mark the dividing lines between
    // framework and CDR regions according to the selected chain type and
    // numbering scheme.
    if (scheme == "imgt"){
        this->cdrBreakpoints = {0, IMGT_CDR_BREAKPOINT_1, IMGT_CDR_BREAKPOINT_2,
                        IMGT_CDR_BREAKPOINT_3, IMGT_CDR_BREAKPOINT_4,
                        IMGT_CDR_BREAKPOINT_5, IMGT_CDR_BREAKPOINT_6,
                        this->numPositions};
    }
    else if (scheme == "kabat"){
        if (chainName == "L" || chainName == "K"){
            this->cdrBreakpoints = {0, KABAT_LIGHT_CDR_BREAKPOINT_1,
                        KABAT_LIGHT_CDR_BREAKPOINT_2, KABAT_LIGHT_CDR_BREAKPOINT_3,
                        KABAT_LIGHT_CDR_BREAKPOINT_4, KABAT_LIGHT_CDR_BREAKPOINT_5,
                        KABAT_LIGHT_CDR_BREAKPOINT_6,
                        this->numPositions};
        }
        else if (chainName == "H"){
            this->cdrBreakpoints = {0, KABAT_HEAVY_CDR_BREAKPOINT_1,
                        KABAT_HEAVY_CDR_BREAKPOINT_2, KABAT_HEAVY_CDR_BREAKPOINT_3,
                        KABAT_HEAVY_CDR_BREAKPOINT_4, KABAT_HEAVY_CDR_BREAKPOINT_5,
                        KABAT_HEAVY_CDR_BREAKPOINT_6,
                        this->numPositions};
        }
    }
    else if (scheme == "martin"){
        if (chainName == "L" || chainName == "K"){
            this->cdrBreakpoints = {0, MARTIN_LIGHT_CDR_BREAKPOINT_1,
                        MARTIN_LIGHT_CDR_BREAKPOINT_2, MARTIN_LIGHT_CDR_BREAKPOINT_3,
                        MARTIN_LIGHT_CDR_BREAKPOINT_4, MARTIN_LIGHT_CDR_BREAKPOINT_5,
                        MARTIN_LIGHT_CDR_BREAKPOINT_6,
                        this->numPositions};
        }
        else if (chainName == "H"){
            this->cdrBreakpoints = {0, MARTIN_HEAVY_CDR_BREAKPOINT_1,
                        MARTIN_HEAVY_CDR_BREAKPOINT_2, MARTIN_HEAVY_CDR_BREAKPOINT_3,
                        MARTIN_HEAVY_CDR_BREAKPOINT_4, MARTIN_HEAVY_CDR_BREAKPOINT_5,
                        MARTIN_HEAVY_CDR_BREAKPOINT_6,
                        this->numPositions};
        }
    }
    else{
        throw std::runtime_error(std::string("Currently IGAligner only recognizes "
                    "schemes 'martin', 'kabat', 'imgt'."));
    }
}






std::tuple<std::vector<std::string>, double, std::string,
    std::string, std::vector<std::string>> IGAligner::align(std::string query_sequence,
                bool retrieve_cdr_labeling){


    std::vector<std::string> finalNumbering;
    std::vector<std::string> cdrLabeling;
    double percentIdentity = 0;
    std::string errorMessage;

    // A list of allowed error codes. These will be mapped to strings that explain
    // in more detail by the IGAligner class.
    enum allowedErrorCodes {noError = 0, invalidSequence = 1, fatalRuntimeError = 2,
            tooManyInsertions  = 3, alignmentWrongLength = 4,
            unacceptableConservedPositions = 5};
    allowedErrorCodes errorCode;

    if (query_sequence.length() < MINIMUM_SEQUENCE_LENGTH){
        errorCode = invalidSequence;
        errorMessage = this->errorCodeToMessage[errorCode];
        return std::tuple<std::vector<std::string>, double, std::string,
                        std::string, std::vector<std::string>>{finalNumbering,
                            percentIdentity, this->chainName, errorMessage,
                            cdrLabeling};
    }

    std::vector<int> positionKey;

    int rowSize = query_sequence.length() + 1;
    int numElements = rowSize * (this->numPositions + 1);

    double *needleScores = new double[ numElements ];
    int *queryAsIdx = new int[ query_sequence.length() ];
    int *pathTrace = new int[ numElements ];
    int unsigned *initNumbering = new unsigned int[ this->numPositions ];


    if (!convert_sequence_to_array(queryAsIdx, query_sequence)){
        delete[] queryAsIdx;
        delete[] needleScores;
        delete[] pathTrace;
        delete[] initNumbering;
        errorCode = invalidSequence;
        errorMessage = this->errorCodeToMessage[errorCode];
        return std::tuple<std::vector<std::string>, double, std::string,
                        std::string, std::vector<std::string>>{finalNumbering,
                            percentIdentity, this->chainName, errorMessage,
                            cdrLabeling};
    }

    // Fill in the scoring table.
    fillNeedleScoringTable(needleScores, pathTrace,
                    query_sequence.length(), rowSize, queryAsIdx);

    // Set all members of initNumbering array to 0.
    for (int i=0; i < this->numPositions; i++)
        initNumbering[i] = 0;

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
                errorCode = fatalRuntimeError;
                errorMessage = this->errorCodeToMessage[errorCode];
                return std::tuple<std::vector<std::string>, double, std::string,
                        std::string, std::vector<std::string>>{finalNumbering,
                            percentIdentity, this->chainName, errorMessage,
                            cdrLabeling};
                break;
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

    
    // Next, let's check to see if there are gaps in positions 1 - 5. The convention
    // in most numbering tools is to push those gaps to the beginning of the
    // sequence. I'm personally not sure I agree with this, but it is what both ANARCI
    // and AbNum do. As a result, we only do this if the user so requests (hence only
    // if compressInitialGaps is true). At this point, if the sequence has gaps in the
    // first 5 numbered positions, we shuffle them around so the gaps are at the beginning --
    // but again only if user so requested.
    if (this->compressInitialGaps && query_sequence.length() > 5){
        bool gapsFilled = false;
        while (!gapsFilled){
            gapsFilled = true;
            for (int i=0; i < 5; i++){
                if (initNumbering[i] ==1 && initNumbering[i+1] == 0){
                        initNumbering[i+1] = initNumbering[i];
                        initNumbering[i] = 0;
                        gapsFilled = false;
                }
            }
        }
    }


    // Build vector of IMGT numbers. Unfortunately the IMGT system adds insertion codes
    // forwards then backwards where there is > 1 insertion at a given position,
    // but ONLY if the insertion is at an expected position in the CDRs!
    // Everywhere else, insertions are recognized by just placing in
    // the usual order (?!!!?) This annoying quirk adds some complications.
    // We also create the positionKey here which we can quickly use to determine
    // percent identity by referencing this->consensusMap, which indicates
    // what AAs are expected at each position.
    if (this->scheme == "imgt"){
        for (i=0; i < this->numPositions; i++){
            if (initNumbering[i] == 0)
                continue;

            // Highly unlikely given the size of the alphabet we've used that
            // we will ever run into a problem where there are too many insertions,
            // but at least possible.
            if (initNumbering[i] > this->alphabet.size()){
                    delete[] queryAsIdx;
                    delete[] needleScores;
                    delete[] pathTrace;
                    delete[] initNumbering;
                    errorCode = tooManyInsertions;
                    errorMessage = this->errorCodeToMessage[errorCode];
                    return std::tuple<std::vector<std::string>, double, std::string,
                            std::string, std::vector<std::string>>{finalNumbering,
                            percentIdentity, this->chainName, errorMessage,
                            cdrLabeling};
            }
            int ceil_cutpoint, floor_cutpoint;

            // These are the positions for which we need to deploy forward-backward
            // lettering (because the IMGT system is weird, but also the default...)
            // positionKey is set to -1 for all insertions to indicate these should
            // not be considered for calculating percent identity.
            switch (i){
                case CDR1_INSERTION_PT:
                case CDR2_INSERTION_PT:
                    ceil_cutpoint = initNumbering[i] / 2;
                    floor_cutpoint = (initNumbering[i] - 1) / 2;
                    for (j=0; j < floor_cutpoint; j++){
                        finalNumbering.push_back(std::to_string(i) + this->alphabet[j]);
                        positionKey.push_back(-1);
                    }
                    for (j=ceil_cutpoint; j > 0; j--){
                        finalNumbering.push_back(std::to_string(i+1) + this->alphabet[j-1]);
                        positionKey.push_back(-1);
                    }
    
                    finalNumbering.push_back(std::to_string(i+1));
                    positionKey.push_back(i);
                    break;

                case CDR3_INSERTION_PT:
                    finalNumbering.push_back(std::to_string(i+1));
                    positionKey.push_back(i);
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
                    finalNumbering.push_back(std::to_string(i+1));
                    positionKey.push_back(i);
                    for (size_t k=1; k < initNumbering[i]; k++){
                        finalNumbering.push_back(std::to_string(i+1) + this->alphabet[k-1]);
                        positionKey.push_back(-1);
                    }
                    break;
            }
        }
    }
    // Build vector of numbers for other schemes (much simpler).
    else{
        for (i=0; i < this->numPositions; i++){
            if (initNumbering[i] == 0)
                continue;
            
            if (initNumbering[i] > this->alphabet.size()){
                    delete[] queryAsIdx;
                    delete[] needleScores;
                    delete[] pathTrace;
                    delete[] initNumbering;
                    errorCode = tooManyInsertions;
                    errorMessage = this->errorCodeToMessage[errorCode];
                    return std::tuple<std::vector<std::string>, double, std::string,
                            std::string, std::vector<std::string>>{finalNumbering,
                            percentIdentity, this->chainName, errorMessage,
                            cdrLabeling};
            }
            finalNumbering.push_back(std::to_string(i+1));
            positionKey.push_back(i);
            for (size_t k=1; k < initNumbering[i]; k++){
                finalNumbering.push_back(std::to_string(i+1) + this->alphabet[k-1]);
                positionKey.push_back(-1);
            }
        }
    }


    // Add gaps at the end of the numbering corresponding to the
    // number of C-terminal insertions. This ensures the numbering has
    // the same length as the input sequence.
    for (int k=0; k < numCTermGaps; k++){
        finalNumbering.push_back("-");
        positionKey.push_back(-1);
    }

    // If the user wants to retrieve cdr labeling, we now populate a second vector
    // which indicates what framework or CDR region each position belongs to. Note
    // that this uses whatever numbering scheme is used to do the numbering. In some
    // cases this may be undesired (e.g. user may want to extract Kabat CDRs from
    // an IMGT-numbered sequence), but in the most common use-case they will not
    // need to do this.
    if (retrieve_cdr_labeling){
        for (int k=0; k < numNTermGaps; k++)
            cdrLabeling.push_back("-");

        for (int k=0; k < 7; k++){
            for (int m=this->cdrBreakpoints[k]; m < this->cdrBreakpoints[k+1]; m++){
                for (size_t p=0; p < initNumbering[m]; p++)
                    cdrLabeling.push_back(this->cdrRegionLabels[k]);
            }
        }

        for (int k=0; k < numCTermGaps; k++)
            cdrLabeling.push_back("-");
    }

    // positionKey is now the same length as finalNumbering and indicates at
    // each position whether that position maps to a standard 
    // position and if so what. Now use this to calculate percent identity,
    // excluding those positions at which the numbering system tolerates any
    // amino acid (i.e. CDRs). At the same time, we can check whether the
    // expected amino acids are present at the highly conserved residues.
    // The consensus map will have only one possible amino acid at those positions.

    int numRequiredPositionsFound = 0;

    for (size_t k=0; k < positionKey.size(); k++){
        int schemeStdPosition = positionKey[k];
        if (schemeStdPosition < 0)
            continue;
        
        if (this->consensusMap[schemeStdPosition].empty())
            continue;

        if (this->consensusMap[schemeStdPosition].find(query_sequence[k]) !=
                this->consensusMap[schemeStdPosition].end()){
            percentIdentity += 1;
            if(std::binary_search(this->highlyConservedPositions.begin(),
                this->highlyConservedPositions.end(), schemeStdPosition)){
                    numRequiredPositionsFound += 1;
            }
        }
    }

    percentIdentity /= this->numRestrictedPositions;

    delete[] queryAsIdx;
    delete[] needleScores;
    delete[] pathTrace;
    delete[] initNumbering;

    errorCode = noError;
    // Check to make sure query sequence and numbering are same length. If not,
    // this is a fatal error. (This should be extremely rare -- we have not
    // encountered it so far).
    if (query_sequence.length() != finalNumbering.size()){
        errorCode = alignmentWrongLength;
        errorMessage = this->errorCodeToMessage[errorCode];
        return std::tuple<std::vector<std::string>, double, std::string,
                            std::string, std::vector<std::string>>{finalNumbering,
                            percentIdentity, this->chainName, errorMessage,
                            cdrLabeling};
    }
    // If the user requested cdr labeling, check to make sure the cdr labeling is
    // the same length as the query sequence (again, anything else would be
    // very unusual -- have not yet encountered such an issue).
    if (retrieve_cdr_labeling){
        if (query_sequence.length() != cdrLabeling.size()){
            errorCode = alignmentWrongLength;
            errorMessage = this->errorCodeToMessage[errorCode];
            return std::tuple<std::vector<std::string>, double, std::string,
                            std::string, std::vector<std::string>>{finalNumbering,
                            percentIdentity, this->chainName, errorMessage,
                            cdrLabeling};
        }
    }

    if (numRequiredPositionsFound != 6)
        errorCode = unacceptableConservedPositions;

    errorMessage = this->errorCodeToMessage[errorCode];
    return std::tuple<std::vector<std::string>, double, std::string,
                            std::string, std::vector<std::string>>{finalNumbering,
                            percentIdentity, this->chainName, errorMessage,
                            cdrLabeling};
}








// Fill in the scoring table created by caller, using the position-specific
// scores for indels and substitutions, and add the appropriate
// pathway traces.
void IGAligner::fillNeedleScoringTable(double *needleScores, int *pathTrace,
                    int querySeqLen, int rowSize, int *queryAsIdx){
    double lscore, uscore, dscore;
    int i, j, gridPos, diagNeighbor, upperNeighbor;
    auto scoreItr = this->scoreArray.unchecked<2>();

    // Fill in the first row of the tables. We use a default score here
    // to ensure that insertions in the template only are highly tolerated
    // at the beginning of the sequence.
    needleScores[0] = 0;
    pathTrace[0] = 0;
    for (i=0; i < querySeqLen; i++){
        needleScores[i+1] = needleScores[i] + this->terminalTemplateGapPenalty;
        pathTrace[i+1] = LEFT_TRANSFER;
    }




    // This first loop goes up to the last row only (the last row)
    // is handled separately since it requires special logic).
    // Handling the last row separately is about 10 - 15% faster than
    // using if else statements to determine when the last row is
    // reached inside this loop.
    for (i=1; i < this->numPositions; i++){

        gridPos = i * rowSize;
        diagNeighbor = (i - 1) * rowSize;
        upperNeighbor = diagNeighbor + 1;
        // The first column assigns a modified penalty for gaps so that n-terminal
        // deletions if encountered are accepted.
        needleScores[gridPos] = needleScores[diagNeighbor] + scoreItr(i-1,QUERY_GAP_COLUMN);
        pathTrace[gridPos] = UP_TRANSFER;
        gridPos++;

        for (j=0; j < (querySeqLen - 1); j++){
            dscore = needleScores[diagNeighbor] + scoreItr(i-1,queryAsIdx[j]);
            lscore = needleScores[gridPos - 1] + scoreItr(i-1,TEMPLATE_GAP_COLUMN);
            uscore = needleScores[upperNeighbor] + scoreItr(i-1,QUERY_GAP_COLUMN);

            // This is mildly naughty -- we don't consider the possibility
            // of a tie, which could lead to a branched alignment. Realistically,
            // ties are going to be rare and low-impact in this situation because
            // of the way the positions are scored. For now we are defaulting to
            // "match" if there is a tie.
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

        j = querySeqLen - 1;
        dscore = needleScores[diagNeighbor] + scoreItr(i-1,queryAsIdx[j]);
        lscore = needleScores[gridPos - 1] + scoreItr(i-1,TEMPLATE_GAP_COLUMN);
        // We use a default score for the last column, so that c-terminal
        // deletions if encountered are well-tolerated. We exclude however
        // highly conserved positions, for which a large gap penalty should
        // always be assigned.
        if (std::binary_search(this->highlyConservedPositions.begin(),
                    this->highlyConservedPositions.end(), i))
            uscore = needleScores[upperNeighbor] + scoreItr(i-1,QUERY_GAP_COLUMN);
        else
            uscore = needleScores[upperNeighbor] + this->CterminalQueryGapPenalty;

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
    }



    // Now handle the last row.
    i = this->numPositions;
    gridPos = i * rowSize;
    diagNeighbor = (i - 1) * rowSize;
    upperNeighbor = diagNeighbor + 1;
    // The first column assigns a low penalty for gaps so that n-terminal
    // deletions if encountered are accepted, UNLESS we are at a highly
    // conserved position.
    needleScores[gridPos] = needleScores[diagNeighbor] + scoreItr(i-1,QUERY_GAP_COLUMN);
    pathTrace[gridPos] = UP_TRANSFER;
    gridPos++;

    for (j=0; j < (querySeqLen - 1); j++){
        double lscore, uscore, dscore;
        dscore = needleScores[diagNeighbor] + scoreItr(i-1,queryAsIdx[j]);
        lscore = needleScores[gridPos - 1] + this->terminalTemplateGapPenalty;
        uscore = needleScores[upperNeighbor] + scoreItr(i-1,QUERY_GAP_COLUMN);

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



    // And, finally, the last column of the last row.
    j = querySeqLen - 1;
    dscore = needleScores[diagNeighbor] + scoreItr(i-1,queryAsIdx[j]);
    lscore = needleScores[gridPos - 1] + scoreItr(i-1,TEMPLATE_GAP_COLUMN);
    uscore = needleScores[upperNeighbor] + this->CterminalQueryGapPenalty;
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
}
