#include "vj_match_counter.h"
#include "utilities.h"


VJMatchCounter::VJMatchCounter(
                 std::vector<std::string> geneSeqs,
                 std::vector<std::string> geneNames
):
    geneSeqs(geneSeqs),
    geneNames(geneNames)
{
    // Note that exceptions thrown here go back to Python via
    // PyBind as long as this constructor is used within the wrapper.
    if (geneSeqs.size() != geneNames.size()){
        throw std::runtime_error(std::string("The number of sequences passed "
                    "to VJMatchCounter must match the number of names."));
    }
    for (size_t i=0; i < geneSeqs.size(); i++){
        if (geneSeqs[i].length() != REQUIRED_SEQUENCE_LENGTH){
            throw std::runtime_error(std::string("All sequences passed to VJMatchCounter "
                        "must have the correct length."));
        }
        this->namesToPositions[geneNames[i]] = i;
    }
}



std::tuple<std::string, double> VJMatchCounter::vjMatch(std::string query_sequence){

    // Note that exceptions thrown here go back to Python via
    // PyBind as long as this constructor is used within the wrapper.
    if (query_sequence.length() != REQUIRED_SEQUENCE_LENGTH){
        throw std::runtime_error(std::string("All sequences passed to VJMatchCounter "
                        "must have the correct length."));
    }

    int matchingPositions, nonzeroPositions, closestId = 0;
    double bestIdentity = 0, currentIdentity = 0;
    std::string assignedVJ;

    for (size_t i=0; i < this->geneSeqs.size(); i++){
        matchingPositions = 0;
        nonzeroPositions = 0;

        for (std::string::size_type j=0; j < REQUIRED_SEQUENCE_LENGTH; j++){
            if (this->geneSeqs[i][j] == '-')
                continue;
            nonzeroPositions += 1;
            if (this->geneSeqs[i][j] == query_sequence[j])
                matchingPositions += 1;
        }
        if (nonzeroPositions == 0)
            nonzeroPositions = 1;
        currentIdentity = (static_cast<double>(matchingPositions)) / 
            (static_cast<double>(nonzeroPositions));

        if (currentIdentity > bestIdentity){
            closestId = i;
            bestIdentity = currentIdentity;
        }
    }

    assignedVJ = this->geneNames[closestId];
    return std::tuple<std::string, double>{assignedVJ, bestIdentity};
}


std::string VJMatchCounter::findVJSequenceByName(std::string query_name){

    size_t matchID;
    std::string output;

    if (this->namesToPositions.find(query_name) != this->namesToPositions.end()){
        matchID = this->namesToPositions[query_name];
        output = this->geneSeqs[matchID];
    }
    else{
        output = "";
    }
    return output;
}


// This function just returns the vectors stored by the object,
// which pybind stl can convert to Python lists. This is used
// only for testing.
std::tuple<std::vector<std::string>,
   std::vector<std::string>>  VJMatchCounter::getSeqLists(){

    return std::tuple<std::vector<std::string>,
            std::vector<std::string>>{this->geneSeqs, this->geneNames};
}
