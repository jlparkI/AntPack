/* Contains the implementation of the generic chain annotator tool for
 * sequences containing an unknown number and type of chains. This is
 * the most generic annotator but also the slowest since it can afford
 * to make fewer assumptions. */
#include "chain_annotator.h"





ChainAnnotatorCpp::ChainAnnotatorCpp(
        std::string scheme, std::string consensus_filepath
):
    AnnotatorBaseClassCpp(scheme, consensus_filepath),
    scheme(scheme)
{
    std::vector<std::string> chains = {"H", "K", "L"};
    this->analyzer = std::make_unique<SingleChainAnnotatorCpp>(chains,
            scheme, false, consensus_filepath);
}



// ChainAnnotatorCpp function which numbers an input sequence.
// A list of alignments is returned, one for each chain found.
std::vector<std::tuple<std::vector<std::string>, double, std::string, std::string>>
    ChainAnnotatorCpp::analyze_seq(std::string sequence, double pidentity_threshold){

    std::vector<std::pair<size_t,size_t>> subregions;
    this->split_sequence_into_subregions(subregions, sequence, 100);

    // Now align all the identified segments and return these as a list.
    std::vector<std::tuple<std::vector<std::string>, double, std::string,
        std::string>> output_results;

    for (auto &subregion : subregions){
        std::string subsequence = sequence.substr(subregion.first,
                subregion.second - subregion.first);
        std::tuple<std::vector<std::string>, double, std::string, std::string> result;
        std::vector<std::string> output_numbering(subregion.first, "-");

        result = this->analyzer->analyze_seq(subsequence);
        for (auto &token : std::get<0>(result))
            output_numbering.push_back(token);

        if (std::get<1>(result) < pidentity_threshold)
            continue;

        for (size_t i=subregion.second; i < sequence.length(); i++)
            output_numbering.push_back("-");

        std::get<0>(result) = output_numbering;
        output_results.push_back(result);
    }

    return output_results;
}
