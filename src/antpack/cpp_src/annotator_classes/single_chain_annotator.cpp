#include "single_chain_annotator.h"


// Important scheme-specific defaults. Do not change unless
// necessary.
#define IMGT_DEFAULT_TERMINAL_TEMPLATE_GAP_PENALTY -1
#define IMGT_DEFAULT_C_TERMINAL_QUERY_GAP_PENALTY -1
#define AHO_DEFAULT_TERMINAL_TEMPLATE_GAP_PENALTY -1
#define AHO_DEFAULT_C_TERMINAL_QUERY_GAP_PENALTY -1
#define MARTIN_DEFAULT_TERMINAL_TEMPLATE_GAP_PENALTY -1
#define MARTIN_DEFAULT_C_TERMINAL_QUERY_GAP_PENALTY -1
#define KABAT_DEFAULT_TERMINAL_TEMPLATE_GAP_PENALTY -1
#define KABAT_DEFAULT_C_TERMINAL_QUERY_GAP_PENALTY -1


// Special error code for cases where a rare but specific type
// of alignment error may have occurred.
#define POSSIBLE_FWGXG_ERROR_ON_ALIGNMENT 2



SingleChainAnnotatorCpp::SingleChainAnnotatorCpp(
        std::vector<std::string> chains,
                std::string scheme, bool compress_init_gaps,
                std::string consensus_filepath
):
    AnnotatorBaseClassCpp(scheme, consensus_filepath),
    chains(chains),
    scheme(scheme),
    compress_init_gaps(compress_init_gaps)
{
    // Note that exceptions thrown here go back to Python via
    // PyBind as long as this constructor is used within the wrapper.
    if (chains.size() < 1 || chains.size() > 3){
        throw std::runtime_error(std::string("There must be at least one chain specified, "
                    "and no more than 3."));
    }


    for (auto & chain : this->chains){
        if (chain != "H" && chain != "K" && chain != "L")
            throw std::runtime_error(std::string("Valid chains must be one of 'H', 'K', 'L'."));
        else{

            if (scheme == "imgt"){
                this->scoring_tools.push_back(std::make_unique<IGAligner>(consensus_filepath,
                            chain, scheme,
                            IMGT_DEFAULT_TERMINAL_TEMPLATE_GAP_PENALTY,
                            IMGT_DEFAULT_C_TERMINAL_QUERY_GAP_PENALTY,
                            compress_init_gaps)
                        );
            }
            else if (scheme == "aho"){
                this->scoring_tools.push_back(std::make_unique<IGAligner>(consensus_filepath,
                            chain, scheme,
                            AHO_DEFAULT_TERMINAL_TEMPLATE_GAP_PENALTY,
                            AHO_DEFAULT_C_TERMINAL_QUERY_GAP_PENALTY,
                            compress_init_gaps)
                        );
            }
            else if (scheme == "martin"){
                this->scoring_tools.push_back(std::make_unique<IGAligner>(consensus_filepath,
                            chain, scheme,
                            MARTIN_DEFAULT_TERMINAL_TEMPLATE_GAP_PENALTY,
                            MARTIN_DEFAULT_C_TERMINAL_QUERY_GAP_PENALTY,
                            compress_init_gaps)
                        );
            }
            else if (scheme == "kabat"){
                this->scoring_tools.push_back(std::make_unique<IGAligner>(consensus_filepath,
                            chain, scheme,
                            KABAT_DEFAULT_TERMINAL_TEMPLATE_GAP_PENALTY,
                            KABAT_DEFAULT_C_TERMINAL_QUERY_GAP_PENALTY,
                            compress_init_gaps)
                        );
            }
            else{
                throw std::runtime_error(std::string("Invalid scheme specified. Please use "
                            "one of 'imgt', 'kabat', 'martin', 'aho'."));
            }
        }
    }
}



// SingleChainAnnotatorCpp function which numbers a single input sequence.
std::tuple<std::vector<std::string>, double, std::string,
            std::string> SingleChainAnnotatorCpp::analyze_seq(std::string sequence){
   
    double best_identity = -1;
    std::vector<std::string> emptyNumbering; 
    std::tuple<std::vector<std::string>, double, std::string,
                        std::string> best_result{ emptyNumbering,
                            0, "", "Invalid sequence supplied -- nonstandard AAs"};


    if (!validate_sequence(sequence))
        return best_result;

    int err_code = this->align_input_subregion(best_result, best_identity,
                sequence);
    if (err_code !=  POSSIBLE_FWGXG_ERROR_ON_ALIGNMENT)
        return best_result;

    // It can (rarely) happen that we have a chain (usually light) where the expected FGxG
    // motif in the J-gene is altered AND there is significant additional sequence beyond
    // the end of the J-gene, e.g. in a paired chain AND this results in an alignment
    // error. In the event this unusual combination of circumstances occurs, it
    // usually manifests in the form of an unusually long CDR that exceeds the maximum
    // number of allowed insertion codes.

    // If this may have occurred, subdivide the sequence into regions just as a generic
    // ChainAnnotator would, align each region and return the one with the highest
    // percent identity. This will work well unless the user has accidentally supplied
    // a paired chain, in which case they should have used either PairedChainAnnotator
    // or ChainAnnotator -- this tool is called SingleChainAnnotator because it extracts
    // a single chain (even if multiple are present).
    std::vector<std::pair<size_t,size_t>> subregions;
    this->split_sequence_into_subregions(subregions, sequence);

    // Now align all the identified segments and return the best match.
    best_identity = -1;

    for (auto &subregion : subregions){
        double percent_identity = -1;
        std::string subsequence = sequence.substr(subregion.first,
                subregion.second - subregion.first);
        std::tuple<std::vector<std::string>, double, std::string, std::string> result;
        std::vector<std::string> output_numbering(subregion.first, "-");

        err_code = this->align_input_subregion(result, percent_identity,
                subsequence);
        if (percent_identity < best_identity || err_code != VALID_SEQUENCE)
            continue;

        best_identity = percent_identity;

        for (auto &token : std::get<0>(result))
            output_numbering.push_back(token);


        for (size_t i=subregion.second; i < sequence.length(); i++)
            output_numbering.push_back("-");

        std::get<0>(result) = std::move(output_numbering);
        best_result = std::move(result);
    }

    return best_result;
}




// Aligns the input sequence, which may be either the full sequence or a subregion of it.
int SingleChainAnnotatorCpp::align_input_subregion(std::tuple<std::vector<std::string>, double,
                std::string, std::string> &best_result, double &best_identity,
                std::string &query_sequence){
    auto queryAsIdx = std::make_unique<int[]>( query_sequence.length() );
    bool fwgxg_error = false;

    if (!convert_sequence_to_array(queryAsIdx.get(), query_sequence))
        return INVALID_SEQUENCE;

    for (size_t i=0; i < this->scoring_tools.size(); i++){
        std::vector<std::string> finalNumbering;
        std::string errorMessage = "";
        double percent_identity = -1;

        finalNumbering.reserve(query_sequence.length() * 2);

        this->scoring_tools[i]->align(query_sequence,
                    queryAsIdx.get(), finalNumbering, percent_identity,
                    errorMessage);
        if (percent_identity > best_identity){
            std::get<0>(best_result) = finalNumbering;
            std::get<1>(best_result) = percent_identity;
            std::get<2>(best_result) = this->scoring_tools[i]->get_chain_name();
            std::get<3>(best_result) = errorMessage;
            best_identity = percent_identity;
        }
        if (errorMessage.substr(0,2) == "> ")
            fwgxg_error = true;
    }
    if (fwgxg_error && best_identity < 0.8)
        return POSSIBLE_FWGXG_ERROR_ON_ALIGNMENT;
    return VALID_SEQUENCE;
}






// SingleChainAnnotatorCpp function which numbers a list of input sequences.
std::vector<std::tuple<std::vector<std::string>, double, std::string,
            std::string>> SingleChainAnnotatorCpp::analyze_seqs(std::vector<std::string> sequences){
    std::vector<std::tuple<std::vector<std::string>, double, std::string,
            std::string>> outputResults;

    for (size_t i=0; i < sequences.size(); i++){
        outputResults.push_back(this->analyze_seq(sequences[i]));
    }
    return outputResults;
}



// Testing function only; used to check how the needle scoring table is being filled.
void SingleChainAnnotatorCpp::_test_needle_scoring(std::string query_sequence,
                    nb::ndarray<double, nb::shape<-1,-1>, nb::device::cpu, nb::c_contig> scoreMatrix,
                    nb::ndarray<uint8_t, nb::shape<-1,-1>, nb::device::cpu, nb::c_contig> pathTraceMat,
                    std::string chain){
    for (size_t i=0; i < this->scoring_tools.size(); i++){
        if (this->scoring_tools[i]->get_chain_name() == chain){
            this->scoring_tools[i]->_test_fill_needle_scoring_table(query_sequence,
                    scoreMatrix, pathTraceMat);
            return;
        }
    }

}
