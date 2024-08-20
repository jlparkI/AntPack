#include "single_chain_annotator.h"


// Important scheme-specific defaults. Do not change unless
// necessary.
#define IMGT_DEFAULT_TERMINAL_TEMPLATE_GAP_PENALTY -1
#define IMGT_DEFAULT_C_TERMINAL_QUERY_GAP_PENALTY -1
#define MARTIN_DEFAULT_TERMINAL_TEMPLATE_GAP_PENALTY -1
#define MARTIN_DEFAULT_C_TERMINAL_QUERY_GAP_PENALTY -1
#define KABAT_DEFAULT_TERMINAL_TEMPLATE_GAP_PENALTY -1
#define KABAT_DEFAULT_C_TERMINAL_QUERY_GAP_PENALTY -1




SingleChainAnnotator::SingleChainAnnotator(
        std::vector<std::string> chains,
                std::string scheme, bool compress_init_gaps,
                bool multithread
):
    chains(chains),
    scheme(scheme),
    compress_init_gaps(compress_init_gaps),
    multithread(multithread)
{
    // Note that exceptions thrown here are go back to Python via
    // PyBind as long as this constructor is used within the wrapper.
    if (chains.size() < 1 || chains.size() > 3){
        throw std::runtime_error(std::string("There must be at least one chain specified, "
                    "and fewer than 3."));
    }

    std::filesystem::path extensionPath = std::filesystem::current_path();
    std::string uppercaseScheme = scheme;
    for (auto & c: uppercaseScheme) c = toupper(c);

    for (auto & chain : this->chains){
        if (chain != "H" && chain != "K" && chain != "L")
            throw std::runtime_error(std::string("Valid chains must be one of 'H', 'K', 'L'."));
        else{
            std::string npyFName = uppercaseScheme + "_CONSENSUS_" + chain + ".npy";
            std::string consFName = uppercaseScheme + "_CONSENSUS_" + chain + ".txt";
            std::filesystem::path npyFPath = extensionPath / npyFName;
            std::filesystem::path consFPath = extensionPath / consFName;

            std::vector<std::vector<std::string>> positionConsensus;
            int errCode = read_consensus_file(consFPath, positionConsensus);
            if (errCode != 1){
                throw std::runtime_error(std::string("The consensus file / library installation "
                            "has an issue."));
            }

            cnpy::NpyArray rawScoreArr = cnpy::npy_load(npyFPath.string());
            py::array_t<float, py::array::c_style> scoreArr({rawScoreArr.shape[0],
                    rawScoreArr.shape[1]});
            if (scheme == "imgt"){
                this->scoring_tools.push_back(std::make_unique<IGAligner>(scoreArr,
                            positionConsensus, chain, scheme,
                            IMGT_DEFAULT_TERMINAL_TEMPLATE_GAP_PENALTY,
                            IMGT_DEFAULT_C_TERMINAL_QUERY_GAP_PENALTY,
                            compress_init_gaps)
                        );
            }
            else if (scheme == "martin"){
                this->scoring_tools.push_back(std::make_unique<IGAligner>(scoreArr,
                            positionConsensus, chain, scheme,
                            MARTIN_DEFAULT_TERMINAL_TEMPLATE_GAP_PENALTY,
                            MARTIN_DEFAULT_C_TERMINAL_QUERY_GAP_PENALTY,
                            compress_init_gaps)
                        );
            }
            else if (scheme == "kabat"){
                this->scoring_tools.push_back(std::make_unique<IGAligner>(scoreArr,
                            positionConsensus, chain, scheme,
                            KABAT_DEFAULT_TERMINAL_TEMPLATE_GAP_PENALTY,
                            KABAT_DEFAULT_C_TERMINAL_QUERY_GAP_PENALTY,
                            compress_init_gaps)
                        );
            }
            else{
                throw std::runtime_error(std::string("Invalid scheme specified. Please use "
                            "one of 'imgt', 'kabat', 'martin'."));
            }
        }
    }
}




// SingleChainAnnotator function which numbers a single input sequence.
std::tuple<std::vector<std::string>, double, std::string,
            std::string> SingleChainAnnotator::analyze_seq(std::string sequence){
    
    std::vector<std::string> finalNumbering;
    double percentIdentity = 0;
    std::string chainName = "";
    std::string errMessage = "Invalid sequence supplied -- nonstandard AAs";
    std::tuple<std::vector<std::string>, double, std::string,
                        std::string> bestResult;

    if (!validate_sequence(sequence)){
        return std::tuple<std::vector<std::string>, double, std::string,
                        std::string>{finalNumbering,
                            percentIdentity, chainName, errMessage}; 
    }

    auto queryAsIdx = std::make_unique<int[]>( sequence.length() );

    if (!convert_sequence_to_array(queryAsIdx.get(), sequence)){
        return std::tuple<std::vector<std::string>, double, std::string,
                        std::string>{finalNumbering,
                            percentIdentity, chainName, errMessage}; 
    }

    if (this->scoring_tools.size() > 1 && this->multithread){
    }
    else{
        std::tuple<std::vector<std::string>, double, std::string,
                        std::string> currentResult;
        for (int i=0; i < this->scoring_tools.size(); i++){
            currentResult = this->scoring_tools[i]->align(sequence,
                    queryAsIdx.get());
            if (std::get<1>(currentResult) > percentIdentity)
                bestResult = currentResult;
        }
    }
    return bestResult;
}
