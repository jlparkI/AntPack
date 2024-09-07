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




SingleChainAnnotatorCpp::SingleChainAnnotatorCpp(
        std::vector<std::string> chains,
                std::string scheme, bool compress_init_gaps,
                bool multithread, std::string consensus_filepath
):
    AnnotatorBaseClassCpp(scheme),
    chains(chains),
    scheme(scheme),
    compress_init_gaps(compress_init_gaps),
    multithread(multithread)
{
    // Note that exceptions thrown here go back to Python via
    // PyBind as long as this constructor is used within the wrapper.
    if (chains.size() < 1 || chains.size() > 3){
        throw std::runtime_error(std::string("There must be at least one chain specified, "
                    "and fewer than 3."));
    }

    std::filesystem::path extensionPath = consensus_filepath;
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

            std::vector<std::vector<std::string>> position_consensus;
            int errCode = read_consensus_file(consFPath, position_consensus);
            if (errCode != 1){
                throw std::runtime_error(std::string("The consensus file / library installation "
                            "has an issue."));
            }

            cnpy::NpyArray raw_score_arr = cnpy::npy_load(npyFPath.string());
            double *raw_score_ptr = raw_score_arr.data<double>();
            if (raw_score_arr.word_size != 8){
                throw std::runtime_error(std::string("The consensus file / library installation "
                            "has an issue."));
            }

            py::array_t<double, py::array::c_style> scoreArr({raw_score_arr.shape[0],
                    raw_score_arr.shape[1]});
            double *scorePtr = static_cast<double*>(scoreArr.request().ptr);

            for (size_t k=0; k < raw_score_arr.shape[0] * raw_score_arr.shape[1]; k++)
                scorePtr[k] = raw_score_ptr[k];

            if (scheme == "imgt"){
                this->scoring_tools.push_back(std::make_unique<IGAligner>(scoreArr,
                            position_consensus, chain, scheme,
                            IMGT_DEFAULT_TERMINAL_TEMPLATE_GAP_PENALTY,
                            IMGT_DEFAULT_C_TERMINAL_QUERY_GAP_PENALTY,
                            compress_init_gaps)
                        );
            }
            else if (scheme == "aho"){
                this->scoring_tools.push_back(std::make_unique<IGAligner>(scoreArr,
                            position_consensus, chain, scheme,
                            AHO_DEFAULT_TERMINAL_TEMPLATE_GAP_PENALTY,
                            AHO_DEFAULT_C_TERMINAL_QUERY_GAP_PENALTY,
                            compress_init_gaps)
                        );
            }
            else if (scheme == "martin"){
                this->scoring_tools.push_back(std::make_unique<IGAligner>(scoreArr,
                            position_consensus, chain, scheme,
                            MARTIN_DEFAULT_TERMINAL_TEMPLATE_GAP_PENALTY,
                            MARTIN_DEFAULT_C_TERMINAL_QUERY_GAP_PENALTY,
                            compress_init_gaps)
                        );
            }
            else if (scheme == "kabat"){
                this->scoring_tools.push_back(std::make_unique<IGAligner>(scoreArr,
                            position_consensus, chain, scheme,
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
   
    double bestIdentity = -1;
    std::vector<std::string> emptyNumbering; 
    std::tuple<std::vector<std::string>, double, std::string,
                        std::string> bestResult{ emptyNumbering,
                            0, "", "Invalid sequence supplied -- nonstandard AAs"};


    if (!validate_sequence(sequence))
        return bestResult;

    auto queryAsIdx = std::make_unique<int[]>( sequence.length() );

    if (!convert_sequence_to_array(queryAsIdx.get(), sequence))
        return bestResult;

    if (this->scoring_tools.size() > 1 && this->multithread){
        std::vector<std::thread> threads;
        std::vector<std::vector<std::string>> threadNumberings(this->scoring_tools.size());
        std::vector<double> threadIdentities(this->scoring_tools.size(), -1);
        std::vector<std::string> threadErrMessages(this->scoring_tools.size(), "");

        auto queryCopies = std::make_unique<int[]>( sequence.length() *
               this->scoring_tools.size() );
        int *copyPtr = queryCopies.get();

        for (size_t i=0; i < this->scoring_tools.size(); i++){
            for (size_t j=0; j < sequence.length(); j++){
                int copyIdx = (i * sequence.length()) + j;
                queryCopies[copyIdx] = queryAsIdx[j];
            }

            threadNumberings[i].reserve(sequence.length() * 2);

            threads.emplace_back(&IGAligner::align, std::ref(*this->scoring_tools[i]),
                        sequence,
                        copyPtr, std::ref(threadNumberings.at(i)),
                        std::ref(threadIdentities.at(i)),
                        std::ref(threadErrMessages.at(i)));
            copyPtr += sequence.length();

        }
        for (auto& th : threads)
            th.join();

        for (size_t i=0; i < this->scoring_tools.size(); i++){
            if (threadIdentities[i] > bestIdentity){
                std::get<0>(bestResult) = threadNumberings[i];
                std::get<1>(bestResult) = threadIdentities[i];
                std::get<2>(bestResult) = this->scoring_tools[i]->get_chain_name();
                std::get<3>(bestResult) = threadErrMessages[i];
                bestIdentity = threadIdentities[i];
            }
        }
    }

    else{

        for (size_t i=0; i < this->scoring_tools.size(); i++){
            std::vector<std::string> finalNumbering;
            std::string errorMessage = "";
            double percentIdentity = -1;

            finalNumbering.reserve(sequence.length() * 2);

            this->scoring_tools[i]->align(sequence,
                    queryAsIdx.get(), finalNumbering, percentIdentity,
                    errorMessage);
            if (percentIdentity > bestIdentity){
                std::get<0>(bestResult) = finalNumbering;
                std::get<1>(bestResult) = percentIdentity;
                std::get<2>(bestResult) = this->scoring_tools[i]->get_chain_name();
                std::get<3>(bestResult) = errorMessage;
                bestIdentity = percentIdentity;
            }
        }
    }
    return bestResult;
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
