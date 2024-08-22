#include "single_chain_annotator.h"


// Important scheme-specific defaults. Do not change unless
// necessary.
#define IMGT_DEFAULT_TERMINAL_TEMPLATE_GAP_PENALTY -1
#define IMGT_DEFAULT_C_TERMINAL_QUERY_GAP_PENALTY -1
#define MARTIN_DEFAULT_TERMINAL_TEMPLATE_GAP_PENALTY -1
#define MARTIN_DEFAULT_C_TERMINAL_QUERY_GAP_PENALTY -1
#define KABAT_DEFAULT_TERMINAL_TEMPLATE_GAP_PENALTY -1
#define KABAT_DEFAULT_C_TERMINAL_QUERY_GAP_PENALTY -1




SingleChainAnnotatorCpp::SingleChainAnnotatorCpp(
        std::vector<std::string> chains,
                std::string scheme, bool compress_init_gaps,
                bool multithread, std::string project_filepath
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

    std::filesystem::path extensionPath = project_filepath;
    extensionPath /= "consensus_data";
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
            double *rawScorePtr = rawScoreArr.data<double>();
            if (rawScoreArr.word_size != 8){
                throw std::runtime_error(std::string("The consensus file / library installation "
                            "has an issue."));
            }

            py::array_t<double, py::array::c_style> scoreArr({rawScoreArr.shape[0],
                    rawScoreArr.shape[1]});
            double *scorePtr = static_cast<double*>(scoreArr.request().ptr);

            for (size_t k=0; k < rawScoreArr.shape[0] * rawScoreArr.shape[1]; k++)
                scorePtr[k] = rawScorePtr[k];

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



// SingleChainAnnotatorCpp function which sorts a list of position codes (essentially
// by calling the corresponding function in utilities). Essentially a wrapper on the
// corresponding utility function that can be accessed by Python or outside callers.
std::vector<std::string> SingleChainAnnotatorCpp::sort_position_codes(std::vector<std::string> position_code_list){
    std::vector<std::string> cleanedCodes, sortedCodes;
    int isValid;

    for (auto &code : position_code_list){
        if (code != "-")
            cleanedCodes.push_back(code);
    }
    isValid = sort_position_codes_utility(cleanedCodes, this->scheme,
            sortedCodes);
    if (isValid == 0){
        throw std::runtime_error(std::string("One or more of the supplied position "
                    "codes is invalid given the specified scheme."));
    }
    return sortedCodes;
}



// SingleChainAnnotatorCpp function which builds a multiple sequence alignment from
// a list of previously numbered sequences. The sequences must all be of the
// same chain type. Essentially a wrapper on the corresponding utility function
// that can be accessed by Python or outside callers.
std::tuple<std::vector<std::string>, std::vector<std::string>> SingleChainAnnotatorCpp::build_msa(std::vector<std::string> sequences,
        std::vector<std::tuple<std::vector<std::string>, double, std::string, std::string>> annotations){
    std::vector<std::string> positionCodes, alignedSeqs;
    int errCode;

    errCode = build_msa_utility(sequences, annotations, positionCodes,
            alignedSeqs, this->scheme);
    if (errCode != 1)
        throw std::runtime_error(std::string("A fatal error occured when building an MSA -- please report."));
    return std::tuple<std::vector<std::string>, std::vector<std::string>>{positionCodes, alignedSeqs};
}


// SingleChainAnnotatorCpp function which trims an alignment to remove gap positions at
// either end. Essentially a wrapper on the corresponding utility function that can
// be accessed by Python or outside callers.
std::tuple<std::string, std::vector<std::string>, int, int> SingleChainAnnotatorCpp::trim_alignment(std::string sequence,
        std::tuple<std::vector<std::string>, double, std::string, std::string> alignment){
    int exstart, exend;
    std::vector<char> trimmedSeq;
    std::vector<std::string> trimmedAlignment;

    int errCode = trim_alignment_utility(sequence, alignment, trimmedAlignment,
            exstart, exend, trimmedSeq);
    if (errCode != 1)
        throw std::runtime_error(std::string("Invalid sequence / alignment pairing supplied."));

    std::string trimmedSeqStr(trimmedSeq.begin(), trimmedSeq.end());

    return std::tuple<std::string, std::vector<std::string>, int, int>{trimmedSeqStr,
        trimmedAlignment, exstart, exend};
}


// SingleChainAnnotatorCpp function that assigns cdr and framework labels based on an
// alignment performed by analyze_seq (or similar). It can use a scheme other than
// the one used by the SingleChainAnnotatorCpp object, so that the user can assign
// cdr labels for an IMGT-numbered scheme using Kabat CDR definitions if desired.
std::vector<std::string> SingleChainAnnotatorCpp::assign_cdr_labels(std::tuple<std::vector<std::string>, 
               double, std::string, std::string> alignment, std::string cdr_scheme){
    std::string currentScheme;
    std::vector<std::string> cdrLabeling;

    if (cdr_scheme == "")
        currentScheme = this->scheme;
    else if (cdr_scheme == "imgt" || cdr_scheme == "kabat" || cdr_scheme == "martin" ||
            cdr_scheme == "aho")
        currentScheme = cdr_scheme;
    else{
        throw std::runtime_error(std::string("Unrecognized scheme supplied. Use '' to use "
                    "the scheme currently selected for this Annotator, or use a currently "
                    "accepted scheme."));
    }

    return cdrLabeling;
}
