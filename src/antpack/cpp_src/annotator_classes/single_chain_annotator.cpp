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


// Offset on identified c-terminals if performing a realignment.
#define CTERMINAL_OFFSET_ALIGNMENT_SHIFT 11


// Special error code for cases where a rare but specific type
// of alignment error may have occurred.
#define POSSIBLE_FWGXG_ERROR_ON_ALIGNMENT 2



SingleChainAnnotatorCpp::SingleChainAnnotatorCpp(
        std::vector<std::string> chains,
                std::string scheme, bool compress_init_gaps,
                std::string consensus_filepath
):
    AnnotatorBaseClassCpp(scheme),
    chains(chains),
    scheme(scheme),
    compress_init_gaps(compress_init_gaps)
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
    int shape0 = raw_score_arr.shape[0], shape1 = raw_score_arr.shape[1];
    if (shape1 != 20){
        throw std::runtime_error(std::string("The consensus file / library installation "
                            "has an issue."));
    }
    py::array_t<double, py::array::c_style> score_array({shape0, shape1, 3});
    auto scoreMatItr = score_array.mutable_unchecked<3>();

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
        if (raw_score_arr.word_size != 8 || ld_shape0 != shape0 || ld_shape1 != shape1)
            throw std::runtime_error(std::string("The consensus file / library installation "
                            "has an issue."));

        double *raw_score_ptr = raw_score_arr.data<double>();

        for (int j=0; j < shape0; j++){
            for (int k=0; k < shape1; k++){
                scoreMatItr(j,k,i) = *raw_score_ptr;
                raw_score_ptr++;
            }
        }
    }

    this->boundary_finder = std::make_unique<CTermFinder>(score_array);

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

    // If this may have occurred, check for the location of likely j-genes,
    // use them to break up the long sequence and align the part that most likely
    // contains the viable variable region.

    // Using 3 here and throughout since cterm finder looks for a cterminal corresponding
    // to any of the three possible chains, even if our aligner is looking for one or two
    // only. Consequently, when analyzing the results, we have to be careful to check that
    // we only use results corresponding to chains selected for this analyzer.
    std::array<double, 3> mscores;
    std::array<int, 3> mpositions;
    err_code = this->boundary_finder->find_c_terminals(sequence, mscores,
                mpositions);
    // Usually cterm finder error results from invalid characters in the input
    // sequence or a sequence that does not have at least two cysteines.
    // If some very unusual problem has occurred, simply return the existing
    // problematic alignment rather than trying to fix an alignment for
    // a sequence that has some serious issues.
    if (err_code != 1)
        return best_result;
    else if (mscores[0] == 0 && mscores[1] == 0 && mscores[2] == 0)
        return best_result;

    int min_starting_position = 0;
    for (size_t i=0; i < std::get<0>(best_result).size(); i++){
        if (std::get<0>(best_result)[i] != "-"){
            min_starting_position = i;
            break;
        }
    }

    // The ordering H, K, L for cterm finder results is set in the
    // SingleChainAnnotator class constructor.
    int best_cterm_position = 100000;
    for (size_t i=0; i < this->chains.size(); i++){
        if (this->chains[i] == "H"){
            if (mpositions[0] < best_cterm_position &&
                    (mpositions[0] > min_starting_position + 40))
                best_cterm_position = mpositions[0];
        }
        else if (this->chains[i] == "K"){
            if (mpositions[1] < best_cterm_position &&
                    (mpositions[1] > min_starting_position + 40))
                best_cterm_position = mpositions[1];
        }
        else if (this->chains[i] == "L"){
            if (mpositions[2] < best_cterm_position &&
                    (mpositions[2] > min_starting_position + 40))
                best_cterm_position = mpositions[2];
        }
    }
    if (best_cterm_position == 100000)
        return best_result;

    // Add an offset from the beginning of the likely c-terminal region.
    best_cterm_position += CTERMINAL_OFFSET_ALIGNMENT_SHIFT;

    std::string sequence_extract = sequence.substr(0, best_cterm_position);
    err_code = this->align_input_subregion(best_result, best_identity,
                sequence_extract);
    this->pad_right(best_result, sequence);

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
    if (fwgxg_error && best_identity < 0.75)
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
