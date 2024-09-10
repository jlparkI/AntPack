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

    if (!this->align_input_subregion(best_result, best_identity,
                sequence))
        return best_result;


    // It can (rarely) happen that we have a chain (usually light) where the expected FGxG
    // motif in the J-gene is altered AND there is significant additional sequence beyond
    // the end of the J-gene, e.g. in a paired chain AND this results in an alignment
    // error. In the event this unusual combination of circumstances occurs, it
    // usually manifests in the form of an unusually long CDR that exceeds the maximum
    // number of allowed insertion codes. In such cases we 1) find regions that likely
    // represent the c-terminus of a variable region, 2) use the first such to split up
    // the input into regions, 3) align each region and 4) keep the one with the best pid. This
    // is significantly more expensive than a simple alignment but since this issue
    // is so rare it has an insignificant effect on average numbering time.
    if (std::get<3>(best_result).substr(0, 2) == "> "){
        printf("\nBOING BOING!!\n");
        // Using 3 here and throughout since cterm finder looks for a cterminal corresponding
        // to any of the three possible chains, even if our aligner is looking for one only.
        std::array<double, 3> mscores;
        std::array<int, 3> mpositions;
        int cterm_err_code = this->boundary_finder->find_c_terminals(sequence, mscores,
                mpositions);
        // Usually cterm finder error results from invalid characters in the input
        // sequence. If this happens for some horribly unlikely alternative reason,
        // abort and return the existing alignment with the existing error code.
        if (cterm_err_code != 1)
            return best_result;
        else if (mscores[0] == 0 && mscores[1] == 0 && mscores[2] == 0)
            return best_result;

        // Find the FIRST likely c-terminal region in the input sequence.
        int first_cterm_position = 100000;

        for (int i=0; i < 3; i++){
            if (first_cterm_position > mpositions[i])
                first_cterm_position = mpositions[i];
        }
        // Add an offset from the beginning of the likely c-terminal region.
        first_cterm_position += CTERMINAL_OFFSET_ALIGNMENT_SHIFT;
        if (first_cterm_position > (int)sequence.length())
            first_cterm_position = sequence.length();

        std::string sequence_extract;
        double second_identity = -1;
        std::vector<std::string> second_numbering; 
        std::tuple<std::vector<std::string>, double, std::string,
                        std::string> second_result{ second_numbering,
                            0, "", ""};

        if (first_cterm_position > 85 && (sequence.length() - first_cterm_position) > 85){
            sequence_extract = sequence.substr(0, first_cterm_position);
            if (!this->align_input_subregion(best_result, best_identity,
                    sequence_extract))
                return best_result;

            // Shift the cterm position back by a couple in case it is off by one or two.
            first_cterm_position -= 2;
            sequence_extract = sequence.substr(first_cterm_position,
                    sequence.length() - first_cterm_position);
            if (!this->align_input_subregion(second_result, second_identity,
                    sequence_extract))
                return best_result;

            if (second_identity > best_identity){
                this->pad_left(second_result, sequence);
                return second_result;
            }
            else
                this->pad_right(best_result, sequence);
        }

        else if (first_cterm_position > 85){
            sequence_extract = sequence.substr(0, first_cterm_position);
            int err_code = this->align_input_subregion(best_result, best_identity,
                    sequence_extract);

            if (err_code != 1)
                return best_result;
            this->pad_right(best_result, sequence);
        }

        else if ((sequence.length() - first_cterm_position) > 85){
            // Shift the cterm position back by a couple in case it is off by one or two.
            first_cterm_position -= 2;
            if (first_cterm_position < 0)
                first_cterm_position = 0;
            sequence_extract = sequence.substr(first_cterm_position,
                    sequence.length() - first_cterm_position);
            int err_code = this->align_input_subregion(best_result, best_identity,
                    sequence_extract);

            if (err_code != 1)
                return best_result;
            this->pad_left(best_result, sequence);
        }
    }
    
    return best_result;
}


// Aligns the input sequence, which may be either the full sequence or a subregion of it.
int SingleChainAnnotatorCpp::align_input_subregion(std::tuple<std::vector<std::string>, double,
                std::string, std::string> &best_result, double &best_identity,
                std::string &query_sequence){
    auto queryAsIdx = std::make_unique<int[]>( query_sequence.length() );

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
        if (errorMessage.substr(0, 2) == "> ")
            printf("\nOREOS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
        if (percent_identity > best_identity){
            std::get<0>(best_result) = finalNumbering;
            std::get<1>(best_result) = percent_identity;
            std::get<2>(best_result) = this->scoring_tools[i]->get_chain_name();
            std::get<3>(best_result) = errorMessage;
            best_identity = percent_identity;
        }
    }
    return VALID_SEQUENCE;
}



// Pads input numbering on the left to make it the same length as the
// input sequence. Only used when realigning to fix rare alignment errors.
void SingleChainAnnotatorCpp::pad_left(std::tuple<std::vector<std::string>, double,
            std::string, std::string> &alignment,
        std::string &query_sequence){

    std::vector<std::string> &numbering = std::get<0>(alignment);

    int num_gaps = (query_sequence.length() - numbering.size());
    if (num_gaps <= 0)
        return;

    std::vector<std::string> updated_numbering(num_gaps, "-");

    for (size_t i=0; i < numbering.size(); i++)
        updated_numbering.push_back(numbering[i]);

    std::get<0>(alignment) = updated_numbering;
}


// Pads input numbering on the right to make it the same length as the
// input sequence. Only used when realigning to fix rare alignment errors.
void SingleChainAnnotatorCpp::pad_right(std::tuple<std::vector<std::string>, double,
            std::string, std::string> &alignment,
        std::string &query_sequence){

    std::vector<std::string> &numbering = std::get<0>(alignment);
    int num_gaps = (query_sequence.length() - numbering.size());
    if (num_gaps <= 0)
        return;

    for (int i=0; i < num_gaps; i++)
        numbering.push_back("-");
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
