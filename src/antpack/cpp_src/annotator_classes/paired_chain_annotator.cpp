/* Contains the implementation of the paired chain annotator tool for
 * sequences containing a paired heavy and light chain. */
#include "paired_chain_annotator.h"





PairedChainAnnotatorCpp::PairedChainAnnotatorCpp(
        std::string scheme, std::string consensus_filepath
):
    AnnotatorBaseClassCpp(scheme),
    scheme(scheme)
{
    std::vector<std::string> chains = {"K", "L"};
    this->light_chain_analyzer = std::make_unique<SingleChainAnnotatorCpp>(chains,
            scheme, false, consensus_filepath);

    chains = {"H"};
    this->heavy_chain_analyzer = std::make_unique<SingleChainAnnotatorCpp>(chains,
            scheme, false, consensus_filepath);

    chains = {"H", "K", "L"};
    this->analyzer = std::make_unique<SingleChainAnnotatorCpp>(chains,
            scheme, false, consensus_filepath);


    std::filesystem::path extensionPath = consensus_filepath;
    
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

    chains = {"H", "K", "L"};
    int shape0 = raw_score_arr.shape[0], shape1 = raw_score_arr.shape[1];
    if (shape1 != 20){
        throw std::runtime_error(std::string("The consensus file / library installation "
                            "has an issue."));
    }
    py::array_t<double, py::array::c_style> score_array({shape0, shape1, 3});
    auto scoreMatItr = score_array.mutable_unchecked<3>();

    for (size_t i=0; i < chains.size(); i++){
        std::string npyFName = "CTERMFINDER_CONSENSUS_" + chains[i] + ".npy";
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



// PairedChainAnnotatorCpp function which numbers a list of input sequences.
// An alignment is returned for each input.
std::pair<std::tuple<std::vector<std::string>, double, std::string, std::string>,
    std::tuple<std::vector<std::string>, double, std::string, std::string>>
    PairedChainAnnotatorCpp::analyze_seq(std::string sequence){
    
    // Using 3 here and throughout since there are three possible chains.
    std::array<double, 3> mscores;
    std::array<int, 3> mpositions;

    // Find likely c-terminal regions in the input sequence.
    int err_code = this->boundary_finder->find_c_terminals(sequence, mscores,
            mpositions);

    if (err_code != 1){
        std::tuple<std::vector<std::string>, double, std::string, std::string> vecres =
            {{}, 0., "", "Invalid sequence supplied -- nonstandard AAs"};
        std::pair<std::tuple<std::vector<std::string>, double, std::string, std::string>,
            std::tuple<std::vector<std::string>, double, std::string, std::string>> output_result =
            {vecres, vecres};
        return output_result;
    }
    else if (mscores[0] == 0 && mscores[1] == 0 && mscores[2] == 0){
        std::tuple<std::vector<std::string>, double, std::string, std::string> vecres =
            {{}, 0., "", "Invalid sequence supplied -- nonstandard AAs"};
        std::pair<std::tuple<std::vector<std::string>, double, std::string, std::string>,
            std::tuple<std::vector<std::string>, double, std::string, std::string>> output_result =
            {vecres, vecres};
        return output_result;
    }

    // Find the FIRST likely c-terminal region in the input sequence.
    int start_position = 100000;

    for (int i=0; i < 3; i++){
        if (start_position > mpositions[i])
            start_position = mpositions[i];
    }

    // Add an offset from the beginning of the likely c-terminal region.
    start_position += CTERMINAL_STANDARD_OFFSET;

    // If the start position is less than 85, this is a severely truncated sequence
    // (variable chains are normally > 100 in length). Likewise, if the start position
    // is close to the end of the sequence. In either case, do not attempt to number,
    // because results are not likely to be useful.
    if (start_position < 85 || (sequence.length() - start_position) < 85){
        std::tuple<std::vector<std::string>, double, std::string, std::string> vecres =
            {{}, 0., "", "Could not find c-terminal of first variable chain -- "
                "sequence may not be a paired chain sequence"};
        std::pair<std::tuple<std::vector<std::string>, double, std::string, std::string>,
            std::tuple<std::vector<std::string>, double, std::string, std::string>> output_result =
            {vecres, vecres};
        return output_result;
    }

    std::string sequence_extract = sequence.substr(start_position,
            sequence.length() - start_position);

    std::tuple<std::vector<std::string>, double, std::string, std::string> init_results =
        this->analyzer->analyze_seq(sequence_extract);
    
    std::vector<std::string> padded_init_numbering((sequence.length() -
                std::get<0>(init_results).size()), "-");

    for (size_t i=0; i < std::get<0>(init_results).size(); i++)
        padded_init_numbering.push_back(std::get<0>(init_results).at(i));

    std::get<0>(init_results) = std::move(padded_init_numbering);

    size_t exstart = 0;

    for (size_t i=0; i < std::get<0>(init_results).size(); i++){
        if (std::get<0>(init_results)[i] != "-"){
            exstart = i;
            break;
        }
    }

    // Unlikely since we already checked that the start position for the first
    // alignment was at least 85 from start of sequence, so exstart should always
    // be > 85. Nonetheless, check just in case.
    if (exstart < 85){
        std::tuple<std::vector<std::string>, double, std::string, std::string> vecres =
            {{}, 0., "", "Could not find c-terminal of first variable chain -- "
                "sequence may not be a paired chain sequence"};
        std::pair<std::tuple<std::vector<std::string>, double, std::string, std::string>,
            std::tuple<std::vector<std::string>, double, std::string, std::string>> output_result =
            {vecres, vecres};
        return output_result;
    }

    sequence_extract = sequence.substr(0, exstart);
    std::tuple<std::vector<std::string>, double, std::string, std::string> second_result;
    

    if (std::get<2>(init_results) == "H")
        second_result = this->light_chain_analyzer->analyze_seq(sequence_extract);
    else
        second_result = this->heavy_chain_analyzer->analyze_seq(sequence_extract);

    size_t padding_amount = (sequence.length() - std::get<0>(second_result).size());

    for (size_t i=0; i < padding_amount; i++)
        std::get<0>(second_result).push_back("-");

    std::pair<std::tuple<std::vector<std::string>, double, std::string, std::string>,
            std::tuple<std::vector<std::string>, double, std::string, std::string>> output_result;

    if (std::get<2>(init_results) == "H")
        output_result = {init_results, second_result};

    else
        output_result = {second_result, init_results};

    return output_result;
}
