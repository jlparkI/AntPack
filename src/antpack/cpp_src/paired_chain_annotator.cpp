/* Contains the implementation of the paired chain annotator tool for
 * sequences containing a paired heavy and light chain. */
#include "paired_chain_annotator.h"





PairedChainAnnotatorCpp::PairedChainAnnotatorCpp(
        std::string scheme, bool multithread,
        std::string consensus_filepath
):
    AnnotatorBaseClassCpp(scheme),
    scheme(scheme),
    multithread(multithread),
    project_filepath(project_filepath)
{
    std::vector<std::string> chains = {"K", "L"};
    this->light_chain_analyzer = std::make_unique<SingleChainAnnotatorCpp>(chains,
            scheme, false, multithread, project_filepath);

    chains = {"H"};
    this->heavy_chain_analyzer = std::make_unique<SingleChainAnnotatorCpp>(chains,
            scheme, false, multithread, project_filepath);

    chains = {"H", "K", "L"};
    this->analyzer = std::make_unique<SingleChainAnnotatorCpp>(chains,
            scheme, false, multithread, project_filepath);


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
    double *scorePtr = static_cast<double*>(score_array.request().ptr);

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
        if (raw_score_arr.word_size != 8 || raw_score_arr.shape[0] != shape0 ||
                raw_score_arr.shape[1] != shape1){
            throw std::runtime_error(std::string("The consensus file / library installation "
                            "has an issue."));
        }

        double *raw_score_ptr = raw_score_arr.data<double>();
        int copy_start = i * shape0 * shape1;

        for (size_t k=0; k < (shape0 * shape1); k++)
            scorePtr[k + copy_start] = raw_score_ptr[k];
    }

    this->boundary_finder = std::make_unique<CTermFinder>(score_array);
}



// PairedChainAnnotatorCpp function which numbers a list of input sequences.
// An alignment is returned for each input.
std::vector<std::tuple<std::vector<std::string>, double, std::string,
            std::string>> PairedChainAnnotatorCpp::analyze_seqs(std::vector<std::string> sequences){
    std::vector<std::tuple<std::vector<std::string>, double, std::string,
            std::string>> outputResults;

    for (size_t i=0; i < sequences.size(); i++){
        outputResults.push_back(this->analyze_seq(sequences[i]));
    }
    return outputResults;
}
