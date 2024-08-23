/* Contains the implementation of the paired chain annotator tool for
 * sequences containing a paired heavy and light chain. */
#include "paired_chain_annotator.h"



PairedChainAnnotatorCpp::PairedChainAnnotatorCpp(
        std::string scheme, bool multithread,
        std::string project_filepath
):
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
}
