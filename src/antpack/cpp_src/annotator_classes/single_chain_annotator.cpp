#include "single_chain_annotator.h"



// Special error code for cases where a rare but specific type
// of alignment error may have occurred.
#define POSSIBLE_FWGXG_ERROR_ON_ALIGNMENT 2



namespace NumberingTools{


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
                    this->scoring_tools.push_back(IGAligner(consensus_filepath,
                            chain, scheme,
                            IMGT_DEFAULT_TERMINAL_TEMPLATE_GAP_PENALTY,
                            IMGT_DEFAULT_C_TERMINAL_QUERY_GAP_PENALTY,
                            compress_init_gaps)
                        );
                }
                else if (scheme == "aho"){
                    this->scoring_tools.push_back(IGAligner(consensus_filepath,
                            chain, scheme,
                            AHO_DEFAULT_TERMINAL_TEMPLATE_GAP_PENALTY,
                            AHO_DEFAULT_C_TERMINAL_QUERY_GAP_PENALTY,
                            compress_init_gaps)
                        );
                }
                else if (scheme == "martin"){
                    this->scoring_tools.push_back(IGAligner(consensus_filepath,
                            chain, scheme,
                            MARTIN_DEFAULT_TERMINAL_TEMPLATE_GAP_PENALTY,
                            MARTIN_DEFAULT_C_TERMINAL_QUERY_GAP_PENALTY,
                            compress_init_gaps)
                        );
                }
                else if (scheme == "kabat"){
                    this->scoring_tools.push_back(IGAligner(consensus_filepath,
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



    /// SingleChainAnnotatorCpp function which numbers a single input sequence.
    std::tuple<std::vector<std::string>, double, std::string,
            std::string> SingleChainAnnotatorCpp::analyze_seq(std::string sequence){
   
        double best_identity = -1;
        std::vector<std::string> emptyNumbering; 
        std::tuple<std::vector<std::string>, double, std::string,
                        std::string> best_result{ emptyNumbering,
                            0, "", "Invalid sequence supplied -- nonstandard AAs"};


        if (!SequenceUtilities::validate_x_sequence(sequence))
            return best_result;

        int err_code = this->align_input_subregion(best_result, best_identity,
                sequence);
        if (err_code !=  POSSIBLE_FWGXG_ERROR_ON_ALIGNMENT)
            return best_result;

        // It can (rarely) happen that we have a chain (usually light) where the expected FGxG
        // motif in the J-gene is altered AND there is significant additional sequence beyond
        // the end of the J-gene, e.g. in a paired chain, AND this results in an alignment
        // error. In the event this unusual combination of circumstances occurs, it
        // usually manifests in the form of an unusually long CDR that exceeds the maximum
        // number of allowed insertion codes.

        // If this may have occurred, subdivide the sequence into regions,
        // align each region and return the one with the highest
        // percent identity. This is not ideal if the sequence contains
        // multiple chains because the other chains will go unreported,
        // but a SingleChainAnnotator by definition is not ideal for
        // this scenario.
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
            if (percent_identity < best_identity || err_code != SequenceUtilities::VALID_SEQUENCE)
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




    /// Aligns the input sequence, which may be either the full sequence or a subregion of it.
    int SingleChainAnnotatorCpp::align_input_subregion(std::tuple<std::vector<std::string>, double,
                    std::string, std::string> &best_result, double &best_identity,
                    std::string &query_sequence){
        auto queryAsIdx = std::make_unique<int[]>( query_sequence.length() );
        bool fwgxg_error = false;

        if (!SequenceUtilities::convert_x_sequence_to_array(queryAsIdx.get(), query_sequence))
            return INVALID_SEQUENCE;

        for (size_t i=0; i < this->scoring_tools.size(); i++){
            std::vector<std::string> final_numbering;
            std::string error_message = "";
            double percent_identity = -1;

            final_numbering.reserve(query_sequence.length() * 2);

            this->scoring_tools[i].align(query_sequence,
                        queryAsIdx.get(), final_numbering, percent_identity,
                        error_message);
            if (percent_identity > best_identity){
                std::get<0>(best_result) = final_numbering;
                std::get<1>(best_result) = percent_identity;
                std::get<2>(best_result) = this->scoring_tools[i].get_chain_name();
                std::get<3>(best_result) = error_message;
                best_identity = percent_identity;
            }
            if (error_message.length() > 2){
                if (error_message.substr(0,2) == "> ")
                    fwgxg_error = true;
            }
        }
        if (fwgxg_error && best_identity < 0.8)
            return POSSIBLE_FWGXG_ERROR_ON_ALIGNMENT;
        return VALID_SEQUENCE;
    }






    // SingleChainAnnotatorCpp function which numbers a list of input sequences.
    std::vector<std::tuple<std::vector<std::string>, double, std::string,
                std::string>> SingleChainAnnotatorCpp::analyze_seqs(std::vector<std::string> sequences){
        std::vector<std::tuple<std::vector<std::string>, double, std::string,
                std::string>> output_results;

        for (size_t i=0; i < sequences.size(); i++)
            output_results.push_back(this->analyze_seq(sequences[i]));
        
        return output_results;
    }

}
