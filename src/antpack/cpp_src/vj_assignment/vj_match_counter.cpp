#include "vj_match_counter.h"


VJMatchCounter::VJMatchCounter(
                std::map<std::string, std::vector<std::string>> gene_names,
                std::map<std::string, std::vector<std::string>> gene_seqs,
                py::array_t<int16_t, py::array::c_style> blosum_matrix,
                std::string scheme
):
    gene_seqs(gene_seqs),
    gene_names(gene_names),
    blosum_matrix(blosum_matrix),
    scheme(scheme)
{
    
    // Note that exceptions thrown here go back to Python via
    // PyBind as long as this constructor is used within the wrapper.
    py::buffer_info blosum_info = blosum_matrix.request();
    if (blosum_info.shape.size() != 2){
        throw std::runtime_error(std::string("Error with package installation; "
                    "incorrect BLOSUM matrix supplied."));
    }
    if (blosum_info.shape[0] != 22 || blosum_info.shape[1] != 22){
        throw std::runtime_error(std::string("Error with package installation; "
                    "incorrect BLOSUM matrix supplied."));
    }

    for ( const auto &gene_seq_element : gene_seqs ) {
        std::map<std::string, std::vector<std::string>>::iterator gene_name_element =
            gene_names.find(gene_seq_element.first);

        if (gene_name_element == gene_names.end()){
            throw std::runtime_error(std::string("One or more keys in the gene_seqs "
                        "dictionary is not in the gene_names dictionary."));
        }
        if (gene_name_element->second.size() != gene_seq_element.second.size()){
            throw std::runtime_error(std::string("The gene name lists and gene sequence "
                        "lists must be of the same sizes."));
        }
        for (const auto &gene_sequence : gene_seq_element.second){
            if (gene_sequence.length() != REQUIRED_SEQUENCE_LENGTH){
                throw std::runtime_error(std::string("All sequences passed to VJMatchCounter "
                        "must have the correct length."));
            }
        }
    }

    for ( const auto &gene_name_element : gene_names ) {
        if (gene_seqs.find(gene_name_element.first) == gene_seqs.end()){
            throw std::runtime_error(std::string("One or more keys in the "
                        "gene_names dictionary is not in the gene_seqs dictionary."));
        }
        std::map<std::string, int> name_submap;
        for (size_t i=0; i < gene_name_element.second.size(); i++)
            name_submap[gene_name_element.second[i]] = i;

        this->names_to_positions[gene_name_element.first] = std::move(name_submap);
    }

    // Initialize the set of expected imgt positions.
    for (size_t i=1; i < REQUIRED_SEQUENCE_LENGTH + 1; i++)
        this->essential_imgt_map.insert(std::to_string(i));

    // If we are using a scheme OTHER than IMGT, we will map the key positions needed for
    // VJ matching in this other scheme to IMGT using a function from utilities. For now,
    // check to make sure scheme is accepted.
    if (this->scheme != "imgt" && this->scheme != "aho"
            && this->scheme != "kabat" && this->scheme != "martin"){
        throw std::runtime_error(std::string("Unrecognized scheme supplied. "
                    "Scheme must be one of imgt, martin, kabat or aho."));
    }
}



std::tuple<std::string, std::string,
    double, double> VJMatchCounter::assign_vj_genes(std::tuple<std::vector<std::string>,
        double, std::string, std::string> alignment, std::string sequence,
        std::string species, std::string mode){

    if (species != "human" && species != "mouse"){
        throw std::runtime_error(std::string("Species for VJ gene assignment must be one of "
                    "'human', 'mouse'."));
    }

    if (std::get<2>(alignment) != "H" && std::get<2>(alignment) != "K" &&
            std::get<2>(alignment) != "L"){
        throw std::runtime_error(std::string("Chain must be one of 'H', 'K', 'L'."));
    }

    if (sequence.length() != std::get<0>(alignment).size()){
        throw std::runtime_error(std::string("Numbering and sequence must have the "
                    "same length."));
    }

    // Notice that for now we assume immunoglobulin (not TCR).
    std::string vkey = species + "_IG" + std::get<2>(alignment) + "V";
    std::string jkey = species + "_IG" + std::get<2>(alignment) + "J";
    
    std::map<std::string, std::vector<std::string>>::iterator vname_dict_finder =
            this->gene_names.find(vkey);
    std::map<std::string, std::vector<std::string>>::iterator jname_dict_finder =
            this->gene_names.find(jkey);
    std::map<std::string, std::vector<std::string>>::iterator vgene_dict_finder =
            this->gene_seqs.find(vkey);
    std::map<std::string, std::vector<std::string>>::iterator jgene_dict_finder =
            this->gene_seqs.find(jkey);
        
    if (vname_dict_finder == gene_names.end() || jname_dict_finder == gene_names.end()
            || vgene_dict_finder == gene_seqs.end() || jgene_dict_finder == gene_seqs.end()){
        throw std::runtime_error(std::string("There was an error in the construction of the "
                    "vj gene tool; necessary gene databases were not found. Please report."));
    }

    std::vector<std::string> &vgenes = vgene_dict_finder->second;
    std::vector<std::string> &jgenes = jgene_dict_finder->second;
    std::vector<std::string> &vnames = vname_dict_finder->second;
    std::vector<std::string> &jnames = jname_dict_finder->second;

    // This string will store the positions in the sequence that correspond to IMGT 1-128.
    // If the scheme is not IMGT, we can convert the numbering to IMGT temporarily for
    // this purpose.
    std::string prepped_sequence(REQUIRED_SEQUENCE_LENGTH, '-');
    this->prep_sequence(prepped_sequence, sequence,
            std::get<0>(alignment));

    if (vgenes.size() != vnames.size() || jgenes.size() != jnames.size()){
        throw std::runtime_error(std::string("There was an error in the construction of the "
                    "vj gene tool; necessary gene databases were not found. Please report."));
    }

    std::string vgene_name, jgene_name;
    double videntity, jidentity;

    if (mode == "identity"){
        this->assign_gene_by_identity(vgenes, vnames, prepped_sequence,
                videntity, vgene_name, 'v');
        this->assign_gene_by_identity(jgenes, jnames, prepped_sequence,
                jidentity, jgene_name, 'j');
    }
    else if (mode == "evalue"){
        auto encoded_query = std::make_unique<int[]>( REQUIRED_SEQUENCE_LENGTH );
        int err_code = convert_sequence_to_array(encoded_query.get(), prepped_sequence);

        if (err_code != 1){
            throw std::runtime_error(std::string("The input sequence contains invalid "
                        "amino acids and could not be encoded."));
        }

        this->assign_gene_by_evalue(vgenes, vnames, encoded_query.get(),
                videntity, vgene_name, 'v');
        this->assign_gene_by_evalue(jgenes, jnames, encoded_query.get(),
                jidentity, jgene_name, 'j');
    }
    else{
        throw std::runtime_error(std::string("Unrecognized mode was supplied. "
                    "Mode should be one of 'identity', 'evalue'."));
    }

    return std::tuple<std::string, std::string,
           double, double>{vgene_name, jgene_name,
                videntity, jidentity};
}



void VJMatchCounter::assign_gene_by_identity(std::vector<std::string> &gene_seqs,
                std::vector<std::string> &gene_names,
                std::string &prepped_sequence,
                double &best_identity,
                std::string &best_gene_name,
                char gene_type){

    best_identity = 0;
    int start_letter = 0, end_letter = 0;

    if (gene_type == 'j'){
        start_letter = 105;
        end_letter = REQUIRED_SEQUENCE_LENGTH;
    }
    else{
        start_letter = 0;
        end_letter = 108;
    }

    int matching_positions, nonzero_positions, closest_id = 0;
    double current_identity = 0;

    for (size_t i=0; i < gene_seqs.size(); i++){
        matching_positions = 0;
        nonzero_positions = 0;

        for (int j=start_letter; j < end_letter; j++){
            if (gene_seqs[i][j] == '-')
                continue;
            nonzero_positions += 1;
            if (gene_seqs[i][j] == prepped_sequence[j])
                matching_positions += 1;
        }
        if (nonzero_positions == 0)
            nonzero_positions = 1;
        current_identity = (static_cast<double>(matching_positions)) / 
            (static_cast<double>(nonzero_positions));

        if (current_identity > best_identity){
            closest_id = i;
            best_identity = current_identity;
        }
    }

    best_gene_name = gene_names[closest_id];
}




int VJMatchCounter::assign_gene_by_evalue(std::vector<std::string> &gene_seqs,
                std::vector<std::string> &gene_names,
                int *encoded_sequence,
                double &best_identity,
                std::string &best_gene_name,
                char gene_type){

    best_identity = 21 * REQUIRED_SEQUENCE_LENGTH;
    int start_letter = 0, end_letter = 0;

    if (gene_type == 'j'){
        start_letter = 105;
        end_letter = REQUIRED_SEQUENCE_LENGTH;
    }
    else{
        start_letter = 0;
        end_letter = 108;
    }
    
    auto blosum_itr = this->blosum_matrix.unchecked<2>();
    int blosum_distance, closest_id = 0;

    for (size_t i=0; i < gene_seqs.size(); i++){
        blosum_distance = 0;

        for (int j=start_letter; j < end_letter; j++){
            switch (gene_seqs[i][j]){
                case 'A':
                    blosum_distance += blosum_itr(0,encoded_sequence[j]);
                    break;
                case 'C':
                    blosum_distance += blosum_itr(1,encoded_sequence[j]);
                    break;
                case 'D':
                    blosum_distance += blosum_itr(2,encoded_sequence[j]);
                    break;
                case 'E':
                    blosum_distance += blosum_itr(3,encoded_sequence[j]);
                    break;
                case 'F':
                    blosum_distance += blosum_itr(4,encoded_sequence[j]);
                    break;
                case 'G':
                    blosum_distance += blosum_itr(5,encoded_sequence[j]);
                    break;
                case 'H':
                    blosum_distance += blosum_itr(6,encoded_sequence[j]);
                    break;
                case 'I':
                    blosum_distance += blosum_itr(7,encoded_sequence[j]);
                    break;
                case 'K':
                    blosum_distance += blosum_itr(8,encoded_sequence[j]);
                    break;
                case 'L':
                    blosum_distance += blosum_itr(9,encoded_sequence[j]);
                    break;
                case 'M':
                    blosum_distance += blosum_itr(10,encoded_sequence[j]);
                    break;
                case 'N':
                    blosum_distance += blosum_itr(11,encoded_sequence[j]);
                    break;
                case 'O':
                    blosum_distance += blosum_itr(12,encoded_sequence[j]);
                    break;
                case 'P':
                    blosum_distance += blosum_itr(13,encoded_sequence[j]);
                    break;
                case 'Q':
                    blosum_distance += blosum_itr(14,encoded_sequence[j]);
                    break;
                case 'R':
                    blosum_distance += blosum_itr(15,encoded_sequence[j]);
                    break;
                case 'S':
                    blosum_distance += blosum_itr(16,encoded_sequence[j]);
                    break;
                case 'T':
                    blosum_distance += blosum_itr(17,encoded_sequence[j]);
                    break;
                case 'W':
                    blosum_distance += blosum_itr(18,encoded_sequence[j]);
                    break;
                case 'Y':
                    blosum_distance += blosum_itr(19,encoded_sequence[j]);
                    break;
                case '-':
                    break;
                default:
                    return INVALID_SEQUENCE;
                    break;

            }
        }

        if (blosum_distance < best_identity){
            closest_id = i;
            best_identity = blosum_distance;
        }
    }

    best_gene_name = gene_names[closest_id];
    return VALID_SEQUENCE;
}






void VJMatchCounter::prep_sequence(std::string &prepped_sequence, std::string &sequence,
        std::vector<std::string> &numbering){
    if (this->scheme == "imgt"){
        for (size_t i=0; i < sequence.length(); i++){
            if (this->essential_imgt_map.find(numbering[i]) !=
                    this->essential_imgt_map.end())
                prepped_sequence[std::stoi(numbering[i]) - 1] = sequence[i];
        }
    }
    else{
        std::vector<std::string> prepped_numbering;
        convert_numbering_to_imgt(numbering, prepped_numbering, this->scheme);
        for (size_t i=0; i < sequence.length(); i++){
            if (this->essential_imgt_map.find(prepped_numbering[i]) !=
                    this->essential_imgt_map.end())
                prepped_sequence[std::stoi(prepped_numbering[i]) - 1] = sequence[i];
        }
    }
}



std::string VJMatchCounter::get_vj_gene_sequence(std::string query_name,
        std::string species){

    size_t matchID;
    std::string output = "";

    if (query_name.length() < 4)
        return output;

    std::string dict_key = species + "_" + query_name.substr(0, 4);

    // We assume here that if the key is in names to positions, it is also in gene_seqs,
    // and that the vectors stored in both are the same length, because of the way
    // we initialized names to positions (see class constructor).
    if (this->names_to_positions.find(dict_key) != this->names_to_positions.end()){
        std::map<std::string, int> &name_to_position_map = this->names_to_positions[dict_key];
        if (name_to_position_map.find(query_name) != name_to_position_map.end()){
            matchID = name_to_position_map[query_name];
            output = this->gene_seqs[dict_key][matchID];
        }
    }
    return output;
}



// This function just returns the maps stored by the object,
// which are converted to python dicts if used through the
// wrapper. This is used only for testing.
std::tuple<std::map<std::string, std::vector<std::string>>,
  std::map<std::string, std::vector<std::string>>>  VJMatchCounter::get_seq_lists(){

    return std::tuple<std::map<std::string, std::vector<std::string>>,
        std::map<std::string, std::vector<std::string>>>  {this->gene_seqs, this->gene_names};
}
