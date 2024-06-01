"""Contains the SingleChainAnnotator class, which provides the tools needed
to parse single chains (either heavy or light) or a sequence which may contain
a heavy or light chain. If you want to extract the heavy and light chain variable
regions from a sequence that probably contains both, use MultiChainAnnotator."""
from .annotator_base_class import AnnotatorBaseClass
from antpack_cpp_ext import validate_sequence, sort_position_codes_cpp




class SingleChainAnnotator(AnnotatorBaseClass):
    """This class contains the tools needed to parse and number
    either a single sequence or a list of sequences that are antibody
    heavy or light chain variable regions.
    """

    def __init__(self, chains=["H", "K", "L"], scheme="imgt",
            compress_init_gaps=False):
        """Class constructor.

        Args:
            chains (list): A list of chains. Each must be one of "H", "K", "L".
                If ["H", "K", "L"] (default), the annotator will automatically
                determine the most appropriate chain type for each input
                sequence.
            scheme (str): The numbering scheme. Must be one of "imgt",
                "martin", "kabat".
            compress_init_gaps (bool): If True, rearrange gaps in the first 5
                positions post-alignment so that gaps are at the beginning of
                the sequence wherever possible. Defaults to False.

        Raises:
            ValueError: A ValueError is raised if unacceptable inputs are
                supplied.
        """
        super().__init__(chains, scheme, compress_init_gaps)



    def analyze_seqs(self, sequences, get_region_labels = False):
        """Numbers and scores a list of input sequences.

        Args:
            sequences (list): A list of strings, each of which is a sequence
                containing the usual 20 amino acids.
            get_region_labels (bool): If True, get a list of labels: "-", "fmwk1",
                "cdr1", "fmwk2", "cdr2", "fmwk3", "cdr3", "fmwk4" to indicate
                to which region each numbered amino acid belongs. The cdr definitions
                that are used are the same as those for the numbering scheme (i.e.
                if using IMGT numbering IMGT CDR definitions are used). If False
                this list is not generated.

        Returns:
            sequence_results (list): A list of tuples of (sequence numbering, percent_identity,
                chain_name, error_message) if get_region_labels is False or
                (sequence_numbering, percent_identity, chain_name, error_message, region_labels)
                if get_region_labels is True. If no error was encountered, the error
                message is "". An alignment with low percent identity (e.g. < 0.85)
                may indicate a sequence that is not really an antibody, that contains
                a large deletion, or is not of the selected chain type.
        """
        if not isinstance(sequences, list):
            raise ValueError("sequences should be a list of strings.")


        sequence_results = []
        for sequence in sequences:
            if not validate_sequence(sequence):
                sequence_results.append(([], 0.0, "",
                    "Invalid sequence supplied -- nonstandard AAs"))
                continue
            results = [scoring_tool.align(sequence, get_region_labels) for
                    scoring_tool in self.scoring_tools]
            results = sorted(results, key=lambda x: x[1])
            if get_region_labels:
                sequence_results.append(results[-1])
            else:
                sequence_results.append(results[-1][:-1])

        return sequence_results


    def analyze_seq(self, sequence, get_region_labels = False):
        """Numbers and scores a single input sequence.

        Args:
            sequence (str): A string which is a sequence
                containing the usual 20 amino acids.
            get_region_labels (bool): If True, get a list of labels: "-", "fmwk1",
                "cdr1", "fmwk2", "cdr2", "fmwk3", "cdr3", "fmwk4" to indicate
                to which region each numbered amino acid belongs. The cdr definitions
                that are used are the same as those for the numbering scheme (i.e.
                if using IMGT numbering IMGT CDR definitions are used). If False
                this list is not generated.

        Returns:
            sequence_results (tuple): A tuple of (sequence numbering, percent_identity,
                chain_name, error_message) if get_region_labels is False or
                (sequence_numbering, percent_identity, chain_name, error_message, region_labels)
                if get_region_labels is True. If no error was encountered, the error
                message is "". An alignment with low percent identity (e.g. < 0.85)
                may indicate a sequence that is not really an antibody, that contains
                a large deletion, or is not of the selected chain type.
        """
        if not validate_sequence(sequence):
            return ([], 0.0, "", "Invalid sequence supplied -- nonstandard AAs")

        results = [scoring_tool.align(sequence, get_region_labels) for
                scoring_tool in self.scoring_tools]
        results = sorted(results, key=lambda x: x[1])
        sequence_results = results[-1]
        if get_region_labels:
            return sequence_results
        return sequence_results[:-1]



    def sort_position_codes(self, position_code_list, scheme):
        """Convenience function that takes as input a list of position codes
        in any order and sorts them in a way that respects the properties of
        the numbering scheme. This is not as trivial as it appears, since
        in the IMGT scheme for example '111', '111A', '112B', '112A', '112" is
        correct ordering (but '86A', '87B', '87A' would not be). This is useful
        if you have a list of all of the unique position codes observed in a large
        dataset and would like to sort it so you can then encode all of these
        sequences as fixed-length arrays.

        Args:
            position_code_list (list): A list of strings. Each must be a valid position
                code for the scheme in question.
            scheme (str): Must be a currently supported scheme.

        Returns:
            sorted_code_list (list): A list of sorted position codes.

        Raises:
            RuntimeError: A runtime error is raised if invalid codes are supplied.
        """
        if scheme not in ("imgt", "martin", "kabat"):
            raise RuntimeError("An invalid scheme was supplied.")

        sorted_codes, is_valid = sort_position_codes_cpp(position_code_list, scheme)
        if not is_valid:
            raise RuntimeError("One or more of the position codes supplied was rejected as "
                    "invalid by the sorter.")
        return sorted_codes
