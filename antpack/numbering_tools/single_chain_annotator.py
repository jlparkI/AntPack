"""Contains the SingleChainAnnotator class, which provides the tools needed
to parse single chains (either heavy or light) or a sequence which may contain
a heavy or light chain. If you want to extract the heavy and light chain variable
regions from a sequence that probably contains both, use MultiChainAnnotator."""
from .annotator_base_class import AnnotatorBaseClass
from antpack_cpp_ext import validate_sequence




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
                the sequence wherever possible. This is more consistent with
                results from some other tools although it is debatable
                if this is more correct. Defaults to False.

        Raises:
            ValueError: A ValueError is raised if unacceptable inputs are
                supplied.
        """
        super().__init__(chains, scheme, compress_init_gaps)



    def analyze_seqs(self, sequences):
        """Numbers and scores a list of input sequences.

        Args:
            sequences (list): A list of strings, each of which is a sequence
                containing the usual 20 amino acids.

        Returns:
            sequence_results (list): A list of tuples of (sequence numbering, percent_identity,
                chain_name, error_message). If no error was encountered, the error
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
            results = [scoring_tool.align(sequence) for scoring_tool in self.scoring_tools]
            results = sorted(results, key=lambda x: x[1])
            sequence_results.append(results[-1])

        return sequence_results


    def analyze_seq(self, sequence):
        """Numbers and scores a single input sequence.

        Args:
            sequence (str): A string which is a sequence
                containing the usual 20 amino acids.

        Returns:
            sequence_results (tuple): A tuple of (sequence numbering, percent_identity,
                chain_name, error_message). If no error was encountered, the error
                message is "". An alignment with low percent identity (e.g. < 0.85)
                may indicate a sequence that is not really an antibody, that contains
                a large deletion, or is not of the selected chain type.
        """
        if not validate_sequence(sequence):
            return ([], 0.0, "", "Invalid sequence supplied -- nonstandard AAs")

        results = [scoring_tool.align(sequence) for scoring_tool in self.scoring_tools]
        results = sorted(results, key=lambda x: x[1])
        sequence_results = results[-1]
        return sequence_results
