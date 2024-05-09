"""Contains the MultiChainAnnotator class, which extracts variable region
heavy and light chains from an input sequence and returns them together with
appropriate numbering (if requested)."""
from .single_chain_annotator import SingleChainAnnotator



class MultiChainAnnotator:
    """This class contains the tools needed to parse and number
    either a single sequence, a list of sequences, or a fasta file containing
    many sequences that belong to one of several possible chain types.
    """

    def __init__(self, scheme = "imgt"):
        """Class constructor.

        Args:
            scheme (str): The numbering scheme. Must be one of "imgt",
                "martin", "kabat".

        Raises:
            RuntimeError: A RuntimeError is raised if unacceptable inputs are
                supplied.
        """
        if scheme not in ["imgt", "martin", "kabat"]:
            raise RuntimeError("Unsupported scheme supplied.")

        self.light_aligner = SingleChainAnnotator(chains=["K", "L"],
                scheme=scheme)
        self.heavy_aligner = SingleChainAnnotator(chains=["H"], scheme=scheme)



    def analyze_seq(self, sequence):
        """Extracts the variable regions from a sequence. Also
        returns numbering if requested.

        Args:
            sequence (str): A string which is a sequence
                containing the usual 20 amino acids.

        Returns:
            extracted_seqs (list): A list of tuples of (subset, chain_type,
                percent_identity, error_message, (start, end)), where (start, end)
                mark the start and end of the extracted region. If no error was encountered,
                the error message is "". An alignment with low percent identity (e.g. < 0.85)
                may indicate a sequence that is not really an antibody, that contains
                a large deletion, or is not of the selected chain type.
        """
        return
