"""Contains the MultiChainAnnotator class, which extracts variable region
heavy and light chains from an input sequence and returns them together with
appropriate numbering (if requested)."""
import copy
from .annotator_base_class import AnnotatorBaseClass
from antpack_cpp_ext import validate_sequence



class MultiChainAnnotator(AnnotatorBaseClass):
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
        super().__init__(scheme = scheme)

        self.chain_maps = {"H":"H", "K":"L", "L":"L"}


    def analyze_seq(self, sequence, identity_threshold = 0.8,
            max_chains_to_extract = 2):
        """Extracts the variable regions from a sequence. Also
        returns numbering if requested.

        Args:
            sequence (str): A string which is a sequence
                containing the usual 20 amino acids.
            identity_threshold (float): A minimum percent identity to assign a
                segment to a chain.
            max_chains_to_extract (int): Up to this number of chains can be
                extracted.

        Returns:
            extracted_seqs (list): A list of tuples of (subset, percent_identity,
                chain_name, error_message, (start, end)), where (start, end)
                mark the start and end of the extracted region. If no error was encountered,
                the error message is "". An alignment with low percent identity (e.g. < 0.85)
                may indicate a sequence that is not really an antibody, that contains
                a large deletion, or is not of the selected chain type.
        """
        if not validate_sequence(sequence):
            return [ ("", 0.0, "", "Invalid sequence supplied -- nonstandard AAs", (0,0)) ]

        extracted_seqs = []

        new_fragments = [sequence]
        k = 0

        while len(new_fragments) > 0 and k < max_chains_to_extract:
            fragments = copy.deepcopy(new_fragments)
            new_fragments = []

            for fragment in fragments:
                assignment = self._assign_fragment(fragment, identity_threshold)
                if assignment is None:
                    continue

                span = assignment[-1]

                extracted_seqs.append( (fragment[span[0]:span[1]], assignment[0], assignment[1],
                        assignment[2], assignment[3], span) )

                if span[0] > 0:
                    new_fragments.append(fragment[:span[0]])
                if span[1] < len(fragment):
                    new_fragments.append(fragment[span[1]:])
            k += 1

        return extracted_seqs



    def _assign_fragment(self, sequence, identity_threshold=0.8):
        """Assigns a fragment of a validated sequence to a chain type.

        Args:
            sequence (str): A string containing the standard amino acids.
            identity_threshold (float): A minimum percent identity to assign a
                segment to a chain.

        Returns:
            results (tuple): Either None or a tuple containing (numbering, percent identity,
                chain name, error message, (start, end). It is None if the chain cannot be
                assigned or if the minimum percent identity threshold is not reached.
        """
        if len(sequence) < 25:
            return None

        results = [scoring_tool.align(sequence, False) for scoring_tool in self.scoring_tools]
        results = sorted(results, key=lambda x: x[1])
        results = results[-1]
        if results[1] < identity_threshold:
            return None

        span = [0, len(sequence)]

        for i, imgt_num in enumerate(results[0]):
            if imgt_num != "-":
                span[0] = i
                break

        for i, imgt_num in reversed(list(enumerate(results[0]))):
            if imgt_num != "-":
                span[1] = i + 1
                break

        return (results[0][span[0]:span[1]], results[1], results[2], results[3], span)
