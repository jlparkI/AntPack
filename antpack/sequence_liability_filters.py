"""Contains the LiabilitySearchTool, which searches input
sequences for certain motifs known to correspond to
developability liabilities."""
import re
from ant_ext import validate_sequence
from .single_chain_annotator import SingleChainAnnotator


class LiabilitySearchTool:
    """This class provides a simple toolkit for searching
    for some common motifs which may correspond to
    possible developability liabilities.

    Attributes:
        filters (list): A list of filters to search for.
        filter_types (list): A list of filter types, indicates
            the type of liability corresponding to each filter.
            Must be of the same length as filters.
        aligner (SingleChainAnnotator): A tool for aligning input chains.
    """

    def __init__(self):
        """Class constructor."""
        input_filters = ["n[^p][st]", "n[gs]",
                "d[dghst]", "n[ahnt]",
                "[nd]p"]
        self.filters = [re.compile(input_filter, re.IGNORECASE)
                for input_filter in input_filters]
        self.filter_types = ["N-glycosylation", "Deamidation (elevated risk)",
                "Isomerization (elevated risk)", "Deamidation (moderate risk)",
                "pH-dependent hydrolysis (moderate risk)"]
        self.aligner = SingleChainAnnotator(chains = ["H", "K", "L"],
                    scheme = "imgt", compress_init_gaps = False)



    def analyze_seq(self, sequence):
        """Runs all the stored filters against the input sequence
        and determines which if any are matches and what type of
        liabilities correspond.

        Args:
            sequence (str): A string. May be either lower or upper case but
                must contain only valid amino acids.

        Returns:
            liabilities (list): A list of tuples. The first element of each
                tuple is a 2-tuple of (starting position, ending position)
                numbered with the start of the sequence as 0, indicating the
                start and end of the liability. The second element of each tuple
                is a string describing the type of liability found. If the list
                is empty, no liabilities were found. If sequence numbering fails,
                the list contains a single tuple indicating the cause of the failure.
                (This can occur if the expected cysteines in the chain are not present
                at the expected positions, which is also a liability.)
        """
        numbering, percent_ident, chain_name, err = self.aligner.analyze_seq(sequence)
        if err != "":
            return [((0,0), err)]

        liabilities = []

        for i, search_filter in enumerate(self.filters):
            for match_st in search_filter.finditer(sequence):
                liabilities.append( (match_st, self.filter_types[i]) )

        return liabilities
