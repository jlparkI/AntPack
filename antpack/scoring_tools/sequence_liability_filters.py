"""Contains the LiabilitySearchTool, which searches input
sequences for certain motifs known to correspond to
developability liabilities."""
import re
from ..numbering_tools import SingleChainAnnotator
from .scoring_constants import allowed_imgt_pos as ahip


class LiabilitySearchTool:
    """This class provides a simple toolkit for searching
    for some common motifs which may correspond to
    possible developability liabilities.

    Attributes:
        filters (list): A list of tuples. The first element of each tuple
            is the compiled regex filter to search for. The second element is
            a bool indicating whether the filter should only be applied to CDRs
            or to the whole variable sequence. The third element is a string
            describing the liability to which the filter corresponds.
        aligner (SingleChainAnnotator): A tool for aligning input chains.
        cdr_sets (dict): A dictionary mapping chain name to a set
            of IMGT-numbered positions present in the CDRs of that chain.
    """

    def __init__(self):
        """Class constructor."""
        self.search_filters = [
                (re.compile("N[^P][ST]"), False, "N-glycosylation"),
                (re.compile("N[GS]"), True, "Deamidation (elevated risk)"),
                (re.compile("D[DGHST]"), True, "Isomerization (elevated risk)"),
                (re.compile("N[AHNT]"), True, "Deamidation (moderate risk)"),
                (re.compile("[ND]P"), True, "pH-dependent hydrolysis (moderate risk)"),
                (re.compile("TS"), True, "pH-dependent hydrolysis risk"),
                (re.compile("[MW]"), True, "Methionine / Tryptophan oxidation moderate risk"),
                (re.compile("[STK]N"), True, "Deamidation (low risk)"),
                ]

        self.aligner = SingleChainAnnotator(chains = ["H", "K", "L"],
                    scheme = "imgt", compress_init_gaps = False)

        self.cdr_sets = {"H": set(ahip.imgt_cdrs["H"]["1"] + ahip.imgt_cdrs["H"]["2"] +
                        ahip.imgt_cdrs["H"]["3"]),
                        "L": set(ahip.imgt_cdrs["L"]["1"] + ahip.imgt_cdrs["L"]["2"] +
                            ahip.imgt_cdrs["L"]["3"]),
                        "K": set(ahip.imgt_cdrs["L"]["1"] + ahip.imgt_cdrs["L"]["2"] +
                            ahip.imgt_cdrs["L"]["3"])
                        }



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
        numbering, _, chain_name, err = self.aligner.analyze_seq(sequence)

        # An alignment error may (often) indicate a liability. Return it with no
        # further analysis.
        if err != "":
            return [((0,0), err)]

        # Set up a mask indicating for each aa whether it is or is not part of a
        # CDR and of the variable region.
        cdr_mask = [k in self.cdr_sets[chain_name] for k in numbering]
        variable_mask = [True for k in numbering]
        for i, ntag in enumerate(numbering):
            if ntag != "-":
                break
            variable_mask[i] = False

        for i, ntag in reversed(list(enumerate(numbering))):
            if ntag != "-":
                break
            variable_mask[i] = False

        liabilities = []

        for search_filter, cdr_only, description in self.search_filters:
            for match_st in search_filter.finditer(sequence):
                span = match_st.span()
                if span[0] == span[1]:
                    continue
                if sum(variable_mask[span[0]:span[1]]) == 0:
                    continue
                if cdr_only:
                    if sum(cdr_mask[span[0]:span[1]]) == 0:
                        continue
                liabilities.append( (span, description) )

        # A last (very unusual) liability is the presence of cysteines OUTSIDE of
        # the standard positions. We check for this here.
        for i, (aa, imgt_pos) in enumerate(zip(sequence, numbering)):
            if aa == "C":
                if imgt_pos not in ["23", "104"]:
                    liabilities.append( ((i, i+1), "Unusual cysteine") )

        return liabilities
