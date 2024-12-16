"""Class for storing information about a specific sequence."""
from ..constants import properties



class SeqData:
    """Stores information about a specific sequence."""

    def __init__(self, raw_sequence, metadata, heavy_numbering = [],
            hchain_lim = (0,0), light_numbering = [],
            lchain_lim = (0,0), errors = ""):
        self.sequence = raw_sequence
        self.metadata = metadata
        self.heavy_numbering = heavy_numbering
        self.light_numbering = light_numbering
        self.hchain_lim = hchain_lim
        self.lchain_lim = lchain_lim
        self.errors = errors


    def get_metadata(self):
        """Returns the metadata (usually the sequence name)."""
        return self.metadata


    def get_heavy_chain(self):
        """Returns the subset of the input sequence that is heavy chain."""
        return self.sequence[self.hchain_lim[0]:self.hchain_lim[1]]

    def get_light_chain(self):
        """Returns the subset of the input sequence that is heavy chain."""
        return self.sequence[self.lchain_lim[0]:self.lchain_lim[1]]

    def get_heavy_numbering(self):
        """Returns the numbering (trimmed) for the heavy chain."""
        return self.heavy_numbering

    def get_light_numbering(self):
        """Returns the numbering (trimmed) for the heavy chain."""
        return self.light_numbering

    def get_errors(self):
        """Returns any annotation error messages associated with this sequence."""
        return self.errors


    def get_property(self, selected_property, chain = "heavy"):
        """Returns a selected property for each amino acid
        in this sequence for easy plotting by other tools."""
        if chain == "heavy":
            seqlim = self.hchain_lim
            nmbr = self.heavy_numbering
        elif chain == "light":
            seqlim = self.lchain_lim
            nmbr = self.light_numbering

        seq_subset = self.sequence[seqlim[0]:seqlim[1]]
        if len(seq_subset) == 0:
            return [], []
        if selected_property == "Hydrophobicity":
            props = [properties.AA_HYDROPHOBICITY[a] for a in
                        seq_subset]
            idx = list(range(len(props)))
        return idx, nmbr, props, self.metadata
