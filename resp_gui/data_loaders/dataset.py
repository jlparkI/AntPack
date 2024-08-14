"""Contains the Dataset class which houses tools needed to load
and process data and store specific information about the loaded
data."""
from antpack import SingleChainAnnotator, PairedChainAnnotator



class Dataset:
    """Stores all of the important information about a dataset that has
    been loaded. Contains methods needed to load data from a fasta file
    and a csv."""

    def __init__(self):
        self.light_metadata, self.light_annotation = [], []
        self.heavy_metadata, self.heavy_annotation = [], []
        self.heavy_seqs, self.light_seqs = [], []



    def load_fasta_file(self, fpath, scheme = "imgt",
            dtype = "paired", seqtype = "amino"):
        """Loads a fasta file and processes the input data. If any errors
        are encountered, returns a message that is not None. Caller should
        not attempt to use the dataset / should discard and reload if
        error message is not None."""
        if dtype not in ("paired", "single"):
            return "Unknown data type passed."

        if seqtype not in ("amino"):
            return "Currently only amino acid data is supported."

        if not fpath.endswith(".fa") and not fpath.endswith(".fasta"):
            return "File must be a fasta file."

        try:
            with open(fpath, "r", encoding="utf-8") as _:
                pass
        except:
            return "Invalid filepath specified."

        metadatas, sequences = [], []

        with open(fpath, "r", encoding="utf-8") as fhandle:
            metadata, sequence = None, []

            for l in fhandle:
                if len(l.strip()) == 0:
                    continue
                if l.startswith(">"):
                    if metadata is not None:
                        metadatas.append(metadata)
                        sequences.append("".join(sequence))

                    metadata = l[1:].strip()
                    sequence = []
                    continue

                sequence.append(l.strip())

            if metadata is not None:
                metadatas.append(metadata)
                sequences.append("".join(sequence))

        if len(metadata) != len(sequences) or len(metadata) == 0:
            return "No valid sequences were found."

        if dtype == "paired":
            annotator = PairedChainAnnotator(scheme=scheme)
            for seq, mdata in zip(sequences, metadatas):
                hchain, lchain = annotator.analyze_seq(seq)
                self.heavy_annotation.append(hchain)
                self.light_annotation.append(lchain)
                self.heavy_metadata.append(mdata)
                self.light_metadata.append(mdata)
        else:
            annotator = SingleChainAnnotator(chains=["H", "K", "L"], scheme=scheme)
            for seq, mdata in zip(sequences, metadatas):

                hchain, lchain = annotator.analyze_seq(seq)
                self.heavy_chains.append(hchain)
                self.light_chains.append(lchain)
                self.heavy_metadata.append(mdata)
                self.light_metadata.append(mdata)

        return None
