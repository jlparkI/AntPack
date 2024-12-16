"""Contains the Dataset class which houses tools needed to load
and process data and store specific information about the loaded
data."""
from antpack import SingleChainAnnotator, PairedChainAnnotator
from .seq_data import SeqData


class Dataset:
    """Stores all of the important information about a dataset that has
    been loaded. Contains methods needed to load data from a fasta file
    and a csv."""

    def __init__(self, dtype = "paired", scheme = "imgt"):
        self.seq_data = []
        self.dtype = dtype
        self.scheme = scheme
        self.input_file = ""


    def load_fasta_file(self, fpath, seqtype = "amino"):
        """Loads a fasta file and processes the input data. If any errors
        are encountered, returns a message that is not None. Caller should
        not attempt to use the dataset / should discard and reload if
        error message is not None."""
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

        if len(metadatas) != len(sequences) or len(metadatas) == 0:
            return "No valid sequences were found."

        self.seq_data, err_message = self.align_sequences(sequences, metadatas)
        if err_message is not None:
            self.seq_data = []
            return err_message

        self.input_file = fpath
        return None



    def align_sequences(self, sequences, metadatas):
        """Aligns the list of input sequences and for each
        creates a SeqData object storing its metadata and
        other information."""
        seq_data = []

        if self.dtype == "paired":
            annotator = PairedChainAnnotator(scheme=self.scheme)
            for seq, mdata in zip(sequences, metadatas):
                try:
                    hchain, lchain = annotator.analyze_seq(seq)
                except:
                    return [], ("A fatal error occurred when loading the file; "
                                "please report.")

                if len(hchain[0]) == 0:
                    heavy_num, hlim = [], (0,0)
                else:
                    _, heavy_num, heavy_start, heavy_end = annotator.trim_alignment(seq,
                        hchain)
                    hlim = (heavy_start, heavy_end)
                if len(lchain[0]) == 0:
                    light_num, llim = [], (0,0)
                else:
                    _, light_num, light_start, light_end = annotator.trim_alignment(seq,
                        lchain)
                    hlim, llim = (heavy_start, heavy_end), (light_start, light_end)

                if hchain[3] == "" and lchain[3] == "":
                    errors = ""
                else:
                    errors = hchain[3] + " " + lchain[3]

                seq_data.append(SeqData(seq, mdata, heavy_num, hlim,
                    light_num, llim, errors = errors))

        else:
            annotator = SingleChainAnnotator(chains=["H", "K", "L"], scheme=self.scheme)
            for seq, mdata in zip(sequences, metadatas):

                try:
                    results = annotator.analyze_seq(seq)
                except:
                    return [], ("A fatal error occurred when loading the file; "
                                "please report.")
                _, num, sstart, send = annotator.trim_alignment(seq, results)
                lim = (sstart, send)

                if results[2] == "H":
                    seq_data.append(SeqData(seq, mdata, num, lim, errors = results[-1]))
                elif results[2] in ("K", "L"):
                    seq_data.append(SeqData(seq, mdata, light_numbering = num,
                        lchain_lim = lim, errors = results[-1]))
                else:
                    seq_data.append(SeqData(seq, mdata, errors = results[-1]))

        return seq_data, None



    def get_seq_data(self):
        """Returns the seq data objects stored by this dataset."""
        for seq_data in self.seq_data:
            yield seq_data


    def get_input_file(self):
        """Convenience file which returns the input filepath."""
        return self.input_file


    def get_num_seqs(self):
        """Convenience function which returns the number of sequences."""
        return len(self.seq_data)


    def get_property(self, seqnum, selected_property, selected_chain = "heavy"):
        """Returns the numbering and selected property for
        the indicated sequence (for plotting by other tools)."""
        return self.seq_data[seqnum].get_property(selected_property, selected_chain)
