"""This is the base class for Single and MultiChainAnnotator. Contains routines
shared by both."""
from antpack_cpp_ext import sort_position_codes_cpp




class AnnotatorBaseClass:

    def __init__(self, scheme):
        if scheme not in ["imgt", "martin", "kabat"]:
            raise ValueError("Unsupported scheme supplied.")

        self.scheme = scheme
        self.chain_maps = {"H":"H", "K":"L", "L":"L"}


    def sort_position_codes(self, position_code_list):
        """Convenience function that takes as input a list of position codes
        in any order and sorts them in a way that respects the properties of
        the numbering scheme. This is not as trivial as it appears, since
        in the IMGT scheme for example '111', '111A', '112B', '112A', '112" is
        correct ordering (but '86A', '87B', '87A' would not be). This is useful
        if you have a list of all of the unique position codes observed in a large
        dataset and would like to sort it so you can then encode all of these
        sequences as fixed-length arrays.

        Note that '-' is not a valid code so it is removed if it is present in the
        list passed to this function (you do not have to remove it yourself).

        Args:
            position_code_list (list): A list of strings. Each must be a valid position
                code for the scheme selected for this Annotator.

        Returns:
            sorted_code_list (list): A list of sorted position codes.

        Raises:
            RuntimeError: A runtime error is raised if invalid codes are supplied.
        """
        # Passing '-' to this function is a common mistake, so we check for it first.
        cleaned_codes = [c for c in position_code_list if c != "-"]
        sorted_codes, is_valid = sort_position_codes_cpp(cleaned_codes, self.scheme)
        if not is_valid:
            raise RuntimeError("One or more of the position codes supplied was rejected as "
                    "invalid by the sorter.")
        return sorted_codes


    def build_msa(self, sequences, annotations):
        """Builds a multiple sequence alignment using a list of sequences and a
        list of tuples output by analyze_seq or analyze_seqs (PairedChainAnnotator
        or SingleChainAnnotator).

        Args:
            sequences (list): A list of sequences.
            annotations (list): A list of tuples, each containing (numbering,
                percent_identity, chain_name, error_message). These tuples
                are what you will get as output if you pass sequences to
                the analyze_seq or analyze_seqs methods of SingleChainAnnotator
                or PairedChainAnnotator.

        Returns:
            position_codes (list): A list of position codes from the appropriate numbering
                scheme.
            aligned_seqs (list): A list of strings -- the input sequences all aligned
                to form an MSA.

        Raises:
            RuntimeError: A RuntimeError is raised if 1) sequences of different chaintypes
                are passed, 2) incorrect input is passed.
        """
        if not isinstance(annotations, list) or not isinstance(sequences, list):
            raise RuntimeError("sequences and annotations should both be lists.")

        if len(sequences) != len(annotations):
            raise RuntimeError("The length of sequences and annotations did not match.")

        position_codes = set()
        chain_types = set()

        for annotation in annotations:
            if len(annotation) != 4:
                raise RuntimeError("Each annotation should be a tuple with four elements: "
                        "numbering, percent identity, chain name, error message. "
                        "In general it is easiest to pass the output of analyze_seq() to this "
                        "function.")
            if not isinstance(annotation[0], list) or not isinstance(annotation[1], float) \
                or not isinstance(annotation[2], str) or not isinstance(annotation[3], str):
                raise RuntimeError("Each annotation should be a tuple with four elements: "
                        "numbering, percent identity, chain name, error message. "
                        "In general it is easiest to pass the output of analyze_seq() to this "
                        "function.")
            for position_code in annotation[0]:
                if position_code != "-":
                    position_codes.add(position_code)

            chain_types.add(self.chain_maps[annotation[2]])

        if len(chain_types) > 1:
            raise RuntimeError("You tried to build an MSA out of two different chain types!")

        position_codes = self.sort_position_codes(list(position_codes))

        position_dict = {k:i for i, k in enumerate(position_codes)}
        msa = [['-' for _ in position_codes] for s in sequences]

        for i, (sequence, annotation) in enumerate(zip(sequences, annotations)):
            for position, letter in zip(annotation[0], sequence):
                if position == "-":
                    continue
                msa[i][position_dict[position]] = letter

        msa = ["".join(m) for m in msa]
        return position_codes, msa


    def trim_alignment(self, sequence, alignment):
        """Takes as input a sequence and atuple produced by
        analyze_seq and trims off any gap regions at the end
        that result when there are amino acids on either end
        which are not part of the numbered variable region.
        The output from analyze_seq can be fed directly to this
        function.

        Args:
            sequence (str): The sequence corresponding to the numbering / alignment.
            alignment (tuple): A tuple of any length where the first element is
                the sequence_numbering list returned by analyze_seq.
                sequence_numbering must be of the same length as sequence.

        Returns:
            trimmed_sequence (str): The input sequence trimmed so that non-numbered AAs on
                either end have been removed.
            trimmed_numbering (list): The sequence_numbering from the input trimmed so that
                any '-' assignments on either end have been removed.
            exstart (int): The start of the extracted region represented by trimmed_sequence.
            exend (int): The end of the extracted region represented by trimmed_sequence.
                trimmed_sequence in other words is sequence[start:end].

        Raises:
            RuntimeError: A RuntimeError is raised if unexpected inputs are supplied.
        """
        if not isinstance(alignment, tuple):
            raise RuntimeError("Argument to trim_alignment must be a tuple.")

        numbering = alignment[0]
        if len(numbering) != len(sequence):
            raise RuntimeError("Sequence and numbering should be the same length to "
                    "use this function.")
        exstart = next((i for i in range(len(numbering)) if numbering[i] != '-'), 0)
        exend = next((i for i in range(exstart + 1, len(numbering)) if numbering[i] == '-'),
                len(numbering))
        trimmed_sequence = sequence[exstart:exend]
        trimmed_numbering = numbering[exstart:exend]

        return trimmed_sequence, trimmed_numbering, exstart, exend
