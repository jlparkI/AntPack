"""Simple tools for reading biological data files. This avoids
any need to use Biopython or other similar libraries as dependencies."""
import itertools
import gzip


def read_fasta(filepath):
    """Generator for reading fasta files. By default, reads
    one sequence at a time. Can take either gzipped or non-
    gzipped files as input. Assumes utf-8 encoding."""
    if filepath.endswith(".gz"):
        with gzip.open(filepath, "rt", encoding="utf-8") as fhandle:
            for header, line_group in itertools.groupby(fhandle,
                    lambda x: x.startswith(">")):
                if header:
                    seq_info = line_group.next()[1:]
                else:
                    sequence = ''.join(line.strip() for line in line_group)
                    yield seq_info, sequence
    else:
        with open(filepath, "r", encoding="utf-8") as fhandle:
            for header, line_group in itertools.groupby(fhandle,
                    lambda x: x.startswith(">")):
                if header:
                    seq_info = line_group.next()[1:]
                else:
                    sequence = ''.join(line.strip() for line in line_group)
                    yield seq_info, sequence
