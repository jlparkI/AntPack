"""Simple tools for reading biological data files. This avoids
any need to use Biopython or other similar libraries as dependencies."""
import gzip
import csv


def read_fasta(filepath):
    """Generator for reading fasta files. By default, reads
    one sequence at a time. Can take either gzipped or non-
    gzipped files as input. Assumes utf-8 encoding."""
    if filepath.endswith(".gz"):
        fopener = gzip.open
        fopen_code = "rt"
    else:
        fopener = open
        fopen_code = "r"

    sequence, seqinfo = [], ""

    with fopener(filepath, fopen_code, encoding="utf-8") as fhandle:
        for line in fhandle:
            if line.startswith(">"):
                if len(sequence) > 0 and len(seqinfo) > 0:
                    yield seqinfo, "".join(sequence)
                    sequence, seqinfo = [], line[1:].strip()

                elif len(sequence) == 0 and len(seqinfo) == 0:
                    seqinfo = line[1:].strip()

                else:
                    raise RuntimeError("Incorrect fasta file formatting or missing "
                            f"sequence id / name; please check input file {filepath}.")

            elif len(line.strip()) == 0:
                continue

            else:
                sequence.append(line.strip())

        if len(sequence) > 0 and len(seqinfo) > 0:
            yield seqinfo, "".join(sequence)
        elif len(sequence) != 0 or len(seqinfo) != 0:
            raise RuntimeError("Incorrect fasta file formatting or missing "
                    f"sequence id / name; please check input file {filepath}.")



def read_csv(filepath, skiprows=0, delimiter=','):
    """Generator for reading csv files. By default, reads
    one sequence at a time. Can take either gzipped or non-
    gzipped files as input. Assumes utf-8 encoding."""
    if filepath.endswith(".gz"):
        with gzip.open(filepath, "rt", encoding="utf-8") as fhandle:
            reader = csv.reader(fhandle, delimiter=delimiter)
            for i in range(skiprows):
                _ = next(reader)
            for row in reader:
                yield row
    else:
        with open(filepath, "rt", encoding="utf-8") as fhandle:
            reader = csv.reader(fhandle, delimiter=delimiter)
            for i in range(skiprows):
                _ = next(reader)
            for row in reader:
                yield row
