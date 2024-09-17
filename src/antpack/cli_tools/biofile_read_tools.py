"""Simple tools for reading biological data files. This avoids
any need to use Biopython or other similar libraries as dependencies."""
import gzip


def read_fasta(filepath):
    """Generator for reading fasta files. By default, reads
    one sequence at a time. Can take either gzipped or non-
    gzipped files as input. Assumes utf-8 encoding."""
    if filepath.endswith(".gz"):
        fhandle = gzip.open(filepath, "rt", encoding="utf-8")
    else:
        fhandle = open(filepath, "r", encoding="utf-8")

    sequence, seqinfo = [], ""

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

    fhandle.close()
