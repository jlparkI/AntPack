"""Wraps tools for constructing a searchable database from
various input types."""
from ..utilities import read_fasta
from antpack.antpack_cpp_ext import DatabaseConstructionTool



def build_database_from_fasta(fasta_filepath,
        database_filepath, numbering_scheme="imgt",
        chain_type="single", receptor_type="mab",
        exclude_errs=True, pid_threshold=0.7):
    """Builds a database from a fasta file which may or may
    not be gzipped. The database is constructed so it can be
    searched in sublinear time and the sequence descriptions
    for each sequence in the fasta file are saved as metadata.

    Args:
        fasta_filepath (str): The location of the fasta file.
        database_filepath (str): The desired location and filename
            for the database.
        numbering_scheme (str): One of 'imgt', 'kabat', 'martin' or
            'aho'. If receptor_type is 'tcr', 'imgt' is the only
            allowed option.
        chain_type (str): One of 'single', 'paired', 'unknown'. If
            'paired' each sequence is assumed to be paired. If 'unknown'
            it is assumed each sequence MAY be paired and should be analyzed
            as paired just in case.
        receptor_type (str): One of 'mab', 'tcr'.
        exclude_errs (bool): If True, any sequences for which an alignment
            error is reported (e.g. missing residue at highly conserved
            position) are excluded. If sequences are paired and one of
            the paired sequences generates an error, that member of the pair
            only is excluded (the other member is retained).
        pid_threshold (float): A value between 0 and 1 for percent identity
            threshold. Sequences with pid less than this are assumed to
            represent alignment error and are excluded. If sequences are
            paired and one of the paired sequences has a pid less than this
            threshold, that member of the pair only is excluded (the other
            member is retained).

    Raises:
        RuntimeError: A RuntimeError is raised if invalid arguments are
            supplied.
    """
    db_construct_tool = DatabaseConstructionTool(database_filepath,
            numbering_scheme, chain_type, receptor_type,
            exclude_errs, pid_threshold)


    for seqinfo, seq in read_fasta(fasta_filepath):
        db_construct_tool.add_sequence(seqinfo, seq)

    db_construct_tool.finalize_db_construction()
