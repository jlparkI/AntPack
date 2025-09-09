"""Wraps tools for constructing a searchable database from
various input types."""
import os
import numpy as np
from ..utilities import read_fasta
from ..numbering_tools.cterm_finder import _load_nterm_kmers
from ..antpack_license import get_license_key_info
from ..utilities.vj_utilities import load_vj_gene_consensus_db
from antpack.antpack_cpp_ext import DatabaseConstructionTool



def build_database_from_fasta(fasta_filepath,
        database_filepath, numbering_scheme="imgt",
        cdr_definition_scheme="imgt",
        sequence_type="single", receptor_type="mab",
        pid_threshold=0.7, user_memo="", verbose=True):
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
        cdr_definition_scheme (str): One of 'imgt', 'kabat', 'martin',
            'aho' or 'north'. If receptor_type is 'tcr', 'imgt' only is allowed.
        sequence_type (str): One of 'single', 'paired', 'unknown'. If
            'paired' each sequence is assumed to be paired. If 'unknown'
            it is assumed each sequence MAY be paired and should be analyzed
            as paired just in case.
        receptor_type (str): One of 'mab', 'tcr'.
        pid_threshold (float): A value between 0 and 1 for percent identity
            threshold. If sequence_type is 'single' or 'unknown', sequences
            not meeting this threshold are excluded. If sequence_type is
            'paired' the sequences are still retained as long as one of
            the chains meets this threshold.
        user_memo (str): A string describing the purpose of the database / anything
            important you want your future self or other users to know about the
            contents. Will be saved as part of the database metadata.
        verbose (bool): If True, print regular updates while running.

    Raises:
        RuntimeError: A RuntimeError is raised if invalid arguments are
            supplied.
    """
    if os.path.exists(database_filepath):
        raise RuntimeError("The database already exists.")
    license_key, user_email = get_license_key_info()
    project_path = os.path.join(os.path.abspath(os.path.dirname(__file__)),
            "..")
    consensus_path = os.path.join(project_path,
            "numbering_tools", "consensus_data")
    nterm_kmer_dict = _load_nterm_kmers()

    vj_db_path = os.path.join(project_path, "vj_tools", "consensus_data")
    vj_names, vj_seqs, _ = load_vj_gene_consensus_db(os.getcwd(),
            vj_db_path, "imgt")

    blosum_matrix = np.load(os.path.join(project_path,
        "numbering_tools", "consensus_data", "mabs",
        "blosum_matrix.npy")).astype(np.float64)

    db_construct_tool = DatabaseConstructionTool(database_filepath,
            numbering_scheme, cdr_definition_scheme,
            sequence_type, receptor_type,
            pid_threshold, license_key, user_email,
            consensus_path, nterm_kmer_dict,
            vj_names, vj_seqs, blosum_matrix,
            user_memo)


    db_construct_tool.open_transaction()

    for i, (seqinfo, seq) in enumerate(read_fasta(fasta_filepath)):
        db_construct_tool.add_sequence_new_db(seq, seqinfo)
        if i % 1000 == 0:
            db_construct_tool.close_transaction()
            db_construct_tool.open_transaction()
            if verbose:
                print(f"{i} complete.")

    db_construct_tool.close_transaction()
    db_construct_tool.open_transaction()
    db_construct_tool.finalize_db_construction()
    db_construct_tool.close_transaction()
