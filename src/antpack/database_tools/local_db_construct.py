"""Wraps tools for constructing a searchable database from
various input types."""
import os
import gc
import numpy as np
from ..utilities import read_fasta
from ..numbering_tools.cterm_finder import _load_nterm_kmers
from ..antpack_license import get_license_key_info
from ..utilities.vj_utilities import load_vj_gene_consensus_db
from antpack.antpack_cpp_ext import DatabaseConstructionTool



def build_database_from_fasta(fasta_filepaths:list,
        database_filepath:str, numbering_scheme:str="imgt",
        cdr_definition_scheme:str="imgt",
        sequence_type:str="single", receptor_type:str="mab",
        pid_threshold:float=0.7, user_memo:str="",
        reject_file:str = None,
        verbose:bool=True):
    """Builds a database from a fasta file which may or may
    not be gzipped. The database is constructed so it can be
    searched in sublinear time and the sequence descriptions
    for each sequence in the fasta file are saved as metadata.

    Args:
        fasta_filepath (str): A list of fasta filepath locations, which
            may or may not be gzipped.
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
            not meeting this threshold are rejected. If sequence_type is
            'paired' the sequences are still retained as long as one of
            the chains meets this threshold.
        user_memo (str): A string describing the purpose of the database / anything
            important you want your future self or other users to know about the
            contents. Will be saved as part of the database metadata.
        reject_file (str): Either None or a filepath. If None, any sequences that
            are rejected are silently ignored. If a filepath, rejected sequences
            are written to that filepath which is saved as a fasta file.
        verbose (bool): If True, print regular updates while running.

    Raises:
        RuntimeError: A RuntimeError is raised if invalid arguments are
            supplied.
    """
    if not isinstance(fasta_filepaths, list):
        raise RuntimeError("The fasta filepaths should be a list of filepaths. "
                "If you have only one filepath, you can enclose it in brackets "
                "to make a list.")
    if os.path.exists(database_filepath):
        raise RuntimeError("The database already exists.")

    os.makedirs(database_filepath)

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

    print("Estimating number of sequences.")
    nseqs = 0

    for fasta_filepath in fasta_filepaths:
        for i, (seqinfo, seq) in enumerate(read_fasta(fasta_filepath)):
            nseqs += 1

    if nseqs == 0:
        print("No sequences found.")
        return

    db_construct_tool = DatabaseConstructionTool(database_filepath,
            numbering_scheme, cdr_definition_scheme,
            sequence_type, receptor_type,
            pid_threshold, license_key, user_email,
            consensus_path, nterm_kmer_dict,
            vj_names, vj_seqs, blosum_matrix,
            user_memo, nseqs)

    print("Starting db construction.")
    db_construct_tool.open_transaction()

    if reject_file is None:
        for fasta_filepath in fasta_filepaths:
            for i, (seqinfo, seq) in enumerate(read_fasta(fasta_filepath)):
                db_construct_tool.add_sequence(seq, seqinfo, "", "", 0)
                if i % 10000 == 0:
                    db_construct_tool.close_transaction()
                    db_construct_tool.open_transaction()
                    if verbose:
                        print(f"{i} complete.")
    else:
        with open(reject_file, "w+", encoding="utf-8") as reject_handle:
            for fasta_filepath in fasta_filepaths:
                for i, (seqinfo, seq) in enumerate(read_fasta(fasta_filepath)):
                    rcode = db_construct_tool.add_sequence(seq, seqinfo, "", "", 0)
                    if rcode > 0:
                        reject_handle.write(f">{seqinfo}\n{seq}\n")
                    if i % 10000 == 0:
                        db_construct_tool.close_transaction()
                        db_construct_tool.open_transaction()
                        if verbose:
                            print(f"{i} complete.")

    if verbose:
        print("Now constructing database indices...")

    db_construct_tool.close_transaction()
    db_construct_tool.open_transaction()
    db_construct_tool.finalize_db_construction()
    db_construct_tool.close_transaction()
    gc.collect()
