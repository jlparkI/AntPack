"""Wraps tools for constructing a searchable database from
various input types."""
import os
import gc
import numpy as np
from ..utilities import read_fasta, read_csv
from ..numbering_tools.cterm_finder import _load_nterm_kmers
from ..antpack_license import get_license_key_info
from ..utilities.vj_utilities import load_vj_gene_consensus_db
from antpack.antpack_cpp_ext import DatabaseConstructionTool



def build_database_from_fasta(fasta_filepaths:list,
        database_filepath:str, temp_file:str,
        numbering_scheme:str="imgt",
        cdr_definition_scheme:str="imgt",
        sequence_type:str="single", receptor_type:str="mab",
        pid_threshold:float=0.7, user_memo:str="",
        reject_file:str = None,
        verbose:bool=True):
    """Builds a database from a list of fasta files which may or may
    not be gzipped. The database is constructed so it can be
    searched quickly and the sequence descriptions for each sequence
    in the fasta file are saved as metadata.

    Args:
        fasta_filepath (str): A list of fasta filepath locations, which
            may or may not be gzipped.
        database_filepath (str): The desired location and filename
            for the database.
        temp_file (str): A path to a temporary file where the application
            can temporarily store data needed while initializing the
            database.
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
    if os.path.exists(temp_file):
        raise RuntimeError("The temporary filepath you supplied "
                           "is a file that already exists.")

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

    print("Estimating number of sequences.")
    nseqs = 0

    for fasta_filepath in fasta_filepaths:
        for i, (seqinfo, seq) in enumerate(read_fasta(fasta_filepath)):
            nseqs += 1

    if nseqs == 0:
        print("No sequences found.")
        return

    db_construct_tool = DatabaseConstructionTool(database_filepath,
            numbering_scheme, temp_file, cdr_definition_scheme,
            sequence_type, receptor_type,
            pid_threshold, license_key, user_email,
            consensus_path, nterm_kmer_dict,
            vj_names, vj_seqs, user_memo, nseqs)

    print("Starting db construction.")
    db_construct_tool.open_transaction()

    seqcount = 0

    if reject_file is None:
        for fasta_filepath in fasta_filepaths:
            for (seqinfo, seq) in read_fasta(fasta_filepath):
                db_construct_tool.add_sequence(seq, seqinfo)
                seqcount += 1
                if seqcount % 10000 == 0:
                    db_construct_tool.close_transaction()
                    db_construct_tool.open_transaction()
                    if verbose:
                        print(f"{seqcount} complete.")
    else:
        with open(reject_file, "w+", encoding="utf-8") as reject_handle:
            for fasta_filepath in fasta_filepaths:
                for (seqinfo, seq) in read_fasta(fasta_filepath):
                    rcode = db_construct_tool.add_sequence(seq, seqinfo)
                    if rcode > 0:
                        reject_handle.write(f">{seqinfo}\n{seq}\n")
                    if seqcount % 10000 == 0:
                        db_construct_tool.close_transaction()
                        db_construct_tool.open_transaction()
                        if verbose:
                            print(f"{seqcount} complete.")

    if verbose:
        print("Now constructing database indices...")

    db_construct_tool.finalize_db_construction(verbose)
    os.remove(temp_file)
    gc.collect()



def build_database_from_csv(csv_filepaths:list,
        database_filepath:str, temp_file:str,
        column_selections:dict,
        header_rows:int=1,
        numbering_scheme:str="imgt",
        cdr_definition_scheme:str="imgt",
        receptor_type:str="mab", pid_threshold:float=0.7,
        user_memo:str="", reject_file:str = None,
        verbose:bool=True):
    """Builds a database from a list of csv files which may or may
    not be gzipped. The database is constructed so it can be
    searched quickly. The csv files should already contain heavy
    and light chains split up into appropriate columns. The csv files
    can optionally contain additional columns with numbering or Vgenes
    -- if present these will save time for database construction since
    AntPack will not need to perform numbering or vgene identification.

    Args:
        csv_filepaths (list): A list of csv filepaths, which may or
            may not be gzipped.
        database_filepath (str): The desired location and filename
            for the database.
        temp_file (str): A path to a temporary file where the application
            can temporarily store data needed while initializing the
            database.
        column_selections (dict): A dictionary which should contain at least
            some of the following keys:

            * ``"heavy_chain"``: The number (from 0) of the column in each
              csv file containing the heavy chains. Do not include this
              if there are no heavy chains.
            * ``"light_chain"``: The number (from 0) of the column in each
              csv file containing the light chains. Do not include this
              if there are no light chains.
            * ``"heavy_chaintype"``: Required if there are heavy chains. The
              number (from 0) of the column in each csv file containing the
              heavy chain type. This will be H for antibodies but may be e.g.
              A for TCRs.
            * ``"light_chaintype"``: Required if there are light chains. The
              number (from 0) of the column in each csv file containing the
              light chain type. This could be L or K for antibodies or
              e.g. B for TCRs.
            * ``"heavy_numbering"``: The number (from 0) of the column (if any)
              containing numbering for the heavy chains. The numbering
              should be delimited using the '_' character. If you number
              your sequences using AntPack you can obtain this using
              '_'.join(antpack_numbering) for each sequence. If this is
              not supplied each heavy chain will be numbered by this
              function.
            * ``"light_numbering"``: The number (from 0) of the column (if any)
              containing numbering for the light chains. The numbering
              should be delimited using the '_' character. If you number
              your sequences using AntPack you can obtain this using
              '_'.join(antpack_numbering) for each sequence. If this is
              not supplied each light chain will be numbered by this
              function.
            * ``"heavy_vgene"``: The number (from 0) of the column (if any)
              containing vgene assignments for the heavy chains.
            * ``"light_vgene"``: The number (from 0) of the column (if any)
              containing vgene assignments for the light chains.
            * ``"species"``: The number (from 0) of the column (if any)
              containing species assignments (which are assumed to be the
              same for both chains in a row if heavy and light chains are
              present).
            * ``"metadata"``: The number (from 0) of the column (if any)
              containing metadata for each heavy / light chain or chain pair.

        header_rows (int): The number of header rows. Header rows are
            skipped.
        numbering_scheme (str): One of 'imgt', 'kabat', 'martin' or
            'aho'. If receptor_type is 'tcr', 'imgt' is the only
            allowed option.
        cdr_definition_scheme (str): One of 'imgt', 'kabat', 'martin',
            'aho' or 'north'. If receptor_type is 'tcr', 'imgt' only is allowed.
        receptor_type (str): One of 'mab', 'tcr'.
        pid_threshold (float): A value between 0 and 1 for percent identity
            threshold. If you have not supplied pre-written numbering,
            sequences are numbered during construction, and any not meeting
            this threshold are rejected.
        user_memo (str): A string describing the purpose of the database / anything
            important you want your future self or other users to know about the
            contents. Will be saved as part of the database metadata.
        reject_file (str): Either None or a filepath. If None, any sequences that
            are rejected are silently ignored. If a filepath, rejected sequences
            are written to that filepath which is saved as a csv file.
        verbose (bool): If True, print regular updates while running.

    Raises:
        RuntimeError: A RuntimeError is raised if invalid arguments are
            supplied.
    """
    if not isinstance(csv_filepaths, list):
        raise RuntimeError("The csv filepaths should be a list of filepaths. "
                "If you have only one filepath, you can enclose it in brackets "
                "to make a list.")
    if os.path.exists(database_filepath):
        raise RuntimeError("The database already exists.")
    if os.path.exists(temp_file):
        raise RuntimeError("The temporary filepath you supplied "
                           "is a file that already exists.")

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

    print("Estimating number of sequences.")
    nseqs = 0

    for csv_filepath in csv_filepaths:
        for i, _ in enumerate(read_csv(csv_filepath)):
            if i >= header_rows:
                break
        for i, row in enumerate(read_csv(csv_filepath)):
            nseqs += 1

    if nseqs == 0:
        print("No sequences found.")
        return

    db_construct_tool = DatabaseConstructionTool(database_filepath,
            numbering_scheme, temp_file, cdr_definition_scheme,
            "single", receptor_type,
            pid_threshold, license_key, user_email,
            consensus_path, nterm_kmer_dict,
            vj_names, vj_seqs, user_memo, nseqs)

    settings_list = []
    for expected_key in ["heavy_chain", "light_chain",
            "heavy_numbering", "light_numbering",
            "heavy_vgene", "light_vgene", "species",
            "metadata", "heavy_chaintype", "light_chaintype"]:
        if expected_key in column_selections:
            settings_list.append(column_selections[expected_key])
        else:
            settings_list.append(-1)

    if max(settings_list) < 0:
        raise RuntimeError("There must be at least some column "
                "selections supplied. Please refer to the docs "
                "for the build_database_from_csv function.")
    if max(settings_list[:2]) < 0:
        raise RuntimeError("Either a heavy column, a light column "
                "or both should be supplied under column_selections.")
    if settings_list[4] >= 0 or settings_list[5] >= 0:
        if settings_list[6] < 0:
            raise RuntimeError("If csv files contain either a "
                    "heavy or light vgene, they should contain "
                    "a species as well.")

    print("Starting db construction.")
    db_construct_tool.open_transaction()

    seqcount = 0

    if reject_file is None:
        for csv_filepath in csv_filepaths:
            for row in read_csv(csv_filepath, skiprows=header_rows):
                if len(row) == 0:
                    continue
                db_construct_tool.add_csv_sequence(row, settings_list)
                seqcount += 1
                if seqcount % 10000 == 0:
                    db_construct_tool.close_transaction()
                    db_construct_tool.open_transaction()
                    if verbose:
                        print(f"{seqcount} complete.")
    else:
        with open(reject_file, "w+", encoding="utf-8") as reject_handle:
            for csv_filepath in csv_filepaths:
                for row in read_csv(csv_filepath, skiprows=header_rows):
                    if len(row) == 0:
                        continue
                    rcode = db_construct_tool.add_csv_sequence(row, settings_list)
                    seqcount += 1
                    if rcode > 0:
                        reject_handle.write(f"{','.join(row)}\n")
                    if seqcount % 10000 == 0:
                        db_construct_tool.close_transaction()
                        db_construct_tool.open_transaction()
                        if verbose:
                            print(f"{seqcount} complete.")

    if verbose:
        print("Now constructing database indices...")

    db_construct_tool.finalize_db_construction(verbose)
    os.remove(temp_file)
    gc.collect()
