"""Test database construction for errors.
Use sqlite for testing and prototyping."""
import os
import gzip
import sys
import struct
import pytest
import lmdb
import numpy as np
from antpack import (build_database_from_fasta,
        build_database_from_csv,
        build_tcr_database_from_csv,
        SingleChainAnnotator, VJGeneTool)
from antpack.utilities import read_fasta
from .database_utilities import get_vgene_code, setup_canonical_numbering
from ..conftest import (get_test_data_filepath, get_test_base_filepath)



def test_local_db_construct(build_local_mab_lmdb):
    """Check that local database construction yields a
    database containing the information we expect."""
    seqs, seqinfos, db_filepath, params = build_local_mab_lmdb
    env = lmdb.Environment(db_filepath, readonly=True, max_dbs=9)
    with env.begin() as txn:
        subdb = env.open_db(b"metadata", txn, create=False)
        cursor = lmdb.Cursor(subdb, txn)
        assert (cursor.get(b"numbering_scheme").decode()==
            params["nmbr_scheme"])
        assert cursor.get(b"receptor_type").decode()=="mab"
        assert (cursor.get(b"sequence_type").decode()==
                params["sequence_type"])
        assert (cursor.get(b"cdr_definition_scheme").decode()==
                params["cdr_scheme"])
        assert (cursor.get(b"user_memo").decode()==
                params["memo"])
        pid = struct.unpack('@d', cursor.get(b"pid_threshold"))[0]
        assert pid==0.7

        # Next, check that the sequences and seqinfo in
        # the db match what we loaded.
        subdb = env.open_db(b"main_seq_table", txn, create=False)
        cursor = lmdb.Cursor(subdb, txn)
        for i, seq in enumerate(seqs):
            key = struct.pack('@I', i)
            test_seq = cursor.get(key).decode()
            assert seq==test_seq

        subdb = env.open_db(b"seq_metadata_table", txn, create=False)
        cursor = lmdb.Cursor(subdb, txn)
        for i, seqinfo in enumerate(seqinfos):
            key = struct.pack('@I', i)
            assert seqinfo==cursor.get(key).decode()

        sca = SingleChainAnnotator(scheme=params["nmbr_scheme"])
        annotations = sca.analyze_seqs(seqs)
        vj_tool = VJGeneTool(scheme=params["nmbr_scheme"])

        # Check each set of chain tables for expected info.
        for chain, chain_coding in zip(["heavy", "light"],
                [("H",), ("L", "K")]):
            aligned_seqs, cdrs, unusual_positions, cdr3_region_len, vgenes, vspecies, child_ids = \
                    prep_seqs_for_comparison(params["nmbr_scheme"],
                        params["cdr_scheme"], chain_coding, seqs,
                        annotations, sca, vj_tool)

            env = lmdb.Environment(db_filepath, readonly=True,
                max_dbs=9+cdr3_region_len)
            txn = env.begin()
            subdb = env.open_db(f"{chain}_cdrs".encode(),
                txn, create=False)
            cursor = lmdb.Cursor(subdb, txn)
            for i, (child_id, cdr_grp) in enumerate(
                    zip(child_ids, aligned_seqs)):
                key = struct.pack('@I', child_id)
                value = cursor.get(key)
                vgene_code = get_vgene_code(vgenes[i], vspecies[i])
                assert value[:-8].decode()==cdr_grp[0]
                assert vgene_code[0] == value[-8]
                assert vgene_code[1] == value[-7]
                assert vgene_code[2] == value[-6]
                assert vgene_code[3] == value[-5]
                assert value[-1] == unusual_positions[i]

            if chain == "heavy":
                table_code = 0
            else:
                table_code = 1

            kmer_profile = {}

            for j in range(cdr3_region_len - 1):
                subdb = env.open_db(f"{j}_{table_code}_dimers".encode(),
                    txn, create=False, dupsort=True)
                cursor = lmdb.Cursor(subdb, txn)
                kmer_key_list = list(cursor.iternext(values=False))
                kmer_to_child = {int.from_bytes(kmer[:4], sys.byteorder,
                signed=False):set() for kmer in kmer_key_list}
                for key in kmer_to_child.keys():
                    value = cursor.get(struct.pack('@i', key))
                    kmer_to_child[key].add(convert_value_to_bin(value))
                    for value in cursor.iternext_dup(keys=False):
                        kmer_to_child[key].add(convert_value_to_bin(value))

                kmer_profile = eval_nmbr_table_row_contents(
                        kmer_to_child, cdrs, child_ids,
                        kmer_profile, j)

            subdb = env.open_db(f"{chain}_diversity".encode(),
                    txn, create=False)
            cursor = lmdb.Cursor(subdb, txn)

            cursor = lmdb.Cursor(subdb, txn)
            for key, value in cursor.iternext():
                int_key = struct.unpack('@i', key)[0]
                assert int_key in kmer_profile

                assert len(value) % 8 == 0
                assert (len(value) / 8 ==
                        kmer_profile[int_key].flatten().shape[0])

                loaded_kmer_table = \
                        np.frombuffer(value, dtype=np.int64)
                assert np.allclose(loaded_kmer_table,
                    kmer_profile[int_key].flatten())
                kmer_profile[int_key] = None

            for _, value in kmer_profile.items():
                assert value is None



def test_low_quality_seqs(tmp_path):
    """Test what happens if we try to write low-quality
    sequences to the database."""
    input_file = os.path.join(tmp_path, "temp_data_file.fa")
    db_path = os.path.join(tmp_path, "TEMP_DB")
    temp_file = os.path.join(tmp_path, "TEMP_FILE")
    reject_file = os.path.join(tmp_path, "REJECTS")

    with open(input_file, "w+", encoding="utf-8") as fhandle:
        fhandle.write(">NAME\nAAAAAAAATTTTTTTTT\n")
        fhandle.write(">NAME\ntesting123\n")
        fhandle.write(">NAME\nEVQLEVQLEVQL\n")

    build_database_from_fasta([input_file], db_path,
        temp_file, numbering_scheme="imgt",
        cdr_definition_scheme="imgt",
        sequence_type="paired", receptor_type="mab",
        pid_threshold=0.7, user_memo="test",
        reject_file=reject_file)

    # There should be no error in writing the db;
    # we expect that the sequences that failed the
    # threshold test will be written to the reject file,
    # and the relevant database tables will be empty.
    with open(reject_file, "r", encoding="utf-8") as reject_handle:
        with open(input_file, "r", encoding="utf-8") as \
                fhandle:
            for line1, line2 in zip(reject_handle, fhandle):
                assert line1.strip().split('\t')[0]==line2.strip() 

    input_file = os.path.join(tmp_path, "temp_data_file.csv")
    db_path = os.path.join(tmp_path, "TEMP_DB2")
    temp_file = os.path.join(tmp_path, "TEMP_FILE2")
    reject_file = os.path.join(tmp_path, "REJECTS2")

    with open(input_file, "w+") as fhandle:
        fhandle.write("AAAAAAAATTTTTTTTT,H,VJ123,humanoid\n")
        fhandle.write(">NAME,H,XYZ123,alien\n")
        fhandle.write("$$ALT,H,,unknown\n")

    build_database_from_csv([input_file], db_path,
            temp_file, {"heavy_chain":0, "heavy_chaintype":1,
                "heavy_vgene":2, "species":3},
            numbering_scheme='imgt',
            cdr_definition_scheme='imgt',
            receptor_type="mab", pid_threshold=0.7,
            user_memo="", reject_file=reject_file,
            header_rows=0)

    # There should be no error in writing the db;
    # we expect that the sequences that failed the
    # threshold test will be written to the reject file,
    # and the relevant database tables will be empty.
    with open(reject_file, "r", encoding="utf-8") as reject_handle:
        with open(input_file, "r", encoding="utf-8") as \
                fhandle:
            for line1, line2 in zip(reject_handle, fhandle):
                assert line1==line2


def test_tcr_fmt_loading(get_test_data_filepath, tmp_path):
    """Check that tcr standard format data is loaded
    as expected."""
    csv_filepath = os.path.join(get_test_data_filepath,
            "non_antibody_test_data", "tcr_simpletest.csv.gz")
    db_filepath = os.path.join(tmp_path, "TEMP_DB")
    temp_filepath = os.path.join(tmp_path, "TEMP_FILE")

    expected_cdrs, expected_mainseq = [], []
    build_tcr_database_from_csv(
            [csv_filepath], db_filepath, temp_filepath,
            {"beta_cdr3":0, "beta_vgene":1,
             "beta_jgene":2, "species":4},
            user_memo="", reject_file=None)
    with gzip.open(csv_filepath, "rt") as fhandle:
        _ = fhandle.readline()
        for line in fhandle:
            elements = line.split(",")
            expected_mainseq.append(elements[0] +
                "," + elements[1] + "," +
                elements[2])
            expected_cdrs.append(elements[3])

    env = lmdb.Environment(db_filepath, readonly=True, max_dbs=10)
    with env.begin() as txn:
        subdb = env.open_db(b"heavy_cdrs", txn, create=False)
        cursor = lmdb.Cursor(subdb, txn)
        for i, expected_cdr in enumerate(expected_cdrs):
            key = struct.pack('@I', i)
            test_seq = cursor.get(key)
            assert expected_cdr==test_seq[:-7].decode()

    with env.begin() as txn:
        subdb = env.open_db(b"main_seq_table", txn, create=False)
        cursor = lmdb.Cursor(subdb, txn)
        for i, expected_info in enumerate(expected_mainseq):
            key = struct.pack('@I', i)
            test_seq = cursor.get(key).decode()
            assert expected_info==test_seq


def eval_nmbr_table_row_contents(kmer_to_child,
        cdrs, child_ids, profile_counts, position,
        numbering_scheme="imgt"):
    """Tests the contents of the rows from the numbering
    table."""
    AAMAP = {k:i for i,k in enumerate("ACDEFGHIKLMNPQRSTVWY-")}

    for i, cdr_group in enumerate(cdrs):
        # TODO: For now only considering cdr3 for indexing.
        cdr = cdr_group[2]
        cdr3len = len(cdr.replace('-', ''))

        if cdr3len not in profile_counts:
            profile_counts[cdr3len] = \
                    np.zeros(( len(cdr) - 1, 21*21))

        # Construct the augmented child id.
        if numbering_scheme in ("imgt", "aho"):
            cdr_extract = cdr[:8] + cdr[-8:]
        else:
            cdr_extract = cdr[:16]

        bytestring = ['0' for j in range(32)]
        for j, letter in enumerate(cdr_extract):
            if letter in ('C', 'G', 'P', 'H', 'M', 'X'):
                bytestring[j] = '1'
            elif letter in ('S', 'T', 'N', 'Q', 'Y'):
                bytestring[j+16] = '1'
            elif letter in ('R', 'K', 'D', 'E'):
                bytestring[j] = '1'
                bytestring[j+16] = '1'

        bytestring = (child_ids[i], int(''.join(bytestring)[::-1], 2))

        kmer = cdr[position:position+2]
        codeval = AAMAP[kmer[0]] * 21 + AAMAP[kmer[1]]
        profile_counts[cdr3len][position, codeval] += 1
        if kmer == "--":
            continue
        codeval = position * 441 * 441 + cdr3len * 441 + \
                codeval

        assert codeval in kmer_to_child

        assert bytestring in kmer_to_child[codeval]
        kmer_to_child[codeval].remove(bytestring)

    return profile_counts





def convert_value_to_bin(value):
    """Converts an input value to an
    integer and a binary string."""
    output_int = int.from_bytes(value[4:], sys.byteorder,
            signed=False)
    binstring = int.from_bytes(value[:4], sys.byteorder, signed=False)
    return (output_int, binstring)



def prep_seqs_for_comparison(numbering_scheme,
        cdr_scheme, chain_type, sequences,
        annotations, annotator,
        vj_tool):
    """Prep sequences which are of the specified chain
    type for analysis by aligning to the template,
    extracting cdrs etc."""
    selected_seqs = [(s, a) for (s, a) in zip(
        sequences, annotations) if a[2] in chain_type]
    child_ids = [i for i,a in enumerate(annotations)
            if a[2] in chain_type]
    vgenes, vjspecies = [], []

    for selected_seq in selected_seqs:
        vgene, _, _, _, species = vj_tool.assign_vj_genes(
                selected_seq[1], selected_seq[0], "unknown",
                "identity")
        vgenes.append(vgene)
        vjspecies.append(species)

    template_aligner, canon_nmbr = setup_canonical_numbering(
            numbering_scheme, cdr_scheme, chain_type[0])
    recognized_positions = set(canon_nmbr)
    recognized_positions.add('-')
    cdr_labels = annotator.assign_cdr_labels(canon_nmbr,
            chain_type[0], cdr_scheme)

    cdr3_codes = [a for (a,c) in zip(canon_nmbr, cdr_labels)
            if c == "cdr3"]
    cdr3_region_len = len(cdr3_codes)

    aligned_seqs, cdrs, unusual_positions = \
            [], [], []

    for selected_seq in selected_seqs:
        aligned_seq = template_aligner.align_sequence(
                selected_seq[0], selected_seq[1][0], False)
        unusual_positions.append( len([a for a in selected_seq[1][0]
                if a not in recognized_positions]) )
        cdr1 = ''.join([a for (a,c) in zip(aligned_seq,
            cdr_labels) if c == "cdr1"])
        cdr2 = ''.join([a for (a,c) in zip(aligned_seq,
            cdr_labels) if c == "cdr2"])
        cdr3 = ''.join([a for (a,c) in zip(aligned_seq,
            cdr_labels) if c == "cdr3"])
        cdrs.append((cdr1, cdr2, cdr3))
        aligned_seqs.append((cdr1 + cdr2 + cdr3, cdr3))

    return aligned_seqs, cdrs, unusual_positions, \
        cdr3_region_len, vgenes, vjspecies, child_ids




@pytest.fixture(params=[
    {"filepath":"addtnl_test_data.fasta.gz",
     "nmbr_scheme":"imgt", "cdr_scheme":"imgt",
     "sequence_type":"single", "memo":"testing123"},
    {"filepath":"addtnl_test_data.fasta.gz",
     "nmbr_scheme":"imgt", "cdr_scheme":"north",
     "sequence_type":"unknown", "memo":"testing123"},
    {"filepath":"db_csv_load_testing.csv.gz",
     "nmbr_scheme":"imgt", "cdr_scheme":"imgt",
     "sequence_type":"single", "memo":"testing123"},
    ])
def build_local_mab_lmdb(tmp_path, get_test_data_filepath,
        request):
    """Builds a local db using a set of arguments
    appropriate for testing LMDB-based local dbs
    and returns all the info the test function will need
    to evaluate this local db."""
    db_filepath = os.path.join(tmp_path, "TEMP_DB")
    seqs, seqinfos = [], []

    if "csv" in request.param["filepath"]:
        csv_filepath = os.path.join(get_test_data_filepath,
                        request.param["filepath"])
        # Using hardcoded column selections here (since the
        # test data is always the same)
        build_database_from_csv([csv_filepath], db_filepath,
                os.path.join(tmp_path, "TEMP_FILE"),
                {"heavy_chain":0, "light_chain":3,
                    "heavy_chaintype":2, "light_chaintype":5,
                    "heavy_numbering":6, "light_numbering":7},
                numbering_scheme=request.param["nmbr_scheme"],
                cdr_definition_scheme=request.param["cdr_scheme"],
                receptor_type="mab", pid_threshold=0.7,
                user_memo=request.param["memo"], reject_file=None)
        with gzip.open(csv_filepath, "rt") as fhandle:
            _ = fhandle.readline()
            for line in fhandle:
                seqinfos.append("")
                elements = line.split(",")
                if len(elements[0]) == 0:
                    seqs.append(elements[3])
                else:
                    seqs.append(elements[0])

    else:
        fasta_filepath = os.path.join(get_test_data_filepath,
                    request.param["filepath"])
        build_database_from_fasta([fasta_filepath],
            db_filepath, os.path.join(tmp_path, "TEMP_FILE"),
            numbering_scheme=request.param["nmbr_scheme"],
            cdr_definition_scheme=request.param["cdr_scheme"],
            sequence_type=request.param["sequence_type"],
            receptor_type="mab",
            pid_threshold=0.7, user_memo=request.param["memo"],
            reject_file=None)
        for seqinfo, seq in read_fasta(fasta_filepath):
            seqs.append(seq)
            seqinfos.append(seqinfo)

    return seqs, seqinfos, db_filepath, request.param
