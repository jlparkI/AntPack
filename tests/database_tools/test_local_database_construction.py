"""Test database construction for errors.
Use sqlite for testing and prototyping."""
import os
import gzip
import struct
import sqlite3
import pytest
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
    con = sqlite3.connect(db_filepath)
    cur = con.cursor()

    cur.execute("SELECT * FROM database_info;")
    metadata = cur.fetchone()
    assert metadata[0]==params["nmbr_scheme"]
    assert metadata[1]=="mab"
    assert metadata[2]==params["sequence_type"]
    assert metadata[3]==params["cdr_scheme"]
    assert metadata[4]==params["memo"]
    assert metadata[5]==0.7

    for i, row in enumerate(cur.execute("SELECT * FROM sequences;")):
        assert row[0]==seqs[i]
        assert row[1]==seqinfos[i]

    sca = SingleChainAnnotator(scheme=params["nmbr_scheme"])
    annotations = sca.analyze_seqs(seqs)
    vj_tool = VJGeneTool(scheme=params["nmbr_scheme"])

    # Check each set of chain tables for expected info.
    for table_code, chain_coding in enumerate([("H",),
                                    ("L", "K")]):
        aligned_seqs, cdrs, unusual_positions, cdr3_region_len, vgenes, vspecies, child_ids = \
                prep_seqs_for_comparison(params["nmbr_scheme"],
                    params["cdr_scheme"], chain_coding, seqs,
                    annotations, sca, vj_tool)

        for i, (child_id, cdr_grp) in enumerate(
                zip(child_ids, aligned_seqs)):
            cur.execute(f"SELECT * FROM _{table_code}_cdrs "
                        f"WHERE rowid = {child_id+1};")
            value = cur.fetchone()
            import pdb
            pdb.set_trace()
            vgene_code = get_vgene_code(vgenes[i], vspecies[i])
            assert value[0]==cdr_grp[0]
            assert vgene_code[0]==value[1]
            assert vgene_code[1]==value[2]
            assert vgene_code[2]==value[3]
            assert vgene_code[3]==value[4]
            assert unusual_positions[i]==value[5]
            assert value[10]==child_id+1

        kmer_profile = {}

        for j in range(cdr3_region_len - 1):
            kmer_to_child = {}

            for row in cur.execute("SELECT * FROM "
                                   f"_{table_code}_{j};"):
                if row[0] not in kmer_to_child:
                    kmer_to_child[row[0]] = set()
                kmer_to_child[row[0]].add( row[1] )

            kmer_profile = eval_nmbr_table_row_contents(
                    kmer_to_child, cdrs, child_ids,
                    kmer_profile, j)

        for row in cur.execute("SELECT * FROM "
                f"_{table_code}_column_diversity;"):
            assert row[1] in kmer_profile
            test_arr = np.frombuffer(row[-1], dtype=np.int64)
            test_arr = test_arr.reshape((row[2], row[3]))
            assert np.allclose(test_arr, kmer_profile[row[1]])



def test_low_quality_seqs(tmp_path):
    """Test what happens if we try to write low-quality
    sequences to the database."""
    return
    input_file = os.path.join(tmp_path, "temp_data_file.fa")
    db_path = os.path.join(tmp_path, "TEMP_DB")
    reject_file = os.path.join(tmp_path, "REJECTS")

    with open(input_file, "w+", encoding="utf-8") as fhandle:
        fhandle.write(">NAME\nAAAAAAAATTTTTTTTT\n")
        fhandle.write(">NAME\ntesting123\n")
        fhandle.write(">NAME\nEVQLEVQLEVQL\n")

    build_database_from_fasta([input_file], db_path,
        numbering_scheme="imgt",
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
    reject_file = os.path.join(tmp_path, "REJECTS2")

    with open(input_file, "w+") as fhandle:
        fhandle.write("AAAAAAAATTTTTTTTT,H,VJ123,humanoid\n")
        fhandle.write(">NAME,H,XYZ123,alien\n")
        fhandle.write("$$ALT,H,,unknown\n")

    build_database_from_csv([input_file], db_path,
            {"heavy_chain":0, "heavy_chaintype":1,
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
    return
    csv_filepath = os.path.join(get_test_data_filepath,
            "non_antibody_test_data", "tcr_simpletest.csv.gz")
    db_filepath = os.path.join(tmp_path, "TEMP_DB")

    expected_cdrs, expected_mainseq = [], []
    build_tcr_database_from_csv(
            [csv_filepath], db_filepath,
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

    con = sqlite3.connect(db_filepath)
    cur = con.cursor()

    for i, row in enumerate(cur.execute(
        "SELECT * FROM _0_cdrs;")):
        assert expected_cdrs[i]==row[0]

    for i, row in enumerate(cur.execute(
        "SELECT * FROM sequences;")):
        assert expected_mainseq[i]==row[0]




def eval_nmbr_table_row_contents(kmer_to_child,
        cdrs, child_ids, profile_counts, position,
        numbering_scheme="imgt"):
    """Tests the contents of the rows from the numbering
    table."""
    AAMAP = {k:i for i,k in enumerate("ACDEFGHIKLMNPQRSTVWY-")}
    letter_position_map = {"A":[48], "C":[32],
                           "D":[16], "E":[16],
                           "F":[0], "G":[48,32],
                           "H":[48,16], "I":[48,0],
                           "K":[32,16], "L":[32,0],
                           "M":[32], "N":[16,0],
                           "P":[32], "Q":[16,0],
                           "R":[32,16], "S":[48,32,16],
                           "T":[48,32,0], "V":[48,16,0],
                           "W":[32,16,0], "Y":[48,32,16,0],
                           "-":[], "X":[32]
                           }

    for i, cdr_group in enumerate(cdrs):
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

        bytestring = ['0' for j in range(64)]
        for j, letter in enumerate(cdr_extract):
            for offset_pos in letter_position_map[letter]:
                bytestring[j+offset_pos] = '1'

        unsigned_tag_filter = int(''.join(bytestring)[::-1], 2)
        tag_filter = int.from_bytes(
            unsigned_tag_filter.to_bytes(8, byteorder='little'),
            byteorder='little', signed=True)

        kmer = cdr[position:position+2]
        codeval = AAMAP[kmer[0]] * 21 + AAMAP[kmer[1]]
        profile_counts[cdr3len][position, codeval] += 1
        if kmer == "--":
            continue
        codeval = AAMAP[kmer[0]] * 21 * 50 + AAMAP[kmer[1]] * 50
        codeval += cdr3len

        assert codeval in kmer_to_child

        if tag_filter not in kmer_to_child[codeval]:
            import pdb
            pdb.set_trace()
        assert tag_filter in kmer_to_child[codeval]

    return profile_counts






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
        build_database_from_fasta([fasta_filepath], db_filepath,
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
