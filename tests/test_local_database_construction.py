"""Test database construction for errors.
Use sqlite for testing and prototyping."""
import os
import gzip
import sys
import struct
import shutil
import unittest
import lmdb
import numpy as np
from antpack import (build_database_from_fasta,
        build_database_from_csv,
        build_tcr_database_from_csv,
        SingleChainAnnotator, VJGeneTool)
from antpack.utilities import read_fasta
from antpack.antpack_cpp_ext import (SequenceTemplateAligner,
        return_imgt_canonical_numbering_cpp)



class TestLocalDBConstruction(unittest.TestCase):



    def test_local_db_construct(self):
        """Check that local database construction yields a
        database containing the information we expect."""
        current_dir = os.path.abspath(os.path.dirname(
            __file__))
        fasta_filepath = os.path.join(current_dir,
                "test_data", "addtnl_test_data.fasta.gz")
        csv_filepath = os.path.join(current_dir,
                "test_data", "db_csv_load_testing.csv.gz")

        try:
            shutil.rmtree("TEMP_DB")
        except:
            pass

        for use_csv in [True, False]:
            for (nmbr_scheme, cdr_scheme) in [('imgt', 'north'),
                ('imgt', 'kabat')]:
                if cdr_scheme == "north":
                    sequence_type = "single"
                    memo = ""
                elif nmbr_scheme == "kabat":
                    memo = "testing123"
                    sequence_type = "unknown"
                else:
                    memo="t"
                    sequence_type = "paired"
                if use_csv:
                    sequence_type = "single"

                seqs, seqinfos, vgenes = [], [], []

                if use_csv:
                    # Using hardcoded column selections here (since the
                    # test data is always the same)
                    build_database_from_csv([csv_filepath], "TEMP_DB",
                            "TEMP_FILE",
                            {"heavy_chain":0, "light_chain":3,
                                "heavy_chaintype":2, "light_chaintype":5,
                                "heavy_numbering":6, "light_numbering":7},
                            numbering_scheme=nmbr_scheme,
                            cdr_definition_scheme=cdr_scheme,
                            receptor_type="mab", pid_threshold=0.7,
                            user_memo=memo, reject_file=None)
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
                    build_database_from_fasta([fasta_filepath],
                        "TEMP_DB", "TEMP_FILE",
                        numbering_scheme=nmbr_scheme,
                        cdr_definition_scheme=cdr_scheme,
                        sequence_type=sequence_type, receptor_type="mab",
                        pid_threshold=0.7, user_memo=memo,
                        reject_file=None)
                    for seqinfo, seq in read_fasta(fasta_filepath):
                        seqs.append(seq)
                        seqinfos.append(seqinfo)


                env = lmdb.Environment("TEMP_DB", readonly=True,
                        max_dbs=9)
                with env.begin() as txn:
                    subdb = env.open_db(b"metadata", txn, create=False)
                    cursor = lmdb.Cursor(subdb, txn)
                    self.assertTrue(cursor.get(b"numbering_scheme").decode()==
                            nmbr_scheme)
                    self.assertTrue(cursor.get(b"receptor_type").decode()==
                            "mab")
                    self.assertTrue(cursor.get(b"sequence_type").decode()==
                            sequence_type)
                    self.assertTrue(cursor.get(b"cdr_definition_scheme").decode()==
                            cdr_scheme)
                    self.assertTrue(cursor.get(b"user_memo").decode()==
                            memo)
                    pid = struct.unpack('@d',
                            cursor.get(b"pid_threshold"))[0]
                    self.assertTrue(pid==0.7)

                # Next, check that the sequences and seqinfo in
                # the db match what we loaded.
                txn = env.begin()
                subdb = env.open_db(b"main_seq_table",
                        txn, create=False)
                cursor = lmdb.Cursor(subdb, txn)
                for i, seq in enumerate(seqs):
                    key = struct.pack('@I', i)
                    test_seq = cursor.get(key).decode()
                    self.assertTrue(seq==test_seq)

                subdb = env.open_db(b"seq_metadata_table",
                        txn, create=False)
                cursor = lmdb.Cursor(subdb, txn)
                for i, seqinfo in enumerate(seqinfos):
                    key = struct.pack('@I', i)
                    self.assertTrue(seqinfo==cursor.get(key).decode())

                sca = SingleChainAnnotator(scheme=nmbr_scheme)
                annotations = sca.analyze_seqs(seqs)
                vj_tool = VJGeneTool(scheme=nmbr_scheme)

                # Check each set of chain tables for expected info.
                for chain, chain_coding in zip(["heavy", "light"],
                        [("H",), ("L", "K")]):
                    aligned_seqs, cdrs, _, cdr3_region_len, vgenes, vspecies, child_ids = \
                            prep_seqs_for_comparison(nmbr_scheme,
                            cdr_scheme, chain_coding, seqs,
                            annotations, sca, vj_tool)

                    env = lmdb.Environment("TEMP_DB", readonly=True,
                        max_dbs=9+cdr3_region_len)
                    txn = env.begin()
                    subdb = env.open_db(f"{chain}_cdrs".encode(),
                        txn, create=False)
                    cursor = lmdb.Cursor(subdb, txn)
                    for i, (child_id, cdr_grp) in enumerate(
                            zip(child_ids, aligned_seqs)):
                        key = struct.pack('@I', child_id)
                        value = cursor.get(key)
                        self.assertTrue(value[:-1].decode()==cdr_grp[0])

                    subdb = env.open_db(f"{chain}_vgenes".encode(),
                        txn, create=False)
                    cursor = lmdb.Cursor(subdb, txn)
                    for i, (child_id, cdr_grp) in enumerate(
                            zip(child_ids, aligned_seqs)):
                        key = struct.pack('@I', child_id)
                        value = tuple(cursor.get(key))
                        vgene_code = get_vgene_code(vgenes[i], vspecies[i])
                        self.assertTrue(value==vgene_code)

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

                        kmer_profile = \
                            self.eval_nmbr_table_row_contents(
                                kmer_to_child, cdrs, child_ids,
                                kmer_profile, j)

                    subdb = env.open_db(f"{chain}_diversity".encode(),
                            txn, create=False)
                    cursor = lmdb.Cursor(subdb, txn)

                    cursor = lmdb.Cursor(subdb, txn)
                    for key, value in cursor.iternext():
                        int_key = struct.unpack('@i', key)[0]
                        self.assertTrue(int_key in kmer_profile)

                        self.assertTrue(len(value) % 8 == 0)
                        self.assertTrue(len(value) / 8 ==
                                kmer_profile[int_key].flatten().shape[0])

                        loaded_kmer_table = \
                                np.frombuffer(value, dtype=np.int64)
                        self.assertTrue(np.allclose(loaded_kmer_table,
                            kmer_profile[int_key].flatten()))
                        kmer_profile[int_key] = None

                    for _, value in kmer_profile.items():
                        self.assertTrue(value is None)

                shutil.rmtree("TEMP_DB")



    def test_low_quality_seqs(self):
        """Test what happens if we try to write low-quality
        sequences to the database."""
        with open("temp_data_file.fa", "w+") as fhandle:
            fhandle.write(">NAME\nAAAAAAAATTTTTTTTT\n")
            fhandle.write(">NAME\ntesting123\n")
            fhandle.write(">NAME\nEVQLEVQLEVQL\n")

        build_database_from_fasta(["temp_data_file.fa"],
            "TEMP_DB", "TEMP_FILE", numbering_scheme="imgt",
            cdr_definition_scheme="imgt",
            sequence_type="paired", receptor_type="mab",
            pid_threshold=0.7, user_memo="test",
            reject_file="REJECTS")

        # There should be no error in writing the db;
        # we expect that the sequences that failed the
        # threshold test will be written to the reject file,
        # and the relevant database tables will be empty.
        with open("REJECTS", "r", encoding="utf-8") as reject_handle:
            with open("temp_data_file.fa", "r", encoding="utf-8") as \
                    fhandle:
                for line1, line2 in zip(reject_handle, fhandle):
                    self.assertTrue(line1.strip().split('\t')[0]==line2.strip())

        shutil.rmtree("TEMP_DB")
        os.remove("temp_data_file.fa")
        os.remove("REJECTS")
        

        with open("temp_data_file.csv", "w+") as fhandle:
            fhandle.write("AAAAAAAATTTTTTTTT,H,VJ123,humanoid\n")
            fhandle.write(">NAME,H,XYZ123,alien\n")
            fhandle.write("$$ALT,H,,unknown\n")

        build_database_from_csv(["temp_data_file.csv"],"TEMP_DB",
                "TEMP_FILE",
                {"heavy_chain":0, "heavy_chaintype":1,
                    "heavy_vgene":2, "species":3},
                numbering_scheme='imgt',
                cdr_definition_scheme='imgt',
                receptor_type="mab", pid_threshold=0.7,
                user_memo="", reject_file="REJECTS",
                header_rows=0)

        # There should be no error in writing the db;
        # we expect that the sequences that failed the
        # threshold test will be written to the reject file,
        # and the relevant database tables will be empty.
        with open("REJECTS", "r", encoding="utf-8") as reject_handle:
            with open("temp_data_file.csv", "r", encoding="utf-8") as \
                    fhandle:
                for line1, line2 in zip(reject_handle, fhandle):
                    self.assertTrue(line1==line2)

        shutil.rmtree("TEMP_DB")
        os.remove("temp_data_file.csv")
        os.remove("REJECTS")


    def test_tcr_fmt_loading(self):
        """Check that tcr standard format data is loaded
        as expected."""
        current_dir = os.path.abspath(os.path.dirname(
            __file__))
        csv_filepath = os.path.join(current_dir,
                "test_data", "non_antibody_test_data",
                "tcr_simpletest.csv.gz")

        try:
            shutil.rmtree("TEMP_DB")
        except:
            pass

        expected_cdrs, expected_mainseq = [], []
        build_tcr_database_from_csv(
                [csv_filepath], "TEMP_DB", "TEMP_FILE",
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

        env = lmdb.Environment("TEMP_DB", readonly=True,
                        max_dbs=10)
        with env.begin() as txn:
            subdb = env.open_db(b"heavy_cdrs", txn, create=False)
            cursor = lmdb.Cursor(subdb, txn)
            for i, expected_cdr in enumerate(expected_cdrs):
                key = struct.pack('@I', i)
                test_seq = cursor.get(key).decode()
                self.assertTrue(expected_cdr==test_seq[:-1])

        with env.begin() as txn:
            subdb = env.open_db(b"main_seq_table", txn, create=False)
            cursor = lmdb.Cursor(subdb, txn)
            for i, expected_info in enumerate(expected_mainseq):
                key = struct.pack('@I', i)
                test_seq = cursor.get(key).decode()
                self.assertTrue(expected_info==test_seq)

        shutil.rmtree("TEMP_DB")


    def eval_nmbr_table_row_contents(self, kmer_to_child,
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

            self.assertTrue(codeval in kmer_to_child)

            self.assertTrue(bytestring in kmer_to_child[codeval])
            kmer_to_child[codeval].remove(bytestring)

        return profile_counts



def setup_canonical_numbering(numbering_scheme,
        cdr_scheme, chain_type):
    """Sets up a template aligner for the numbering
    scheme we have chosen."""
    if numbering_scheme == "imgt":
        numbering = return_imgt_canonical_numbering_cpp()

    sta = SequenceTemplateAligner(numbering,
                        chain_type, numbering_scheme,
                        cdr_scheme)
    return sta, numbering



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
    vgenes, jgenes = [], []
    vjspecies = []

    for selected_seq in selected_seqs:
        vgene, jgene, _, _, species = vj_tool.assign_vj_genes(
                selected_seq[1], selected_seq[0], "unknown",
                "identity")

        vgenes.append(vgene)
        jgenes.append(jgene)
        if species == "human":
            vjspecies.append(0)
        elif species == "mouse":
            vjspecies.append(1)
        elif species == "alpaca":
            vjspecies.append(2)
        elif species == "rabbit":
            vjspecies.append(3)
        else:
            vjspecies.append(255)

    template_aligner, canon_nmbr = setup_canonical_numbering(
            numbering_scheme, cdr_scheme, chain_type[0])
    recognized_positions = set(canon_nmbr)
    recognized_positions.add('-')
    cdr_labels = annotator.assign_cdr_labels(canon_nmbr,
            chain_type[0], cdr_scheme)

    cdr1_codes = [a for (a,c) in zip(canon_nmbr, cdr_labels)
            if c == "cdr1"]
    cdr2_codes = [a for (a,c) in zip(canon_nmbr, cdr_labels)
            if c == "cdr2"]
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


def get_vgene_code(vgene, species_code):
    """Converts the input vgene and species to a code."""
    chain_code = 255
    if vgene[2] == "A":
        chain_code = 0
    elif vgene[2] == "G":
        chain_code = 1
    elif vgene[2] == "H":
        chain_code = 2
    elif vgene[2] == "L":
        chain_code = 3
    elif vgene[2] == "K":
        chain_code = 4
    elif vgene[2] == "B":
        chain_code = 5
    elif vgene[2] == "D":
        chain_code = 6
    else:
        raise RuntimeError("Incorrect chain code found in test data.")

    extracted_nums = []
    read_now = False
    current_num = ""

    for i in range(4, len(vgene)):
        if vgene[i] == "*":
            break
        if vgene[i].isnumeric():
            if not read_now:
                read_now = True
                current_num = ""
            current_num += vgene[i]
        elif read_now:
            read_now = False
            if len(current_num) > 0:
                extracted_nums.append(int(current_num))
            current_num = ""

    if len(current_num) > 0:
        extracted_nums.append(int(current_num))
    if len(extracted_nums) < 2:
        raise RuntimeError("Invalid vgene found in test data.")

    return (chain_code, species_code, extracted_nums[0],
            extracted_nums[1])


if __name__ == "__main__":
    unittest.main()
