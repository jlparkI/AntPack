"""Test database construction for errors.
Use sqlite for testing and prototyping."""
import os
import sys
import struct
import shutil
import sqlite3
import unittest
import lmdb
import numpy as np
from antpack import (build_database_from_fasta,
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
        data_filepath = os.path.join(current_dir,
                "test_data", "addtnl_test_data.fasta.gz")

        try:
            shutil.rmtree("TEMP_DB")
        except:
            pass

        for nmbr_scheme in ['imgt']:
            for cdr_scheme in ['north', 'kabat']:
                if cdr_scheme == "north":
                    sequence_type = "single"
                    memo = ""
                elif nmbr_scheme == "kabat":
                    memo = "testing123"
                    sequence_type = "unknown"
                else:
                    memo="t"
                    sequence_type = "paired"

                build_database_from_fasta([data_filepath],
                    "TEMP_DB", numbering_scheme=nmbr_scheme,
                    cdr_definition_scheme=cdr_scheme,
                    sequence_type=sequence_type, receptor_type="mab",
                    pid_threshold=0.7, user_memo=memo,
                    reject_file=None)

                seqs, seqinfos = [], []

                for seqinfo, seq in read_fasta(data_filepath):
                    seqs.append(seq)
                    seqinfos.append(seqinfo)

                env = lmdb.Environment("TEMP_DB", readonly=True,
                        max_dbs=10)
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
                    key = struct.pack('@I', i+1)
                    self.assertTrue(seq==cursor.get(key).decode())

                subdb = env.open_db(b"seq_metadata_table",
                        txn, create=False)
                cursor = lmdb.Cursor(subdb, txn)
                for i, seqinfo in enumerate(seqinfos):
                    key = struct.pack('@I', i+1)
                    self.assertTrue(seqinfo==cursor.get(key).decode())

                sca = SingleChainAnnotator(scheme=nmbr_scheme)
                annotations = sca.analyze_seqs(seqs)
                vj_tool = VJGeneTool(scheme=nmbr_scheme)

                # Check each set of chain tables for expected info.
                for chain, chain_coding in zip(["heavy", "light"],
                        [("H",), ("L", "K")]):
                    aligned_seqs, cdrs, unusual_positions, codes, vgenes, _, vspecies, child_ids = \
                            prep_seqs_for_comparison(nmbr_scheme,
                            cdr_scheme, chain_coding, seqs,
                            annotations, sca, vj_tool)
                    subdb = env.open_db(f"{chain}_cdrs".encode(),
                        txn, create=False)
                    cursor = lmdb.Cursor(subdb, txn)
                    for i, (child_id, cdr_grp) in enumerate(
                            zip(child_ids, aligned_seqs)):
                        key = struct.pack('@I', child_id)
                        value = cursor.get(key)
                        self.assertTrue(value[:-4].decode()==cdr_grp[0])
                        for j, vcode in enumerate(
                                get_vgene_code(vgenes[i], vspecies[i])[1]):
                            self.assertTrue(int(value[-4:][j])==int(vcode))

                    subdb = env.open_db(f"{chain}_dimers".encode(),
                            txn, create=False, dupsort=True)
                    cursor = lmdb.Cursor(subdb, txn)
                    kmer_to_child = {struct.unpack('@i', key)[0]:set()
                            for key in cursor.iternext(values=False)}
                    for key in kmer_to_child.keys():
                        value = cursor.get(struct.pack('@i', key))
                        kmer_to_child[key].add(
                                struct.unpack('@I', value)[0])
                        for value in cursor.iternext_dup(keys=False):
                            kmer_to_child[key].add(
                                    struct.unpack('@I', value)[0])
                    
                    kmer_profile = \
                        self.eval_nmbr_table_row_contents(kmer_to_child,
                        cdrs, child_ids)

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
            "TEMP_DB", numbering_scheme="imgt",
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
                    self.assertTrue(line1==line2)

        shutil.rmtree("TEMP_DB")
        os.remove("temp_data_file.fa")
        os.remove("REJECTS")


    def eval_nmbr_table_row_contents(self, kmer_to_child,
            cdrs, child_ids):
        """Tests the contents of the rows from the numbering
        table."""
        profile_counts = {}

        AAMAP = {k:i for i,k in enumerate("ACDEFGHIKLMNPQRSTVWY-")}

        for i, cdr_group in enumerate(cdrs):
            # TODO: For now only considering cdr3 for indexing.
            cdr = cdr_group[2]
            cdr3len = len(cdr.replace('-', ''))
            child_id = child_ids[i]

            if cdr3len not in profile_counts:
                profile_counts[cdr3len] = \
                        np.zeros(( len(cdr) - 1, 21*21))
                profile_counts[cdr3len + 1000] = \
                        np.zeros(( len(cdr) - 2, 21*21*21))


            for j in range(0, len(cdr) - 1):
                kmer = cdr[j:j+2]
                codeval = AAMAP[kmer[0]] * 21 + AAMAP[kmer[1]]
                profile_counts[cdr3len][j, codeval] += 1
                codeval = j * 10500 + cdr3len * 250 * 10500 + \
                        codeval

                self.assertTrue(codeval in kmer_to_child)
                self.assertTrue(child_id in kmer_to_child[codeval])
                kmer_to_child[codeval].remove(child_id)


            for j in range(0, len(cdr) - 2):
                kmer = cdr[j:j+3]
                codeval = AAMAP[kmer[0]] * 21 * 21 + \
                        AAMAP[kmer[1]] * 21 + AAMAP[kmer[2]]
                profile_counts[cdr3len + 1000][j, codeval] += 1
                codeval = j * 10500 + cdr3len * 250 * 10500 + \
                        codeval + 475

                self.assertTrue(codeval in kmer_to_child)
                self.assertTrue(child_id in kmer_to_child[codeval])
                kmer_to_child[codeval].remove(child_id)

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




def prep_seqs_for_comparison(numbering_scheme,
        cdr_scheme, chain_type, sequences,
        annotations, annotator,
        vj_tool):
    """Prep sequences which are of the specified chain
    type for analysis by aligning to the template,
    extracting cdrs etc."""
    selected_seqs = [(s, a) for (s, a) in zip(
        sequences, annotations) if a[2] in chain_type]
    child_ids = [i+1 for i,a in enumerate(annotations)
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
            [cdr1_codes, cdr2_codes, cdr3_codes], \
            vgenes, jgenes, vjspecies, child_ids


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
        return -1

    family_number = 255
    for i in range(4, len(vgene)):
        if not vgene[i].isnumeric():
            break
    try:
        family_number = int(vgene[4:i])
    except:
        return -1

    if species_code == 255:
        return -1

    return family_number * 500 * 500 + chain_code * 500 + species_code,\
            (chain_code, species_code, family_number)


if __name__ == "__main__":
    unittest.main()
