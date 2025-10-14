"""Test database construction for errors.
Use sqlite for testing and prototyping."""
import os
import math
import sqlite3
import unittest
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
                    "TEMP_DB.db", numbering_scheme=nmbr_scheme,
                    cdr_definition_scheme=cdr_scheme,
                    sequence_type=sequence_type, receptor_type="mab",
                    pid_threshold=0.7, user_memo=memo)

                seqs, seqinfos = [], []

                for seqinfo, seq in read_fasta(data_filepath):
                    seqs.append(seq)
                    seqinfos.append(seqinfo)

                con = sqlite3.connect("TEMP_DB.db")
                cursor = con.cursor()

                # First, check that the database metadata is what
                # was expected.
                rows = cursor.execute("SELECT * from database_info").fetchall()
                self.assertTrue(rows[0][0]==nmbr_scheme)
                self.assertTrue(rows[0][1]=="mab")
                self.assertTrue(rows[0][2]==sequence_type)
                self.assertTrue(rows[0][3]==cdr_scheme)
                self.assertTrue(rows[0][4]==memo)
                self.assertTrue(rows[0][5]==0.7)
                del rows

                # Next, check that the sequences and seqinfo in
                # the db match what we loaded.
                rows = cursor.execute("SELECT * from sequences").fetchall()
                self.assertTrue(seqs==[r[0] for r in rows])
                self.assertTrue(seqinfos==[r[1] for r in rows])
                del rows

                sca = SingleChainAnnotator(scheme=nmbr_scheme)
                annotations = sca.analyze_seqs(seqs)
                vj_tool = VJGeneTool(scheme=nmbr_scheme)

                # Check heavy chains first. Make sure numbering table
                # contains expected info.
                aligned_seqs, cdrs, unusual_positions, codes, vgenes, _, vspecies = \
                        prep_seqs_for_comparison(nmbr_scheme,
                        cdr_scheme, ("H",), seqs,
                        annotations, sca, vj_tool)
                rows = cursor.execute("SELECT * from heavy_numbering").fetchall()
                for i, row in enumerate(rows):
                    self.assertTrue(row[0]==aligned_seqs[i][0])
                    self.assertTrue(row[1]==aligned_seqs[i][1])
                    vcode = get_vgene_code(vgenes[i], vspecies[i], vj_tool)
                    self.assertTrue(row[3]==vcode)

                dimer_profile, trimer_profile = \
                        self.eval_nmbr_table_row_contents(rows,
                        cdrs, unusual_positions)

                del rows

                rows = cursor.execute("SELECT * from "
                    "heavy_column_diversity").fetchall()

                for row in rows:
                    test_arr = np.frombuffer(row[-1], dtype=np.int64)
                    cdr_idx = int(row[0][-1]) - 1
                    if row[1] == "dimer":
                        profile_table = dimer_profile[cdr_idx]
                    else:
                        profile_table = trimer_profile[cdr_idx]

                    self.assertTrue(row[2]==profile_table.shape[0])

                    test_arr = test_arr.reshape((row[2], row[3]))
                    self.assertTrue(np.allclose(test_arr, profile_table))

                del codes, aligned_seqs, cdrs, unusual_positions,\
                        dimer_profile, trimer_profile


                # Now do the same for light chains.
                aligned_seqs, cdrs, unusual_positions, codes, vgenes, _, vspecies = \
                        prep_seqs_for_comparison(nmbr_scheme,
                        cdr_scheme, ("L","K"), seqs,
                        annotations, sca, vj_tool)
                rows = cursor.execute("SELECT * from light_numbering").fetchall()
                for i, row in enumerate(rows):
                    self.assertTrue(row[0]==aligned_seqs[i][0])
                    self.assertTrue(row[1]==aligned_seqs[i][1])
                    vcode = get_vgene_code(vgenes[i], vspecies[i], vj_tool)
                    self.assertTrue(row[3]==vcode)

                dimer_profile, trimer_profile = \
                        self.eval_nmbr_table_row_contents(rows,
                        cdrs, unusual_positions)

                del rows

                rows = cursor.execute("SELECT * from "
                    "light_column_diversity").fetchall()

                for row in rows:
                    test_arr = np.frombuffer(row[-1], dtype=np.int64)
                    cdr_idx = int(row[0][-1]) - 1
                    if row[1] == "dimer":
                        profile_table = dimer_profile[cdr_idx]
                    else:
                        profile_table = trimer_profile[cdr_idx]

                    self.assertTrue(row[2]==profile_table.shape[0])

                    test_arr = test_arr.reshape((row[2], row[3]))
                    self.assertTrue(np.allclose(test_arr, profile_table))

                del codes, aligned_seqs, cdrs, unusual_positions,\
                        dimer_profile, trimer_profile

                con.close()
                os.remove("TEMP_DB.db")


    def test_low_quality_seqs(self):
        """Test what happens if we try to write low-quality
        sequences to the database."""
        with open("temp_data_file.fa", "w+") as fhandle:
            fhandle.write(">NAME\nAAAAAAAATTTTTTTTT\n")
            fhandle.write(">NAME\ntesting123\n")
            fhandle.write(">NAME\nEVQLEVQLEVQL\n")

        build_database_from_fasta(["temp_data_file.fa"],
            "TEMP_DB.db", numbering_scheme="imgt",
            cdr_definition_scheme="imgt",
            sequence_type="paired", receptor_type="mab",
            pid_threshold=0.7, user_memo="test")

        # There should be no error in writing the db;
        # we expect that the sequences that failed the
        # threshold test will still be inserted into the
        # raw sequences table but not the heavy numbering
        # or light numbering tables.
        con = sqlite3.connect("TEMP_DB.db")
        cursor = con.cursor()
        rows = cursor.execute("SELECT * from database_info").fetchall()
        self.assertTrue(rows[0][0]=="imgt")
        self.assertTrue(rows[0][1]=="mab")
        self.assertTrue(rows[0][2]=="paired")
        self.assertTrue(rows[0][3]=="imgt")
        self.assertTrue(rows[0][4]=="test")
        del rows

        rows = cursor.execute("SELECT * from heavy_numbering").fetchall()
        self.assertTrue(len(rows)==0)
        rows = cursor.execute("SELECT * from light_numbering").fetchall()
        self.assertTrue(len(rows)==0)

        rows = cursor.execute("SELECT * from sequences").fetchall()
        self.assertTrue(len(rows)==3)

        con.close()
        os.remove("TEMP_DB.db")
        os.remove("temp_data_file.fa")

    def eval_nmbr_table_row_contents(self, rows, cdrs,
            unusual_positions):
        """Tests the contents of the rows from the numbering
        table."""
        dimer_profile_counts = [np.zeros(( int((len(cdr)+1)/2), 21*21 ))
                for cdr in cdrs[0]]
        trimer_profile_counts = [np.zeros(( int((len(cdr)+2)/3), 21*21*21 ))
                for cdr in cdrs[0]]

        AAMAP = {k:i for i,k in enumerate("ACDEFGHIKLMNPQRSTVWY-")}

        for i, cdr_group in enumerate(cdrs):
            self.assertTrue(rows[i][2] == unusual_positions[i])
            cdr3len = len(cdr_group[2].replace('-', ''))
            self.assertTrue(rows[i][6] == cdr3len)
            # Counts where we are in order to skip things like vj
            # genes not checked in this test.
            row_counter = 7

            # TODO: For now only considering cdr3 for indexing.
            for cdr, dimer_count, trimer_count in \
                    zip(cdr_group[2:], dimer_profile_counts[2:],
                            trimer_profile_counts[2:]):

                for j in range(0, len(cdr), 2):
                    dimer = cdr[j:j+2]
                    if len(dimer) == 1:
                        dimer += '-'

                    codeval = AAMAP[dimer[0]] * 21 + AAMAP[dimer[1]]
                    dimer_count[int(j/2), codeval] += 1
                    self.assertTrue(rows[i][row_counter] == codeval)
                    row_counter += 1

                for j in range(0, len(cdr), 3):
                    trimer = cdr[j:j+3]
                    while len(trimer) < 3:
                        trimer += '-'

                    codeval = AAMAP[trimer[0]] * 21 * 21 + \
                        AAMAP[trimer[1]] * 21 + AAMAP[trimer[2]]
                    trimer_count[int(j/3), codeval] += 1
                    self.assertTrue(rows[i][row_counter] == codeval)
                    row_counter += 1

        return dimer_profile_counts, trimer_profile_counts



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
            vgenes, jgenes, vjspecies


def get_vgene_code(vgene, species_code, vj_tool):
    """Converts the input vgene and species to a code."""
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

    family_number = -1
    for i in range(4, len(vgene)):
        if not vgene[i].isnumeric():
            break
    try:
        family_number = int(vgene[4:i])
    except:
        return -1

    if species_code == -1:
        return -1

    return family_number * 500 * 500 + chain_code * 500 + species_code


if __name__ == "__main__":
    unittest.main()
