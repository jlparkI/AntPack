"""Test database construction for errors.
Use sqlite for testing and prototyping."""
import os
import sqlite3
import unittest
from antpack import build_database_from_fasta, SingleChainAnnotator
from antpack.utilities import read_fasta
from antpack.antpack_cpp_ext import (SequenceTemplateAligner,
        return_imgt_canonical_numbering_cpp)



class TestLocalDBConstruction(unittest.TestCase):



    def test_local_db_construct(self):
        """Check that local database construction yields a
        database containing the information we expect."""
        AAMAP = {k:i for i,k in enumerate("ACDEFGHIKLMNPQRSTVWY-")}

        current_dir = os.path.abspath(os.path.dirname(
            __file__))
        data_filepath = os.path.join(current_dir,
                "test_data", "addtnl_test_data.fasta.gz")
        try:
            os.remove("TEMP_DB.db")
        except:
            pass
        for nmbr_scheme in ['imgt']:
            for cdr_scheme in ['kabat', 'north', 'kabat']:
                if cdr_scheme == "north":
                    sequence_type = "single"
                    memo = ""
                elif nmbr_scheme == "kabat":
                    memo = "testing123"
                    sequence_type = "unknown"
                else:
                    memo="t"
                    sequence_type = "paired"

                build_database_from_fasta(data_filepath,
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
                del rows

                # Next, check that the sequences and seqinfo in
                # the db match what we loaded.
                rows = cursor.execute("SELECT * from sequences").fetchall()
                self.assertTrue(seqs==[r[0] for r in rows])
                self.assertTrue(seqinfos==[r[1] for r in rows])
                del rows

                sca = SingleChainAnnotator(scheme=nmbr_scheme)
                annotations = sca.analyze_seqs(seqs)

                # Check heavy chains first. Make sure numbering table
                # contains expected info and that cdr only columns
                # have been extracted.
                aligned_seqs, cdr1, cdr2, cdr3, unusual_positions, codes = \
                        prep_seqs_for_comparison(nmbr_scheme,
                        cdr_scheme, ("H",), seqs,
                        annotations, sca)
                rows = cursor.execute("SELECT * from heavy_numbering").fetchall()
                self.assertTrue([r[0] for r in rows] == aligned_seqs)
                self.assertTrue([r[1] for r in rows] == cdr1)
                self.assertTrue([r[2] for r in rows] == cdr2)
                self.assertTrue([r[3] for r in rows] == cdr3)
                del rows

                rows = cursor.execute("SELECT * from "
                    "heavy_column_diversity").fetchall()
                cdr_labels = sca.assign_cdr_labels(codes,
                        "H", cdr_scheme)
                k = 0
                for i, label in enumerate(cdr_labels):
                    if not label.startswith("cdr"):
                        continue
                    gt_counts = [0]*21
                    for m in aligned_seqs:
                        gt_counts[AAMAP[m[i]]] += 1
                    self.assertTrue(rows[k][0]==codes[i])
                    self.assertTrue(list(rows[k][2:])==gt_counts)
                    k += 1

                del codes, aligned_seqs, cdr1, cdr2, cdr3, unusual_positions


                # Now do the same for light chains.
                aligned_seqs, cdr1, cdr2, cdr3, unusual_positions, codes = \
                        prep_seqs_for_comparison(nmbr_scheme,
                        cdr_scheme, ("L", "K"), seqs,
                        annotations, sca)
                rows = cursor.execute("SELECT * from light_numbering").fetchall()

                self.assertTrue([r[0] for r in rows] == aligned_seqs)
                self.assertTrue([r[1] for r in rows] == cdr1)
                self.assertTrue([r[2] for r in rows] == cdr2)
                self.assertTrue([r[3] for r in rows] == cdr3)

                rows = cursor.execute("SELECT * from "
                    "light_column_diversity").fetchall()
                cdr_labels = sca.assign_cdr_labels(codes,
                        "L", cdr_scheme)
                k = 0
                for i, label in enumerate(cdr_labels):
                    if not label.startswith("cdr"):
                        continue
                    gt_counts = [0]*21
                    for m in aligned_seqs:
                        gt_counts[AAMAP[m[i]]] += 1
                    self.assertTrue(rows[k][0]==codes[i])
                    self.assertTrue(list(rows[k][2:])==gt_counts)
                    k += 1

                con.close()
                os.remove("TEMP_DB.db")


    def test_low_quality_seqs(self):
        """Test what happens if we try to write low-quality
        sequences to the database."""
        with open("temp_data_file.fa", "w+") as fhandle:
            fhandle.write(">NAME\nAAAAAAAATTTTTTTTT\n")
            fhandle.write(">NAME\ntesting123\n")
            fhandle.write(">NAME\nEVQLEVQLEVQL\n")

        build_database_from_fasta("temp_data_file.fa",
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
        self.assertTrue(rows[1][0]=="testing123")
        self.assertTrue(rows[0][1]=="NAME")

        con.close()
        os.remove("TEMP_DB.db")
        os.remove("temp_data_file.fa")





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
        annotations, annotator):
    """Prep sequences which are of the specified chain
    type for analysis by aligning to the template,
    extracting cdrs etc."""
    selected_seqs = [(s, a) for (s, a) in zip(
        sequences, annotations) if a[2] in chain_type]
    template_aligner, canon_nmbr = setup_canonical_numbering(
            numbering_scheme, cdr_scheme, chain_type[0])
    recognized_positions = set(canon_nmbr)
    cdr_labels = annotator.assign_cdr_labels(canon_nmbr,
            chain_type[0], cdr_scheme)

    aligned_seqs, cdr1, cdr2, cdr3, unusual_positions = \
            [], [], [], [], []

    for selected_seq in selected_seqs:
        aligned_seq = template_aligner.align_sequence(
                selected_seq[0], selected_seq[1][0], False)
        unusual_positions = [a for a in selected_seq[1][0]
                if a not in recognized_positions]
        cdr1.append(''.join([a for (a,c) in zip(aligned_seq,
            cdr_labels) if c == "cdr1"]))
        cdr2.append(''.join([a for (a,c) in zip(aligned_seq,
            cdr_labels) if c == "cdr2"]))
        cdr3.append(''.join([a for (a,c) in zip(aligned_seq,
            cdr_labels) if c == "cdr3"]))
        aligned_seqs.append(aligned_seq)

    return aligned_seqs, cdr1, cdr2, cdr3, unusual_positions, \
            canon_nmbr


if __name__ == "__main__":
    unittest.main()
