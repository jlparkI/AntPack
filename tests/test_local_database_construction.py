"""Test database construction for errors.
Use sqlite for testing and prototyping."""
import os
import sqlite3
import unittest
from antpack import build_database_from_fasta, SingleChainAnnotator
from antpack.utilities import read_fasta
import numpy as np



class TestLocalDBConstruction(unittest.TestCase):



    def test_local_db_construct(self):
        """Check that local database construction yields a
        database containing the information we expect."""
        AAMAP = {k:i for i,k in enumerate("ACDEFGHIKLMNPQRSTVWY-")}

        current_dir = os.path.abspath(os.path.dirname(
            __file__))
        data_filepath = os.path.join(current_dir,
                "test_data", "addtnl_test_data.fasta.gz")

        for nmbr_scheme in ['imgt', 'kabat']:
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
                heavy_chains = [(a,s) for (a,s) in zip(annotations,
                    seqs) if a[2]=="H"]
                heavy_codes, msa = sca.build_msa(
                        [t[1] for t in heavy_chains],
                        [t[0] for t in heavy_chains])

                trimmed_numbering = [sca.trim_alignment(s,a)[1]
                        for (a,s) in heavy_chains]
                labels = sca.assign_cdr_labels(heavy_codes,
                        "H", cdr_scheme)
                rows = cursor.execute(f"SELECT * from {nmbr_scheme}"
                    "_heavy_numbering").fetchall()
                self.assertTrue([r[0] for r in rows]==
                        [s.replace('-', '') for s in msa])
                self.assertTrue([r[1].split('_') for r in rows]==
                        trimmed_numbering)

                cdr_labels = [(k,l) for (k,l) in enumerate(labels)
                        if l.startswith("cdr")]
                for i, (k,_) in enumerate(cdr_labels):
                    self.assertTrue([r[i+11] for r in rows]==
                            [m[k] for m in msa])
                del rows

                rows = cursor.execute(f"SELECT * from {nmbr_scheme}"
                    "_heavy_column_diversity").fetchall()
                for i, (k,_) in enumerate(cdr_labels):
                    gt_counts = [0]*21
                    for m in msa:
                        gt_counts[AAMAP[m[k]]] += 1
                    self.assertTrue(rows[i][0]==heavy_codes[k])
                    self.assertTrue(list(rows[i][2:])==gt_counts)

                del heavy_codes, msa, labels


                # Now do the same for light chains.
                light_chains = [(a,s) for (a,s) in zip(annotations,
                    seqs) if a[2] in ("K", "L")]
                trimmed_numbering = [sca.trim_alignment(s,a)[1]
                        for (a,s) in light_chains]
                light_codes, msa = sca.build_msa(
                        [t[1] for t in light_chains],
                        [t[0] for t in light_chains])
                labels = sca.assign_cdr_labels(light_codes,
                        "L", cdr_scheme)
                rows = cursor.execute(f"SELECT * from {nmbr_scheme}"
                    "_light_numbering").fetchall()

                self.assertTrue([r[0] for r in rows]==
                        [s.replace('-', '') for s in msa])
                self.assertTrue([r[1].split('_') for r in rows]==
                        trimmed_numbering)

                cdr_labels = [(k,l) for (k,l) in enumerate(labels)
                        if l.startswith("cdr")]
                for i, (k,_) in enumerate(cdr_labels):
                    self.assertTrue([r[i+11] for r in rows]==
                            [m[k] for m in msa])

                rows = cursor.execute(f"SELECT * from {nmbr_scheme}"
                    "_light_column_diversity").fetchall()
                for i, (k,_) in enumerate(cdr_labels):
                    gt_counts = [0]*21
                    for m in msa:
                        gt_counts[AAMAP[m[k]]] += 1
                    self.assertTrue(rows[i][0]==light_codes[k])
                    self.assertTrue(list(rows[i][2:])==gt_counts)

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

        rows = cursor.execute("SELECT * from imgt"
            "_heavy_numbering").fetchall()
        self.assertTrue(len(rows)==0)
        rows = cursor.execute("SELECT * from imgt"
            "_light_numbering").fetchall()
        self.assertTrue(len(rows)==0)

        rows = cursor.execute("SELECT * from sequences").fetchall()
        self.assertTrue(len(rows)==3)
        self.assertTrue(rows[1][0]=="testing123")
        self.assertTrue(rows[0][1]=="NAME")

        con.close()
        os.remove("TEMP_DB.db")
        os.remove("temp_data_file.fa")



if __name__ == "__main__":
    unittest.main()
