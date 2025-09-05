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

                os.remove("TEMP_DB.db")
                os.remove("TEMP_DB.db-shm")
                os.remove("TEMP_DB.db-wal")




if __name__ == "__main__":
    unittest.main()
