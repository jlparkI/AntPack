"""Test database construction for errors.
Use sqlite for testing and prototyping."""
import os
import math
import sqlite3
import unittest
import numpy as np
from antpack import (build_database_from_fasta,
        SingleChainAnnotator, EMCategoricalMixture)
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
                    "TEMP_DB.db", os.getcwd(),
                    numbering_scheme=nmbr_scheme,
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
                # contains expected info.
                aligned_seqs, cdrs, unusual_positions, codes, cluster_profiles = \
                        prep_seqs_for_comparison(nmbr_scheme,
                        cdr_scheme, ("H",), seqs,
                        annotations, sca)
                rows = cursor.execute("SELECT * from heavy_numbering").fetchall()
                self.assertTrue([r[0] for r in rows] == aligned_seqs)

                for i, cdr_group in enumerate(cdrs):
                    self.assertTrue(rows[i][1] == unusual_positions[i])
                    row_counter = 8
                    for cdr in cdr_group:
                        for aa_code in cdr:
                            self.assertTrue(rows[i][row_counter] ==
                                    aa_code)
                            row_counter += 1
                del rows

                rows = cursor.execute("SELECT * from "
                    "heavy_column_diversity").fetchall()
                ctr = 0
                for cluster_profile in cluster_profiles:
                    for i in range(cluster_profile.shape[0]):
                        for j in range(cluster_profile.shape[1]):
                            self.assertTrue(np.allclose(rows[ctr][2:],
                                cluster_profile[i,j,:]))
                            ctr += 1

                del codes, aligned_seqs, cdrs, unusual_positions


                # Now do the same for light chains.
                aligned_seqs, cdrs, unusual_positions, codes, cluster_profiles = \
                        prep_seqs_for_comparison(nmbr_scheme,
                        cdr_scheme, ("L", "K"), seqs,
                        annotations, sca)
                rows = cursor.execute("SELECT * from light_numbering").fetchall()
                self.assertTrue([r[0] for r in rows] == aligned_seqs)

                for i, cdr_group in enumerate(cdrs):
                    self.assertTrue(rows[i][1] == unusual_positions[i])
                    row_counter = 8
                    for cdr in cdr_group:
                        for aa_code in cdr:
                            self.assertTrue(rows[i][row_counter] ==
                                    aa_code)
                            row_counter += 1
                del rows

                rows = cursor.execute("SELECT * from "
                    "light_column_diversity").fetchall()
                ctr = 0
                for cluster_profile in cluster_profiles:
                    for i in range(cluster_profile.shape[0]):
                        for j in range(cluster_profile.shape[1]):
                            self.assertTrue(np.allclose(rows[ctr][2:],
                                cluster_profile[i,j,:]))
                            ctr += 1

                del codes, aligned_seqs, cdrs, unusual_positions, cluster_profiles

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
            "TEMP_DB.db", os.getcwd(),
            numbering_scheme="imgt",
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
        self.assertTrue(len(rows)==0)

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
    recognized_positions.add('-')
    cdr_labels = annotator.assign_cdr_labels(canon_nmbr,
            chain_type[0], cdr_scheme)

    aligned_seqs, cdrs, unusual_positions = \
            [], [], []

    for i, selected_seq in enumerate(selected_seqs):
        aligned_seq = template_aligner.align_sequence(
                selected_seq[0], selected_seq[1][0], False)
        unusual_positions.append( len([a for a in selected_seq[1][0]
                if a not in recognized_positions]) )
        aligned_seqs.append(aligned_seq)

    cdr_codes, cluster_profiles, cluster_assignments = [], [], []

    for region in ["cdr1", "cdr2", "cdr3"]:
        cdr_codes.append([a for (a,c) in zip(canon_nmbr, cdr_labels)
            if c == region])
        nclusters = int(10 * max(1,
            math.log10(len(selected_seqs) / 100)))
        em_cat = EMCategoricalMixture(n_components=nclusters,
                numbering=canon_nmbr, chain_type=chain_type[0],
                numbering_scheme=numbering_scheme,
                cdr_scheme=cdr_scheme,
                region=region, max_threads=2,
                verbose=False)
        em_cat.fit(aligned_seqs, max_iter=25,
                tol=1e-2, n_restarts=3,
                random_state=123,
                prune_after_fitting=True)
        cluster_profile = em_cat.initialize_cluster_profiles()
        preds = em_cat.update_cluster_profiles(aligned_seqs,
                cluster_profile)
        cluster_profiles.append(cluster_profile)
        cluster_assignments.append(preds)

    for i, aligned_seq in enumerate(aligned_seqs):
        grouped_cdr = []
        for j, region in enumerate(["cdr1", "cdr2", "cdr3"]):
            assn_code = str(cluster_assignments[j][i])
            cdr = [assn_code + a for (a,c) in zip(aligned_seq,
                cdr_labels) if c == region]
            grouped_cdr.append(cdr)

        cdrs.append(grouped_cdr)


    return aligned_seqs, cdrs, unusual_positions, \
            cdr_codes, cluster_profiles


if __name__ == "__main__":
    unittest.main()
