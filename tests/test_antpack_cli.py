"""Tests basic functionality for the antpack cli. Note that this
relies on the functionality of the SingleChainAnnotator,
PairedChainAnnotator and VJGeneTool classes, so if those tests
do not pass, this one should not either."""
import os
import subprocess
import random
import gzip
import unittest
from antpack import PairedChainAnnotator, SingleChainAnnotator, VJGeneTool



class TestAntPackCLI(unittest.TestCase):




    def test_paired_chain_functionality(self):
        """Take test heavy and light chains and pair them,
        then use the CLI to do basic analysis using a couple
        different options (csv, no csv, assign vj, do not assign
        vj etc). Check to make sure the output file contains
        everything we expect."""
        project_path = os.path.abspath(os.path.dirname(__file__))
        current_dir = os.getcwd()
        os.chdir(os.path.join(project_path, "test_data"))
        with gzip.open("test_data.csv.gz", "rt") as fhandle:
            _ = fhandle.readline()
            seqs = [line.strip().split(",")[0] for line in fhandle]
        os.chdir(current_dir)

        sc_aligner = SingleChainAnnotator(scheme="imgt")
        pc_aligner = PairedChainAnnotator(scheme="imgt")

        alignments = [sc_aligner.analyze_seq(seq) for seq in seqs]

        heavy_chains = [(seq, a) for a, seq in zip(alignments, seqs)
                if a[2] == "H"]
        light_chains = [(seq, a) for a, seq in zip(alignments, seqs)
                if a[2] != "H"]

        random.seed(0)

        # Join the chains using a very arbitrary linker.
        merged_seqs = ["".join( [light_chain[0], "YYYYGGGGSSSS", heavy_chain[0]] )
                for (heavy_chain, light_chain) in zip(heavy_chains, light_chains)]

        heavy_pc_annotations, light_pc_annotations = [], []

        for merged_seq in merged_seqs:
            h_annot, l_annot = pc_aligner.analyze_seq(merged_seq)
            heavy_pc_annotations.append(h_annot)
            light_pc_annotations.append(l_annot)

        # Create a temporary fasta input file.
        with open("input_fasta_data.fasta", "w+", encoding="utf-8") as fhandle:
            for i, merged_seq in enumerate(merged_seqs):
                fhandle.write(f">this is sequence # {i}\n{merged_seq}\n")


        _ = subprocess.run(["AntPack", "input_fasta_data.fasta", "output_data",
            "imgt", "--paired"])

        self.assertTrue("output_data_heavy.fasta" in os.listdir())
        self.assertTrue("output_data_light.fasta" in os.listdir())

        os.remove("input_fasta_data.fasta")
        os.remove("output_data_heavy.fasta")
        os.remove("output_data_light.fasta")


    def test_single_chain_functionality(self):
        """Take test heavy and light chains,
        then use the CLI to do basic analysis using a couple
        different options (csv, no csv, assign vj, do not assign
        vj etc). Check to make sure the output file contains
        everything we expect."""
        project_path = os.path.abspath(os.path.dirname(__file__))
        current_dir = os.getcwd()
        os.chdir(os.path.join(project_path, "test_data"))
        with gzip.open("test_data.csv.gz", "rt") as fhandle:
            _ = fhandle.readline()
            seqs = [line.strip().split(",")[0] for line in fhandle]
        os.chdir(current_dir)


if __name__ == "__main__":
    unittest.main()
