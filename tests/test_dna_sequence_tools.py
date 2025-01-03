"""Tests basic functionality for the DNASeqTranslator
class. We do this essentially just by comparing to
Biopython."""
import os
import gzip
import random
import unittest
from Bio.Seq import Seq
from antpack import DNASeqTranslator


class TestDNASeqTranslator(unittest.TestCase):


    def test_error_checking(self):
        """Check that known incorrect inputs raise
        appropriate errors."""
        dna_translator = DNASeqTranslator()

        with self.assertRaises(RuntimeError):
            dna_translator.translate_dna_known_rf("AaTTTC")
        with self.assertRaises(RuntimeError):
            dna_translator.translate_dna_known_rf("ALTTTC")
        with self.assertRaises(RuntimeError):
            dna_translator.translate_dna_known_rf("ATTTC", -1)
        with self.assertRaises(RuntimeError):
            dna_translator.translate_dna_known_rf("ATTTC", 4)

        with self.assertRaises(RuntimeError):
            dna_translator.translate_dna_unknown_rf("AaTTTC")
        with self.assertRaises(RuntimeError):
            dna_translator.translate_dna_unknown_rf("ALTTTC")



    def test_known_rf(self):
        """Check that the known_rf function returns correct /
        expected input."""
        dna_translator = DNASeqTranslator()
        random.seed(123)

        for i in range(1000):
            new_seq = []

            for j in range(random.randint(3,300)):
                new_seq.append(random.choice(('A', 'G', 'T', 'C')))

            new_seq = "".join(new_seq)
            bioseq = str(Seq(new_seq).translate())
            antseq = dna_translator.translate_dna_known_rf(new_seq, 0)
            self.assertTrue(bioseq==antseq)

            bioseq = str(Seq(new_seq[1:]).translate())
            antseq = dna_translator.translate_dna_known_rf(new_seq, 1)
            self.assertTrue(bioseq==antseq)

            bioseq = str(Seq(new_seq[2:]).translate())
            antseq = dna_translator.translate_dna_known_rf(new_seq, 2)
            self.assertTrue(bioseq==antseq)

            bioseq = str(Seq(new_seq).reverse_complement().translate())
            antseq = dna_translator.translate_dna_known_rf(new_seq, 0,
                    True)
            self.assertTrue(bioseq==antseq)

            bioseq = str(Seq(new_seq).reverse_complement()[1:].translate())
            antseq = dna_translator.translate_dna_known_rf(new_seq, 1,
                    True)
            self.assertTrue(bioseq==antseq)

            bioseq = str(Seq(new_seq).reverse_complement()[2:].translate())
            antseq = dna_translator.translate_dna_known_rf(new_seq, 2,
                    True)
            self.assertTrue(bioseq==antseq)


    def test_unknown_rf(self):
        """Check that the unknown_rf function can
        correctly find the heavy / light chain sequence
        even when reading frame & orientation are not
        known."""
        dna_translator = DNASeqTranslator()
        random.seed(123)

        project_path = os.path.abspath(os.path.dirname(__file__))
        current_dir = os.getcwd()
        os.chdir(os.path.join(project_path, "test_data"))

        with gzip.open("dna_translate_test_data.csv.gz", "rt") as fhandle:
            _ = fhandle.readline()
            dna_heavy_seqs, dna_light_seqs = [], []
            aa_heavy_seqs, aa_light_seqs = [], []
            for line in fhandle:
                elements = line.strip().split(",")
                # Add some arbitrary additional noise onto the ends.
                # Randomly choose to use the reverse complement
                # instead of the actual sequence.
                dna_heavy = "AAGT" + elements[0] + "TCGA"
                if random.randint(0,1) == 1:
                    dna_heavy = str(Seq(dna_heavy).reverse_complement())
                dna_heavy_seqs.append(dna_heavy)
                aa_heavy_seqs.append(elements[1])

                dna_light = "AAGC" + elements[2] + "TCTA"
                if random.randint(0,1) == 1:
                    dna_light = str(Seq(dna_light).reverse_complement())
                dna_light_seqs.append(dna_light)
                aa_light_seqs.append(elements[3])

        for hdna, haa, ldna, laa in zip(
                dna_heavy_seqs, aa_heavy_seqs,
                dna_light_seqs, aa_light_seqs):
            anth = dna_translator.translate_dna_unknown_rf(hdna, True)
            antl = dna_translator.translate_dna_unknown_rf(ldna, True)
            self.assertTrue(haa in anth)
            self.assertTrue(laa in antl)

        os.chdir(current_dir)





if __name__ == "__main__":
    unittest.main()
