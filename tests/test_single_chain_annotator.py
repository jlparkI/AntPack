"""Tests basic functionality for the SingleChainAnnotator class."""
import os
import re
import random
import gzip
import copy
import unittest
from antpack import SingleChainAnnotator

class TestSingleChainAnnotator(unittest.TestCase):


    def test_error_checking(self):
        """Check that sequences which have known issues are flagged
        as such, and that deliberately invalid inputs are recognized."""
        # Pass dummy sequences with errors.
        aligner = SingleChainAnnotator(chains=["H", "K", "L"])
        results = aligner.analyze_seqs(["YaY"])
        self.assertTrue(results[0][3].startswith("Sequence contains invalid"))
        results = aligner.analyze_seqs(["YBW"])
        self.assertTrue(results[0][3].startswith("Sequence contains invalid"))
        results = aligner.analyze_seqs(["Y K"])
        self.assertTrue(results[0][3].startswith("Sequence contains invalid"))
        results = aligner.analyze_seqs(["Y-K"])
        self.assertTrue(results[0][3].startswith("Sequence contains invalid"))

        results = aligner.analyze_seq("Y-K")
        self.assertTrue(results[3].startswith("Sequence contains invalid"))
        results = aligner.analyze_seq("yAy")
        self.assertTrue(results[3].startswith("Sequence contains invalid"))

        # Repeat these tests with a SingleChainAnnotator for TCRs.
        aligner = SingleChainAnnotator(chains=["A", "B", "D", "G"])
        results = aligner.analyze_seqs(["YaY"])
        self.assertTrue(results[0][3].startswith("Sequence contains invalid"))
        results = aligner.analyze_seqs(["YBW"])
        self.assertTrue(results[0][3].startswith("Sequence contains invalid"))
        results = aligner.analyze_seqs(["Y K"])
        self.assertTrue(results[0][3].startswith("Sequence contains invalid"))
        results = aligner.analyze_seqs(["Y-K"])
        self.assertTrue(results[0][3].startswith("Sequence contains invalid"))

        results = aligner.analyze_seq("Y-K")
        self.assertTrue(results[3].startswith("Sequence contains invalid"))
        results = aligner.analyze_seq("yAy")
        self.assertTrue(results[3].startswith("Sequence contains invalid"))

        # It should be impossible to create an annotator for both TCRs and
        # mAbs -- let's check.
        with self.assertRaises(RuntimeError):
            aligner = SingleChainAnnotator(chains=["A", "B", "D", "H"])
        with self.assertRaises(RuntimeError):
            aligner = SingleChainAnnotator(chains=["H", "K", "A"])

        # The only scheme allowed for TCRs should be IMGT.
        with self.assertRaises(RuntimeError):
            aligner = SingleChainAnnotator(chains=["A", "B", "D", "G"], scheme="martin")



    def test_chain_recognition(self):
        """Ensure that single chain annotator can correctly recognize the
        input chain when supplied with something that could be L or H,
        and ensure it can correctly detect sequences with large deletions
        that remove one or more conserved residues."""
        known_k = ("DIVMTQSPSSLTVTAGEKVTMSCKSSQSLLSSGNQKNYLTWYQQIPGQPPKLLIYWASTR"
                    "ESGVPDRFTGSGSGTDFTLTINSVQAEDLAVYYCQNDYTYPLTFGAGTKLELKRTV")
        known_l = ("QSALTQPASVSGSPGQSITISCTGTTSDVGTYNFVSWYQQHPGKAPKAIIFDVTNRPSGI"
                    "SNRFSGSKFGNTASLTISGLQAEDEADYYCAAYTVASTLLFGGGTKVTVLRQP")
        known_h = ("VHLQQSGAELMKPGASVKISCKASGYTFITYWIEWVKQRPGHGLEWIGDILPGSGSTNYN"
                    "ENFKGKATFTADSSSNTAYMQLSSLTSEDSAVYYCARSGYYGNSGFAYWGQGTLVTVSA")

        for scheme in ["martin", "imgt", "kabat", "aho"]:
            aligner = SingleChainAnnotator(chains=["H", "K", "L"],
                        scheme = scheme)
            results = aligner.analyze_seqs([known_k, known_l, known_h])
            self.assertTrue(results[0][2] == "K")
            self.assertTrue(results[1][2] == "L")
            self.assertTrue(results[2][2] == "H")

            self.assertTrue(aligner.analyze_seq(known_k)[2] == "K")
            self.assertTrue(aligner.analyze_seq(known_l)[2] == "L")
            self.assertTrue(aligner.analyze_seq(known_h)[2] == "H")

            bad_chain = known_h[:100]
            self.assertTrue(aligner.analyze_seqs([bad_chain])[0][3].startswith("Unexpected"))
            self.assertTrue(aligner.analyze_seq(bad_chain)[3].startswith("Unexpected"))

        # Repeat this check, but this time using TCRs; the only allowed scheme for these
        # is IMGT.
        known_a = ("MQQVRQSPQSLTVWEGETAILNCSYENSAFDYFPWYQQFPGEGPALLIAIRSVSDKK"
                "EDGRFTIFFNKREKKLSLHITDSQPGDSATYFCAASATGANTGKLTFGHGTILRVHP")
        known_b = ("DAGVIQSPRHEVTEMGQEVTLRCKPISGHNSLFWYRQTMMRGLELLIYFNNNVPIDD"
                "SGMPEDRFSAKMPNASFSTLKIQPSEPRDSAVYFCASTWGRASTDTQYFGPGTRLTVL")
        known_d = ("AQKVTQAQSSVSMPVRKAVTLNCLYETSWWSYYIFWYKQLPSKEMIFLIRQGSDE"
                "QNAKSGRYSVNFKKAAKSVALTISALQLEDSAKYFCALGDPGGLNTDKLIFGKGTRVTVEP")
        known_g = ("SSNLEGGTKSVTRPTRSSAEITCDLTVINAFYIHWYLHQEGKAPQRLLYYDVSNSKDVLE"
                "SGLSPGKYYTHTPRRWSWILILRNLIENDSGVYYCATWDRGNPKTHYYKKLFGSGTTLVVTD")
        aligner = SingleChainAnnotator(chains=["A", "B", "D", "G"])
        results = aligner.analyze_seqs([known_a, known_b, known_d, known_g])
        self.assertTrue(results[0][2] == "A")
        self.assertTrue(results[1][2] == "B")
        self.assertTrue(results[2][2] == "D")
        self.assertTrue(results[3][2] == "G")



    def test_mab_performance(self):
        """Run a batch of test data (approximately 1600 sequences from the
        PDB) to ensure that numbering is consistent with numbering generated
        by another tool. There will occasionally be small differences in
        cases where there are multiple possible acceptable alignments,
        but in general we expect the numbering to be the same the vast
        majority of the time. Also check that the sequences can be correctly
        formed into an MSA. This test is for mAbs specifically, tcrs are
        tested separately since they use a different alignment workflow."""
        project_path = os.path.abspath(os.path.dirname(__file__))
        current_dir = os.getcwd()
        os.chdir(os.path.join(project_path, "test_data"))

        with gzip.open("test_data.csv.gz", "rt") as fhandle:
            _ = fhandle.readline()
            seqs, martin_num, imgt_num, kabat_num, aho_num = \
                    [], [], [], [], []
            for line in fhandle:
                line_elements = line.strip().split(",")
                seqs.append(line_elements[0])
                martin_num.append(line_elements[1].split("_"))
                imgt_num.append(line_elements[2].split("_"))
                kabat_num.append(line_elements[3].split("_"))
                aho_num.append(line_elements[4].split("_"))

        os.chdir(current_dir)

        numberings = [martin_num, kabat_num, imgt_num, aho_num]
        schemes = ["martin", "kabat", "imgt", "aho"]

        for scheme, numbering in zip(schemes, numberings):
            aligner = SingleChainAnnotator(chains=["H", "K", "L"],
                scheme=scheme)
            aligner_results = aligner.analyze_seqs(seqs)
            total_comparisons, num_correct = compare_results(aligner_results,
                                numbering, seqs, scheme)
            print(f"{scheme}: Total comparisons: {total_comparisons}. "
                    f"Num matching: {num_correct}.")
            self.assertTrue(num_correct / total_comparisons > 0.97)

        # The last one produced is IMGT. Use this to test MSA construction.
        hnumbering, lnumbering, hseqs, lseqs = [], [], [], []

        for seq, numbering in zip(seqs, aligner_results):
            if numbering[2] == "H":
                hnumbering.append(numbering)
                hseqs.append(seq)
            elif numbering[2] in ["K", "L"]:
                lnumbering.append(numbering)
                lseqs.append(seq)

        observed_positions = set()
        for n in hnumbering:
            for pos in n[0]:
                observed_positions.add(pos)

        hpositions, hmsa = aligner.build_msa(hseqs, hnumbering, False)
        lpositions, lmsa = aligner.build_msa(lseqs, lnumbering, False)

        for position_set, msa in [(hpositions, hmsa), (lpositions, lmsa)]:
            for msa_seq in msa:
                self.assertTrue(len(msa_seq) == len(position_set))


    def test_tcr_performance(self):
        """Run a batch of test data (approximately 525 sequences from the
        PDB) to ensure that numbering is consistent with numbering generated
        by another tool. There will occasionally be small differences in
        cases where there are multiple possible acceptable alignments,
        but in general we expect the numbering to be the same the vast
        majority of the time. Also check that the sequences can be correctly
        formed into an MSA. This test is for TCRs specifically, mabs are
        tested separately since they use a different alignment workflow."""
        project_path = os.path.abspath(os.path.dirname(__file__))
        current_dir = os.getcwd()
        os.chdir(os.path.join(project_path, "test_data"))

        with gzip.open("tcr_test_data.csv.gz", "rt") as fhandle:
            _ = fhandle.readline()
            seqs, imgt_num = [], []
            for line in fhandle:
                line_elements = line.strip().split(",")
                seqs.append(line_elements[0])
                imgt_num.append(line_elements[1].split("_"))

        os.chdir(current_dir)

        aligner = SingleChainAnnotator(chains=["A", "B", "D", "G"])
        aligner_results = aligner.analyze_seqs(seqs)
        total_comparisons, num_correct = compare_results(aligner_results,
                                    imgt_num, seqs, "imgt")
        print(f"*TCRS, IMGT: Total comparisons: {total_comparisons}. "
                f"Num matching: {num_correct}.")
        self.assertTrue(num_correct / total_comparisons > 0.97)

        # The last one produced is IMGT. Use this to test MSA construction.
        hnumbering, lnumbering, hseqs, lseqs = [], [], [], []

        for seq, numbering in zip(seqs, aligner_results):
            if numbering[2] in ["B", "D"]:
                lnumbering.append(numbering)
                lseqs.append(seq)
            elif numbering[2] in ["A", "G"]:
                hnumbering.append(numbering)
                hseqs.append(seq)

        observed_positions = set()
        for n in hnumbering:
            for pos in n[0]:
                observed_positions.add(pos)

        hpositions, hmsa = aligner.build_msa(hseqs, hnumbering, False)
        lpositions, lmsa = aligner.build_msa(lseqs, lnumbering, False)

        for position_set, msa in [(hpositions, hmsa), (lpositions, lmsa)]:
            for msa_seq in msa:
                self.assertTrue(len(msa_seq) == len(position_set))



    def test_extended_length_handling(self):
        """Check situations where the sequence is much longer than
        a typical variable chain. For now this test is only for
        mAbs."""
        project_path = os.path.abspath(os.path.dirname(__file__))
        current_dir = os.getcwd()
        os.chdir(os.path.join(project_path, "test_data"))

        AAs = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N",
            "P", "Q", "R", "S", "T", "V", "W", "Y"]

        with gzip.open("test_data.csv.gz", "rt") as fhandle:
            _ = fhandle.readline()
            seqs = [line.strip().split(",")[0] for line in fhandle]

        os.chdir(current_dir)
        random.seed(123)

        aligner = SingleChainAnnotator(chains=["H", "K", "L"], scheme="imgt")
        for seq in seqs:
            alignment = aligner.analyze_seq(seq)
            # Exclude sequences that have gaps at the ends
            if alignment[-1].startswith('Unexpected AA'):
                continue
            if '120' not in alignment[0]:
                continue
            if alignment[2] == "H":
                if '128' not in alignment[0]:
                    continue
                if seq[alignment[0].index('128')] == "G":
                    continue
            else:
                if '127' not in alignment[0]:
                    continue

            muddled_seq = seq + "YYYGGGGGSSSS" + "".join([random.choice(AAs) for i in range(250)])
            alt_alignment = aligner.analyze_seq(muddled_seq)
            self.assertTrue(len(alt_alignment[0]) == len(muddled_seq))


            trimmed_alt = aligner.trim_alignment(muddled_seq, alt_alignment)
            trimmed_gt = aligner.trim_alignment(seq, alignment)
            self.assertTrue(trimmed_alt[0] == trimmed_gt[0])



    def test_x_handling(self):
        """Check situations where one or more letters has
        been replaced with X. This is specifically for mAbs."""
        project_path = os.path.abspath(os.path.dirname(__file__))
        current_dir = os.getcwd()
        os.chdir(os.path.join(project_path, "test_data"))

        AAs = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N",
            "P", "Q", "R", "S", "T", "V", "W", "Y"]

        with gzip.open("test_data.csv.gz", "rt") as fhandle:
            _ = fhandle.readline()
            seqs = [line.strip().split(",")[0] for line in fhandle]

        os.chdir(current_dir)
        random.seed(123)

        aligner = SingleChainAnnotator(chains=["H", "K", "L"], scheme="imgt")
        for seq in seqs:
            pre_x_alignment = aligner.analyze_seq(seq)
            if pre_x_alignment[-1] != "":
                continue
            if '1' not in pre_x_alignment[0] or '73' not in pre_x_alignment:
                continue

            xposition = random.randint(0, len(seq) - 1)
            seq_list = list(seq)
            seq_list[xposition] = "X"
            alignment = aligner.analyze_seq("".join(seq_list))

            self.assertTrue(alignment[0] == pre_x_alignment[0])



def compare_results(results, comparator_numbering, seqs, scheme):
    """Compares the numbering generated by AntPack with a comparator,
    and returns the number correct vs total comparisons. Also contains
    some (commented-out) code for writing non-matching results to
    a temporary file for closer inspection."""
    total_comparisons, num_correct = 0, 0
    #outhandle = open(f"temp_{scheme}.csv", "w+", encoding="utf-8")
    for result, comparator, seq in zip(results, comparator_numbering, seqs):
        if result[3] != '':
            continue
        total_comparisons += 1
        if result[0] == comparator:
            num_correct += 1
        #This code writes results to file in a format which is easy to look
        #at. Comment it out normally.
        #else:
        #    for i, resnum in enumerate(result[0]):
        #        if resnum != comparator[i]:
        #            result[0][i] = result[0][i] + "!!"
        #    outhandle.write(f"Sequence,{','.join(list(seq))}\n")
        #    outhandle.write(f"Antpack_result,{','.join(result[0])}\n")
        #    outhandle.write(f"Comparator_result,{','.join(comparator)}\n\n")

    #outhandle.close()
    return total_comparisons, num_correct


if __name__ == "__main__":
    unittest.main()
