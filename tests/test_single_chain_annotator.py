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
        return
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

        with self.assertRaises(RuntimeError):
            _ = aligner.sort_position_codes(["a", "1"])
        with self.assertRaises(RuntimeError):
            _ = aligner.sort_position_codes(["1", "2", "C3"])

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
        return
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
        return
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
        return
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



    def test_alignment_trimming(self):
        """Make sure the alignment trimming procedure yields correct results.
        Since this procedure is the same for mabs and TCRs, we do not need
        to check separately."""
        return
        project_path = os.path.abspath(os.path.dirname(__file__))
        current_dir = os.getcwd()
        os.chdir(os.path.join(project_path, "test_data"))

        with gzip.open("test_data.csv.gz", "rt") as fhandle:
            _ = fhandle.readline()
            seqs = [line.strip().split(",")[0] for line in fhandle]

        padded_seqs = []
        AAs = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N",
            "P", "Q", "R", "S", "T", "V", "W", "Y"]

        random.seed(123)

        for seq in seqs:
            left_padding = random.randint(0,10)
            right_padding = random.randint(0,10)
            padded_seq = "".join([random.choice(AAs) for i in range(left_padding)]) + \
                    seq + "".join([random.choice(AAs) for i in range(right_padding)])
            padded_seqs.append(padded_seq)

        os.chdir(current_dir)
        aligner = SingleChainAnnotator(chains=["H", "K", "L"],
                    scheme="imgt")
        aligner_results = aligner.analyze_seqs(padded_seqs)

        for sequence, alignment in zip(padded_seqs, aligner_results):
            numbering = alignment[0]
            exstart = next((i for i in range(len(numbering)) if numbering[i] != '-'), 0)
            exend = next((i for i in range(exstart + 1, len(numbering)) if numbering[i] == '-'),
                    len(numbering))
            trimmed_sequence = sequence[exstart:exend]
            trimmed_numbering = numbering[exstart:exend]

            eval_seq, eval_numbering, eval_start, eval_end = aligner.trim_alignment(sequence,
                    alignment)

            self.assertTrue(eval_seq == trimmed_sequence)
            self.assertTrue(eval_seq == sequence[eval_start:eval_end])
            self.assertTrue(eval_numbering == trimmed_numbering)



    def test_region_labeling(self):
        """Ensure that the region labels assigned by the region labeling
        procedure correspond to our expectations, using a fairly
        inefficient procedure to determine ground-truth labeling."""
        return
        regex = re.compile(r"^(?P<numbers>\d*)(?P<letters>\w*)$")

        project_path = os.path.abspath(os.path.dirname(__file__))
        os.chdir(os.path.join(project_path, "test_data"))

        for receptor, testfile in [("mab", "test_data.csv.gz"),
                ("tcr", "tcr_test_data.csv.gz")]:
            with gzip.open(testfile, "rt") as fhandle:
                _ = fhandle.readline()
                seqs = [line.strip().split(",")[0] for line in fhandle]

            imgt_labels = [(str(i), "fmwk1") for i in range(1,27)] + \
                [(str(i), "cdr1") for i in range(27,40)] + \
                [(str(i), "fmwk2") for i in range(39,56)] + \
                [(str(i), "cdr2") for i in range(56,66)] + \
                [(str(i), "fmwk3") for i in range(66,105)] + \
                [(str(i), "cdr3") for i in range(105,118)] + \
                [(str(i), "fmwk4") for i in range(118,129)]
            imgt_labels = {a:k for (a,k) in imgt_labels}

            martin_heavy = [(str(i), "fmwk1") for i in range(1,26)] + \
                [(str(i), "cdr1") for i in range(26,33)] + \
                [(str(i), "fmwk2") for i in range(33,52)] + \
                [(str(i), "cdr2") for i in range(52,57)] + \
                [(str(i), "fmwk3") for i in range(57,95)] + \
                [(str(i), "cdr3") for i in range(95,103)] + \
                [(str(i), "fmwk4") for i in range(103,114)]
            martin_heavy = {a:k for (a,k) in martin_heavy}
            martin_light = [(str(i), "fmwk1") for i in range(1,26)] + \
                [(str(i), "cdr1") for i in range(26,33)] + \
                [(str(i), "fmwk2") for i in range(33,50)] + \
                [(str(i), "cdr2") for i in range(50,53)] + \
                [(str(i), "fmwk3") for i in range(53,91)] + \
                [(str(i), "cdr3") for i in range(91,97)] + \
                [(str(i), "fmwk4") for i in range(97,108)]
            martin_light = {a:k for (a,k) in martin_light}

            kabat_heavy = [(str(i), "fmwk1") for i in range(1,31)] + \
                [(str(i), "cdr1") for i in range(31,36)] + \
                [(str(i), "fmwk2") for i in range(36,50)] + \
                [(str(i), "cdr2") for i in range(50,66)] + \
                [(str(i), "fmwk3") for i in range(66,95)] + \
                [(str(i), "cdr3") for i in range(95,103)] + \
                [(str(i), "fmwk4") for i in range(103,114)]
            kabat_heavy = {a:k for (a,k) in kabat_heavy}
            kabat_light = [(str(i), "fmwk1") for i in range(1,24)] + \
                [(str(i), "cdr1") for i in range(24,35)] + \
                [(str(i), "fmwk2") for i in range(35,50)] + \
                [(str(i), "cdr2") for i in range(50,57)] + \
                [(str(i), "fmwk3") for i in range(57,89)] + \
                [(str(i), "cdr3") for i in range(89,98)] + \
                [(str(i), "fmwk4") for i in range(98,108)]
            kabat_light = {a:k for (a,k) in kabat_light}



            scheme_labels = {"imgt":{k:imgt_labels for k in
                ["H", "K", "L", "A", "B", "D", "G"]},
                "martin":{"H":martin_heavy, "K":martin_light,
                    "L":martin_light},
                "kabat":{"H":kabat_heavy, "K":kabat_light,
                    "L":kabat_light}
                }


            def get_gt_regions(numbering, label_map):
                gt_reg = []
                for n in numbering:
                    if n == "-":
                        gt_reg.append("-")
                    else:
                        gt_reg.append(label_map[regex.search(n).groups()[0]])
                return gt_reg

            for scheme in ["imgt", "martin", "kabat"]:
                if receptor == "tcr":
                    if scheme != "imgt":
                        continue
                    aligner = SingleChainAnnotator(
                        chains=["A", "B", "D", "G"],
                        scheme=scheme)
                else:
                    aligner = SingleChainAnnotator(chains=["H", "K", "L"],
                        scheme=scheme)
                num_err = 0

                for seq in seqs:
                    numbering = aligner.analyze_seq(seq)
                    labels = aligner.assign_cdr_labels(numbering[0],
                            numbering[2])

                    gt_regions = get_gt_regions(numbering[0],
                            scheme_labels[scheme][numbering[2]])
                    if gt_regions != labels:
                        num_err += 1
                self.assertTrue(num_err == 0)



    def test_cross_scheme_cdr_assignment(self):
        """For antibodies, there is the possibility of numbering
        using one scheme then assigning cdr labels using another.
        This module tests whether this procedure works as
        expected, by comparing it to a slow but easily
        debugged ground-truth assignment. Since there are
        many such possible conversions, we test the most
        popular rather than exhaustively testing all of
        them."""
        regex = re.compile(r"^(?P<numbers>\d*)(?P<letters>\w*)$")

        project_path = os.path.abspath(os.path.dirname(__file__))
        os.chdir(os.path.join(project_path, "test_data"))

        # In the following, the key is
        # numbering-scheme_chain_cdr-scheme
        scheme_labels = {
                "kabat_H_imgt":dict(
                [(str(i), "fmwk1") for i in range(1,26)] + \
                [(str(i), "cdr1") for i in range(26,36)] + \
                [(str(i), "fmwk2") for i in range(36,51)] + \
                [(str(i), "cdr2") for i in range(51,58)] + \
                [(str(i), "fmwk3") for i in range(58,93)] + \
                [(str(i), "cdr3") for i in range(93,103)] + \
                [(str(i), "fmwk4") for i in range(103,114)]
                ),
                "kabat_L_imgt":dict(
                [(str(i), "fmwk1") for i in range(1,27)] + \
                [(str(i), "cdr1") for i in range(27,33)] + \
                [(str(i), "fmwk2") for i in range(33,50)] + \
                [(str(i), "cdr2") for i in range(50,53)] + \
                [(str(i), "fmwk3") for i in range(53,89)] + \
                [(str(i), "cdr3") for i in range(89,98)] + \
                [(str(i), "fmwk4") for i in range(98,108)]
                ),
                "imgt_H_kabat":dict(
                [(str(i), "fmwk1") for i in range(1,32)] + \
                [(str(i), "cdr1") for i in range(32,41)] + \
                [(str(i), "fmwk2") for i in range(41,55)] + \
                [(str(i), "cdr2") for i in range(55,75)] + \
                [(str(i), "fmwk3") for i in range(75,107)] + \
                [(str(i), "cdr3") for i in range(107,118)] + \
                [(str(i), "fmwk4") for i in range(118,129)]
                ),
                "imgt_L_kabat":dict(
                [(str(i), "fmwk1") for i in range(1,24)] + \
                [(str(i), "cdr1") for i in range(24,41)] + \
                [(str(i), "fmwk2") for i in range(41,56)] + \
                [(str(i), "cdr2") for i in range(56,70)] + \
                [(str(i), "fmwk3") for i in range(70,105)] + \
                [(str(i), "cdr3") for i in range(105,118)] + \
                [(str(i), "fmwk4") for i in range(118,129)]
                ),
                "imgt_H_north":dict(
                [(str(i), "fmwk1") for i in range(1,24)] + \
                [(str(i), "cdr1") for i in range(24,41)] + \
                [(str(i), "fmwk2") for i in range(41,55)] + \
                [(str(i), "cdr2") for i in range(55,67)] + \
                [(str(i), "fmwk3") for i in range(67,105)] + \
                [(str(i), "cdr3") for i in range(105,118)] + \
                [(str(i), "fmwk4") for i in range(118,129)]
                ),
                "imgt_L_north":dict(
                [(str(i), "fmwk1") for i in range(1,24)] + \
                [(str(i), "cdr1") for i in range(24,41)] + \
                [(str(i), "fmwk2") for i in range(41,55)] + \
                [(str(i), "cdr2") for i in range(55,70)] + \
                [(str(i), "fmwk3") for i in range(70,105)] + \
                [(str(i), "cdr3") for i in range(105,118)] + \
                [(str(i), "fmwk4") for i in range(118,129)]
                ),
                "kabat_H_north":dict(
                [(str(i), "fmwk1") for i in range(1,23)] + \
                [(str(i), "cdr1") for i in range(23,36)] + \
                [(str(i), "fmwk2") for i in range(36,50)] + \
                [(str(i), "cdr2") for i in range(50,59)] + \
                [(str(i), "fmwk3") for i in range(59,93)] + \
                [(str(i), "cdr3") for i in range(93,103)] + \
                [(str(i), "fmwk4") for i in range(103,114)]
                ),
                "kabat_L_north":dict(
                [(str(i), "fmwk1") for i in range(1,24)] + \
                [(str(i), "cdr1") for i in range(24,35)] + \
                [(str(i), "fmwk2") for i in range(35,49)] + \
                [(str(i), "cdr2") for i in range(49,57)] + \
                [(str(i), "fmwk3") for i in range(57,89)] + \
                [(str(i), "cdr3") for i in range(89,98)] + \
                [(str(i), "fmwk4") for i in range(98,108)]
                ),
                "imgt_H_aho":dict(
                [(str(i), "fmwk1") for i in range(1,25)] + \
                [(str(i), "cdr1") for i in range(25,39)] + \
                [(str(i), "fmwk2") for i in range(39,56)] + \
                [(str(i), "cdr2") for i in range(56,76)] + \
                [(str(i), "fmwk3") for i in range(76,107)] + \
                [(str(i), "cdr3") for i in range(107,117)] + \
                [(str(i), "fmwk4") for i in range(117,129)]
                ),
                "aho_H_imgt":dict(
                [(str(i), "fmwk1") for i in range(1,27)] + \
                [(str(i), "cdr1") for i in range(27,41)] + \
                [(str(i), "fmwk2") for i in range(41,58)] + \
                [(str(i), "cdr2") for i in range(58,68)] + \
                [(str(i), "fmwk3") for i in range(68,107)] + \
                [(str(i), "cdr3") for i in range(107,139)] + \
                [(str(i), "fmwk4") for i in range(139,150)]
                ),
            }

        with gzip.open("test_data.csv.gz", "rt") as fhandle:
            _ = fhandle.readline()
            seqs = [line.strip().split(",")[0] for line in fhandle]

            def get_gt_regions(numbering, label_map):
                gt_reg = []
                for n in numbering:
                    if n == "-":
                        gt_reg.append("-")
                    else:
                        gt_reg.append(label_map[regex.search(n).groups()[0]])
                return gt_reg

            for scheme in ["imgt", "kabat", "aho", "martin"]:
                for cdr_scheme in ["aho", "imgt", "north",
                        "martin", "kabat"]:
                    if scheme == cdr_scheme:
                        continue
                    if f"{scheme}_H_{cdr_scheme}" not in scheme_labels:
                        continue

                    aligner = SingleChainAnnotator(chains=["H", "K", "L"],
                            scheme=scheme)
                    num_err = 0

                    for seq in seqs:
                        numbering = aligner.analyze_seq(seq)
                        labels = aligner.assign_cdr_labels(numbering[0],
                            numbering[2], scheme=cdr_scheme)

                        if scheme in ("aho", "imgt") and \
                                cdr_scheme in ("aho", "imgt"):
                            label_key = f"{scheme}_H_{cdr_scheme}"
                        else:
                            if numbering[2] in ("K", "L"):
                                label_key = f"{scheme}_L_{cdr_scheme}"
                            else:
                                label_key = f"{scheme}_H_{cdr_scheme}"

                        gt_regions = get_gt_regions(numbering[0],
                                scheme_labels[label_key])
                        if gt_regions != labels:
                            import pdb
                            pdb.set_trace()
                            num_err += 1
                    self.assertTrue(num_err == 0)



    def test_position_code_sorting(self):
        """Checks the position code sorting function to make
        sure it is sorting positions correctly for different schemes.
        This procedure is the same for mabs and tcrs, so no need
        to test separately."""
        return
        project_path = os.path.abspath(os.path.dirname(__file__))
        current_dir = os.getcwd()
        os.chdir(os.path.join(project_path, "test_data"))

        with gzip.open("test_data.csv.gz", "rt") as fhandle:
            _ = fhandle.readline()
            seqs = [line.strip().split(",")[0] for line in fhandle]

        os.chdir(current_dir)
        random.seed(123)

        num_err = 0

        for scheme in ["martin", "imgt", "kabat", "aho"]:
            aligner = SingleChainAnnotator(chains=["H", "K", "L"], scheme=scheme)
            for seq in seqs:
                numbering = aligner.analyze_seq(seq)[0]
                numbering = [n for n in numbering if n != "-"]
                shuffled_numbering = copy.deepcopy(numbering)
                random.shuffle(shuffled_numbering)
                sorted_numbering = aligner.sort_position_codes(shuffled_numbering)
                if numbering != sorted_numbering:
                    num_err += 1

        self.assertTrue(num_err == 0)



    def test_extended_length_handling(self):
        """Check situations where the sequence is much longer than
        a typical variable chain. For now this test is only for
        mAbs."""
        return
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
        return
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
