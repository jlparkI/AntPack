"""Tests basic functionality for the VJGeneTool."""
import os
import copy
import gzip
import unittest
import numpy as np
import blosum as bl
from antpack import VJGeneTool, SingleChainAnnotator

class TestVJGeneTool(unittest.TestCase):


    def test_error_checking(self):
        """Check that sequences which have known issues are flagged
        as such, and that deliberately invalid inputs are recognized."""
        # Pass dummy sequences with errors.
        vj_tool = VJGeneTool()
        vgene, jgene, vident, jident, species = vj_tool.assign_vj_genes(
                (["1", "2", "3"], 0, "H", ""),
                "AYAYAYA", "human")
        self.assertTrue(vident==0)
        self.assertTrue(jident==0)
        self.assertTrue(vgene=="")
        self.assertTrue(jgene=="")
        self.assertTrue(species=="unknown")

        with self.assertRaises(RuntimeError):
            tcr_tool = VJGeneTool(scheme = "aho")
            tcr_tool.assign_vj_genes(
                (["1", "2", "3"], 0, "A", ""),
                "AYA", "human")

        with self.assertRaises(RuntimeError):
            vj_tool.assign_vj_genes(
                (["1", "2", "3"], 0, "H", ""),
                "AYA", "platypus")

        with self.assertRaises(RuntimeError):
            vj_tool.assign_vj_genes(
                (["1", "2", "3"], 0, "H", ""),
                "AYA", "human", "robot")

        with self.assertRaises(RuntimeError):
            vj_tool.assign_vj_genes(
                (["1", "2", "3"], 0, "L", ""),
                "AYA", "alpaca")

        with self.assertRaises(RuntimeError):
            vj_tool.assign_vj_genes(
                (["1", "2", "3"], 0, "K", ""),
                "AYA", "alpaca", "identity")

        with self.assertRaises(RuntimeError):
            vj_tool.assign_vj_genes(
                (["1", "2", "3"], 0, "A", ""),
                "AYA", "alpaca", "identity")


    def test_gene_retrieval(self):
        """Make sure that we can retrieve genes and gene families
        using the appropriate methods."""
        vj_tool = VJGeneTool()
        seq = vj_tool.get_vj_gene_sequence("IGHV2-26*01", "human")
        self.assertTrue(seq == "QVTLKESGP-VLVKPTETLTLTCTVSGFSLS--NARMGVSWIRQPPGKALEWLAHIFSN---DEKSYSTSLK-SRLTISKDTSKSQVVLTMTNMDPVDTATYYCARI---------------------")

        seq = vj_tool.get_vj_gene_sequence("IGHV1-1*01", "alpaca")
        self.assertTrue(seq == "QVQLVQPGA-ELRKPGALLKVSCKASGYTF----TSYYIDWVRQAPGQGLGWVGRIDPE--DGGTNYAQKFQ-GRVTLTADTSTSTAYVELSSLRSEDTAVCYCVR----------------------")

        seq = vj_tool.get_vj_gene_sequence("IGHV1S26*01", "rabbit")
        self.assertTrue(seq == "Q-SVKESEG-GLFKPTDTLTLTCTVSGFSL----SSYAISWVRQAPGNGLEWIGIINSY---GSTYYASWAK-SRSTITRNTNENTVTLKMTSLTAADTATYFCAR----------------------")

        seq = vj_tool.get_vj_gene_sequence("cow", "human")
        self.assertTrue(seq == "")


    def test_species_recognition(self):
        """Make sure that expected species are all recognized."""
        vj_tool = VJGeneTool()
        test_heavy = "QVQLLESGGGLVQPGGSLRLSCAASGFTFSTAAMSWVRQAPGKGLEWVSGISGSGSSTYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCARELSYLYSGYYFDYWGQGTLVTVSSASTKGPSVFPLAPSSKSTSGGTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLSSVVTVPSSSLGTQTYICNVNHKPSNTKVDKRVEPKSCDKTHTCPPCPAPELLGGPSVFLFPPKPKDTLMISRTPEVTCVVVAVSHEDPEVKFNWYVDGVEVHNAKTKPREEQYNSTYRVVSVLTVLHQDWLNGKEYKCKVSNKALAAPIEKTISKAKGQPREPQVYTLPPSREEMTKNQVSLTCLVKGFYPSDIAVEWESNGQPENNYKTTPPVLDSDGSFFLYSKLTVDKSRWQQGNVFSCSVMHEALHNHYTQKSLSLSPGK"
        test_light = "DIELSQSPAILSASPGEKVTMTCRASSSVSYMHWYQQKPGSSPKPWIYAPSNLASGVPARFSGSGSGTSYSLTISRVEAEDAATYYCQQWSFNPPTFGAGTKLEIKRTVAAPSVFIFPPSDEQLKSGTASVVCLLNNFYPREAKVQWKVDNALQSGNSQESVTEQDSKDSTYSLSSTLTLSKADYEKHKVYACEVTHQGLSSPVTKSFNRGEC"

        sc_annotator = SingleChainAnnotator()
        heavy_annotation = sc_annotator.analyze_seq(test_heavy)
        light_annotation = sc_annotator.analyze_seq(test_light)

        assn = vj_tool.assign_vj_genes(heavy_annotation, test_heavy, "human")
        self.assertTrue(assn[0] != '')
        assn = vj_tool.assign_vj_genes(heavy_annotation, test_heavy, "mouse")
        self.assertTrue(assn[0] != '')
        assn = vj_tool.assign_vj_genes(heavy_annotation, test_heavy, "rabbit")
        self.assertTrue(assn[0] != '')
        assn = vj_tool.assign_vj_genes(heavy_annotation, test_heavy, "alpaca")
        self.assertTrue(assn[0] != '')
        assn = vj_tool.assign_vj_genes(heavy_annotation, test_heavy, "unknown")
        self.assertTrue(assn[0] != '')

        assn = vj_tool.assign_vj_genes(light_annotation, test_light, "human")
        self.assertTrue(assn[0] != '')
        assn = vj_tool.assign_vj_genes(light_annotation, test_light, "mouse")
        self.assertTrue(assn[0] != '')
        assn = vj_tool.assign_vj_genes(light_annotation, test_light, "rabbit")
        self.assertTrue(assn[0] != '')
        assn = vj_tool.assign_vj_genes(light_annotation, test_light,
                "unknown")
        self.assertTrue(assn[0] != '')


    def test_percent_ident_calc(self):
        """Double checks the percent identity calculation
        done by the cpp extension using a simple stupid
        Python version."""
        vj_tool = VJGeneTool()

        project_path = os.path.abspath(os.path.dirname(__file__))
        current_dir = os.getcwd()
        os.chdir(os.path.join(project_path, "test_data"))

        with gzip.open("test_data.csv.gz", "rt") as fhandle:
            _ = fhandle.readline()
            test_seqs = [line.strip().split(",")[0] for line in fhandle]

        germline_seqs, germline_names = vj_tool.get_seq_lists()

        sc_annotator = SingleChainAnnotator()

        os.chdir(current_dir)
        for seq in test_seqs[:100]:
            alignment = sc_annotator.analyze_seq(seq)
            chain = alignment[2]
            fmt_seq = prep_sequence(seq, alignment)

            vpred, jpred, videntity, jidentity, _ = vj_tool.assign_vj_genes(alignment,
                    seq, "human")

            gpreds, gidentities = (vpred, jpred), (videntity, jidentity)

            for recep, gpred, id_pred in zip([f"IG{chain}V", f"IG{chain}J"],
                    gpreds, gidentities):

                best_pid, matchnum = 0, 0
                for i, template in enumerate(germline_seqs[f"human_{recep}"]):
                    matches, ntot = 0, 0

                    for qletter, tletter in zip(fmt_seq, template):
                        if tletter == "-":
                            continue
                        ntot += 1
                        if qletter == tletter:
                            matches += 1

                    true_pid = float(matches) / float(ntot)
                    if true_pid > best_pid:
                        matchnum = i
                        best_pid = copy.deepcopy(true_pid)

                gtrue = germline_names[f"human_{recep}"][matchnum]
                self.assertTrue(gtrue in gpred)
                self.assertTrue(np.allclose(id_pred, best_pid))



    def test_evalue_calc(self):
        """Double checks the evalue calculation
        done by the cpp extension using a
        simple stupid Python version."""
        vj_tool = VJGeneTool()

        blosum_matrix = bl.BLOSUM(62)

        project_path = os.path.abspath(os.path.dirname(__file__))
        current_dir = os.getcwd()
        os.chdir(os.path.join(project_path, "test_data"))

        with gzip.open("test_data.csv.gz", "rt") as fhandle:
            _ = fhandle.readline()
            test_seqs = [line.strip().split(",")[0] for line in fhandle]

        germline_seqs, germline_names = vj_tool.get_seq_lists()

        sc_annotator = SingleChainAnnotator()

        os.chdir(current_dir)
        for seq in test_seqs[:100]:
            alignment = sc_annotator.analyze_seq(seq)
            chain = alignment[2]
            fmt_seq = prep_sequence(seq, alignment)

            vpred, jpred, videntity, jidentity, _ = vj_tool.assign_vj_genes(alignment,
                    seq, "human", "evalue")

            gpreds, gidentities = (vpred, jpred), (videntity, jidentity)

            for recep, gpred, id_pred in zip([f"IG{chain}V", f"IG{chain}J"],
                    gpreds, gidentities):

                best_pid, matchnum = -1e10, 0
                for i, template in enumerate(germline_seqs[f"human_{recep}"]):
                    score = 0

                    for qletter, tletter in zip(fmt_seq, template):
                        if tletter == "-":
                            continue
                        if qletter == "-":
                            continue

                        score += blosum_matrix[qletter][tletter]


                    if score > best_pid:
                        matchnum = i
                        best_pid = copy.deepcopy(score)

                gtrue = germline_names[f"human_{recep}"][matchnum]
                self.assertTrue(gtrue in gpred)
                self.assertTrue(np.allclose(id_pred, best_pid))


    def test_vj_assignment(self):
        """Checks vj assignments against those done by other
        tools to ensure that they are usually the same.
        Also, test whether supplying 'unknown' for species yields
        valid results."""
        project_path = os.path.abspath(os.path.dirname(__file__))
        current_dir = os.getcwd()
        os.chdir(os.path.join(project_path, "test_data"))

        vhmatches, vklmatches, vhtests, vkltests = 0, 0, 0, 0
        jmatches, ntests = 0, 0

        for receptor, fname, species in (("tcr", "tcr_vj_gene_testing.csv.gz", "unknown"),
                ("mab", "vj_gene_testing.csv.gz", "human"),
                ("mab", "vj_gene_testing.csv.gz", "unknown")):
            if receptor == "mab":
                vj_tool = VJGeneTool()
                sc_annotator = SingleChainAnnotator()
            elif receptor == "tcr":
                vj_tool = VJGeneTool()
                sc_annotator = SingleChainAnnotator(["A", "B", "D", "G"])

            with gzip.open(fname, "rt") as fhandle:
                _ = fhandle.readline()

                for line in fhandle:
                    seq, vgene, jgene = line.strip().split(",")
                    # Eliminate problematic sequences (e.g. missing cysteine).
                    alignment = sc_annotator.analyze_seq(seq)
                    pid, chain, err = alignment[1], alignment[2], alignment[3]
                    if pid < 0.8 or err != "":
                        continue

                    pred_vgene, pred_jgene, pidv, pidj, rs = \
                            vj_tool.assign_vj_genes(alignment, seq,
                                    species, "identity")

                    if vgene in pred_vgene:
                        if chain in ["H", "A", "G"]:
                            vhmatches += 1
                        else:
                            vklmatches += 1
                    if jgene in pred_jgene:
                        jmatches += 1
                    if chain in ["H", "A", "G"]:
                        vhtests += 1
                    else:
                        vkltests += 1
                    ntests += 1

            print(f"{receptor}, {species} On {vhtests}, vhgene, {vhmatches} success.")
            print(f"{receptor}, {species} On {vkltests}, vklgene, {vklmatches} success.")
            print(f"{receptor}, {species} On {ntests}, jgene, {jmatches} success.")

            self.assertTrue((vhmatches / vhtests) > 0.9)
            self.assertTrue((vklmatches / vkltests) > 0.9)
            self.assertTrue((jmatches / ntests) > 0.9)



    def test_scheme_switching(self):
        """Checks vj assignments using different schemes to
        ensure they match. We do not need to test this for TCRs
        at this time since at this time TCRs support IMGT only."""
        imgt_tool = VJGeneTool(scheme="imgt")
        imgt_aligner = SingleChainAnnotator(scheme="imgt")

        project_path = os.path.abspath(os.path.dirname(__file__))
        current_dir = os.getcwd()
        os.chdir(os.path.join(project_path, "test_data"))

        for alternate_scheme in ["aho", "kabat", "martin"]:
            alternate_tool = VJGeneTool(scheme=alternate_scheme)
            alternate_aligner = SingleChainAnnotator(scheme=alternate_scheme)

            with gzip.open("vj_gene_testing.csv.gz", "rt") as fhandle:
                _ = fhandle.readline()

                for line in fhandle:
                    seq, _, _ = line.strip().split(",")
                    annotation1 = imgt_aligner.analyze_seq(seq)
                    annotation2 = alternate_aligner.analyze_seq(seq)

                    pred_vgene1, pred_jgene1, pidv1, pidj1, species1 = \
                        imgt_tool.assign_vj_genes(annotation1,
                            seq, "human", "identity")
                    pred_vgene2, pred_jgene2, pidv2, pidj2, species2 = \
                        alternate_tool.assign_vj_genes(
                            annotation2, seq, "human", "identity")

                    self.assertTrue(pred_vgene1==pred_vgene2)
                    self.assertTrue(pred_jgene1==pred_jgene2)
                    self.assertTrue(np.allclose(pidv1, pidv2))
                    self.assertTrue(np.allclose(pidj1, pidj2))
                    self.assertTrue(species1==species2)

        os.chdir(current_dir)





def prep_sequence(sequence, alignment):
    """Converts an imgt-formatted sequence into a query for
    matching up to VJ genes."""
    numbering = alignment[0]
    std_positions = {str(i+1) for i in range(128)}

    fmt_seq = ['-' for i in range(128)]
    for letter, nb_assign in zip(sequence, numbering):
        if nb_assign in std_positions:
            fmt_seq[int(nb_assign) - 1] = letter

    return "".join(fmt_seq)



if __name__ == "__main__":
    unittest.main()
