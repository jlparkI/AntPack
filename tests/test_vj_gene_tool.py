"""Tests basic functionality for the VJGeneTool."""
import os
import copy
import gzip
import unittest
from antpack import VJGeneTool, SingleChainAnnotator

class TestVJGeneTool(unittest.TestCase):


    def test_error_checking(self):
        """Check that sequences which have known issues are flagged
        as such, and that deliberately invalid inputs are recognized."""
        # Pass dummy sequences with errors.
        vj_tool = VJGeneTool()
        with self.assertRaises(RuntimeError):
            vj_tool.assign_vj_genes((["1", "2", "3"], 0, "H", ""),
                    "AYAYAYA", "human", "identity")
        with self.assertRaises(RuntimeError):
            vj_tool.assign_vj_genes((["1", "2", "3"], 0, "H", ""),
                    "AYA", "platypus", "identity")
        with self.assertRaises(RuntimeError):
            vj_tool.assign_vj_genes((["1", "2", "3"], 0, "H", ""),
                    "AYA", "human", "robot")


    def test_gene_retrieval(self):
        """Make sure that we can retrieve genes and gene families
        using the appropriate methods."""
        vj_tool = VJGeneTool()
        seq = vj_tool.get_vj_gene_sequence("IGHV2-26*01", "human")
        self.assertTrue(seq == "QVTLKESGP-VLVKPTETLTLTCTVSGFSLS--NARMGVSWIRQPPGKALEWLAHIFSN---DEKSYSTSLK-SRLTISKDTSKSQVVLTMTNMDPVDTATYYCARI---------------------")

        #family_seqs, family_names = vj_tool.get_vj_gene_family("IGHV1", species="human")
        #self.assertTrue(len(family_seqs) == len(family_names))
        #self.assertTrue(len(family_seqs) > 45)
        #self.assertTrue(len([f for f in family_names if not f.startswith("IGHV1")])==0)

        seq = vj_tool.get_vj_gene_sequence("cow", "human")
        self.assertTrue(seq == "")
        #family_seqs, family_names = vj_tool.get_vj_gene_family("cow", species="human")
        #self.assertTrue(len(family_seqs) == 0)
        #self.assertTrue(len(family_names) == 0)


    def test_percent_ident_calc(self):
        """Double checks the percent identity calculation
        done by the cpp extension using a simple stupid
        Python version, just to make sure it's correct."""
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

            vpred, jpred, videntity, jidentity = vj_tool.assign_vj_genes(alignment,
                    seq, "human", "identity")

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
                self.assertTrue(gtrue == gpred)
                self.assertTrue(best_pid == id_pred)


    def test_vj_assignment(self):
        """Checks vj assignments against those done by other
        tools to ensure that they are usually the same."""
        vj_tool = VJGeneTool()
        sc_annotator = SingleChainAnnotator()

        project_path = os.path.abspath(os.path.dirname(__file__))
        current_dir = os.getcwd()
        os.chdir(os.path.join(project_path, "test_data"))

        vhmatches, vklmatches, vhtests, vkltests = 0, 0, 0, 0
        jmatches, ntests = 0, 0

        with gzip.open("vj_gene_testing.csv.gz", "rt") as fhandle:
            _ = fhandle.readline()

            for line in fhandle:
                seq, vgene, jgene = line.strip().split(",")
                # Eliminate problematic sequences (e.g. missing cysteine).
                alignment = sc_annotator.analyze_seq(seq)
                pid, chain, err = alignment[1], alignment[2], alignment[3]
                if pid < 0.8 or err != "":
                    continue
                pred_vgene, pred_jgene, pidv, pidj = vj_tool.assign_vj_genes(alignment,
                        seq, "human", "evalue")

                if pred_vgene == vgene:
                    if chain == "H":
                        vhmatches += 1
                    else:
                        vklmatches += 1
                if pred_jgene == jgene:
                    jmatches += 1
                if chain == "H":
                    vhtests += 1
                else:
                    vkltests += 1
                ntests += 1

        print(f"On {vhtests}, vhgene, {vhmatches} success.")
        print(f"On {vkltests}, vklgene, {vklmatches} success.")
        print(f"On {ntests}, jgene, {jmatches} success.")

        self.assertTrue((vhmatches / vhtests) > 0.9)
        self.assertTrue((vklmatches / vkltests) > 0.9)
        self.assertTrue((jmatches / ntests) > 0.9)

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
