"""Tests basic functionality for the VJGeneTool."""
import os
import gzip
import copy
import numpy as np
import pytest
import blosum as bl
from antpack import VJGeneTool, SingleChainAnnotator
from ..conftest import (standard_aa_list, get_test_data_filepath,
                    load_mab_nmbr_test_data, load_tcr_nmbr_test_data)


def test_error_checking():
    """Check that sequences which have known issues are flagged
    as such, and that deliberately invalid inputs are recognized."""
    # Pass dummy sequences with errors.
    vj_tool = VJGeneTool()
    vgene, jgene, vident, jident, species = vj_tool.assign_vj_genes(
            (["1", "2", "3"], 0, "H", ""),
            "AYAYAYA", "human")
    assert vident==0
    assert jident==0
    assert vgene==""
    assert jgene==""
    assert species=="unknown"

    with pytest.raises(RuntimeError):
        tcr_tool = VJGeneTool(scheme = "aho")
        tcr_tool.assign_vj_genes(
            (["1", "2", "3"], 0, "A", ""),
            "AYA", "human")

    with pytest.raises(RuntimeError):
        vj_tool.assign_vj_genes(
            (["1", "2", "3"], 0, "H", ""),
            "AYA", "platypus")

    with pytest.raises(RuntimeError):
        vj_tool.assign_vj_genes(
            (["1", "2", "3"], 0, "H", ""),
            "AYA", "human", "robot")

    with pytest.raises(RuntimeError):
        vj_tool.assign_vj_genes(
            (["1", "2", "3"], 0, "L", ""),
            "AYA", "alpaca")

    with pytest.raises(RuntimeError):
        vj_tool.assign_vj_genes(
            (["1", "2", "3"], 0, "K", ""),
            "AYA", "alpaca", "identity")

    with pytest.raises(RuntimeError):
        vj_tool.assign_vj_genes(
            (["1", "2", "3"], 0, "A", ""),
            "AYA", "alpaca", "identity")


def test_gene_retrieval(return_known_vgenes):
    """Make sure that we can retrieve genes and gene families
    using the appropriate methods."""
    vj_tool = VJGeneTool()
    for vgene_test in return_known_vgenes:
        seq = vj_tool.get_vj_gene_sequence(vgene_test[0], vgene_test[1])
        assert seq == vgene_test[2]

    seq = vj_tool.get_vj_gene_sequence("cow", "human")
    assert seq == ""


def test_species_recognition(species_test_sequences):
    """Make sure that expected species are all recognized."""
    vj_tool = VJGeneTool()
    test_heavy, test_light = species_test_sequences

    sc_annotator = SingleChainAnnotator()
    heavy_annotation = sc_annotator.analyze_seq(test_heavy)
    light_annotation = sc_annotator.analyze_seq(test_light)

    for species in ["human", "mouse", "rabbit", "alpaca", "unknown"]:
        assn = vj_tool.assign_vj_genes(heavy_annotation, test_heavy, species)
        assert assn[0] != ''

    for species in ["human", "mouse", "rabbit", "unknown"]:
        assn = vj_tool.assign_vj_genes(light_annotation, test_light, species)
        assert assn[0] != ''


def test_percent_ident_calc(load_mab_nmbr_test_data):
    """Double checks the percent identity calculation
    done by the cpp extension using a simple stupid
    Python version."""
    vj_tool = VJGeneTool()
    test_seqs, _ = load_mab_nmbr_test_data
    germline_seqs, germline_names = vj_tool.get_seq_lists()
    sc_annotator = SingleChainAnnotator()

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
            assert gtrue in gpred
            assert np.allclose(id_pred, best_pid)



def test_evalue_calc(load_mab_nmbr_test_data):
    """Double checks the evalue calculation
    done by the cpp extension using a
    simple stupid Python version."""
    vj_tool = VJGeneTool()
    test_seqs, _ = load_mab_nmbr_test_data
    blosum_matrix = bl.BLOSUM(62)

    germline_seqs, germline_names = vj_tool.get_seq_lists()
    sc_annotator = SingleChainAnnotator()

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

                    if qletter == tletter:
                        score += 4
                    else:
                        score += blosum_matrix[qletter][tletter]


                if score > best_pid:
                    matchnum = i
                    best_pid = copy.deepcopy(score)

            gtrue = germline_names[f"human_{recep}"][matchnum]
            assert gtrue in gpred
            assert np.allclose(id_pred, best_pid)


def test_vj_assignment(get_test_data_filepath):
    """Checks vj assignments against those done by other
    tools to ensure that they are usually the same.
    Also, test whether supplying 'unknown' for species yields
    valid results."""
    vhmatches, vklmatches, vhtests, vkltests = 0, 0, 0, 0
    jmatches, ntests = 0, 0

    for receptor, fname, species in (("tcr", "tcr_vj_gene_testing.csv.gz", "unknown"),
            ("mab", "vj_gene_testing.csv.gz", "human"),
            ("mab", "vj_gene_testing.csv.gz", "unknown")):
        if receptor == "mab":
            vj_tool = VJGeneTool()
            sc_annotator = SingleChainAnnotator()
        else:
            vj_tool = VJGeneTool()
            sc_annotator = SingleChainAnnotator(["A", "B", "D", "G"])

        with gzip.open(os.path.join(get_test_data_filepath, fname),
                       "rt") as fhandle:
            _ = fhandle.readline()

            for line in fhandle:
                seq, vgene, jgene = line.strip().split(",")
                # Eliminate problematic sequences (e.g. missing cysteine).
                alignment = sc_annotator.analyze_seq(seq)
                pid, chain, err = alignment[1], alignment[2], alignment[3]
                if pid < 0.8 or err != "":
                    continue

                pred_vgene, pred_jgene, _, _, _ = \
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

        print(f"{receptor}, {species} On {vhtests}, vhgene, {vhmatches} success.", flush=True)
        print(f"{receptor}, {species} On {vkltests}, vklgene, {vklmatches} success.", flush=True)
        print(f"{receptor}, {species} On {ntests}, jgene, {jmatches} success.", flush=True)

        assert (vhmatches / vhtests) > 0.9
        assert (vklmatches / vkltests) > 0.9
        assert (jmatches / ntests) > 0.9



@pytest.mark.parametrize("alternate_scheme", ["aho", "martin", "kabat"])
def test_scheme_switching(get_test_data_filepath,
            alternate_scheme):
    """Checks vj assignments using different schemes to
    ensure they match. We do not need to test this for TCRs
    at this time since at this time TCRs support IMGT only."""
    imgt_tool = VJGeneTool(scheme="imgt")
    imgt_aligner = SingleChainAnnotator(scheme="imgt")

    alternate_tool = VJGeneTool(scheme=alternate_scheme)
    alternate_aligner = SingleChainAnnotator(scheme=alternate_scheme)

    with gzip.open(os.path.join(get_test_data_filepath,
                        "vj_gene_testing.csv.gz"), "rt") as fhandle:
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

            assert pred_vgene1==pred_vgene2
            assert pred_jgene1==pred_jgene2
            assert np.allclose(pidv1, pidv2)
            assert np.allclose(pidj1, pidj2)
            assert species1==species2



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



@pytest.fixture
def return_known_vgenes():
    """A list of known vgenes to check that the vj tool
    retrieves vgenes correctly."""
    return (
        ("IGHV2-26*01", "human", ("QVTLKESGP-VLVKPTETLTLTCTVSGFSLS-"
        "-NARMGVSWIRQPPGKALEWLAHIFSN---DEKSYSTSLK-SRLTISKDTSKSQVVL"
        "TMTNMDPVDTATYYCARI---------------------")),
        ("IGHV1-1*01", "alpaca", (
        "QVQLVQPGA-ELRKPGALLKVSCKASGYTF----TSYYIDWVRQAPGQGLGWVGRIDPE"
        "--DGGTNYAQKFQ-GRVTLTADTSTSTAYVELSSLRSEDTAVCYCVR----------------------")
         ),
        ("IGHV1S26*01", "rabbit", (
        "Q-SVKESEG-GLFKPTDTLTLTCTVSGFSL----SSYAISWVRQAPGNGLEWIGIINSY"
        "---GSTYYASWAK-SRSTITRNTNENTVTLKMTSLTAADTATYFCAR----------------------")
         ),
        )

@pytest.fixture
def species_test_sequences():
    """A couple of sequences of known origin for testing species
    recognition."""
    test_heavy = "QVQLLESGGGLVQPGGSLRLSCAASGFTFSTAAMSWVRQAPGKGLEWVSGISGSGSSTYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCARELSYLYSGYYFDYWGQGTLVTVSSASTKGPSVFPLAPSSKSTSGGTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLSSVVTVPSSSLGTQTYICNVNHKPSNTKVDKRVEPKSCDKTHTCPPCPAPELLGGPSVFLFPPKPKDTLMISRTPEVTCVVVAVSHEDPEVKFNWYVDGVEVHNAKTKPREEQYNSTYRVVSVLTVLHQDWLNGKEYKCKVSNKALAAPIEKTISKAKGQPREPQVYTLPPSREEMTKNQVSLTCLVKGFYPSDIAVEWESNGQPENNYKTTPPVLDSDGSFFLYSKLTVDKSRWQQGNVFSCSVMHEALHNHYTQKSLSLSPGK"
    test_light = "DIELSQSPAILSASPGEKVTMTCRASSSVSYMHWYQQKPGSSPKPWIYAPSNLASGVPARFSGSGSGTSYSLTISRVEAEDAATYYCQQWSFNPPTFGAGTKLEIKRTVAAPSVFIFPPSDEQLKSGTASVVCLLNNFYPREAKVQWKVDNALQSGNSQESVTEQDSKDSTYSLSSTLTLSKADYEKHKVYACEVTHQGLSSPVTKSFNRGEC"
    return test_heavy, test_light
