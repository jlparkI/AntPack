"""Tests basic functionality for the PairedChainAnnotator class.
Notice that this test relies on the correct functioning of
SingleChainAnnotator, so if test_single_chain_annotator does
not pass, this one will not either."""
import random
import pytest
from antpack import PairedChainAnnotator, SingleChainAnnotator
from antpack.scoring_tools.scoring_constants import scoring_constants as SCCONST

from ..conftest import (standard_aa_list, load_mab_nmbr_test_data,
                        load_tcr_nmbr_test_data)


@pytest.mark.parametrize("receptor_type", ["mab", "tcr"])
def test_error_checking(receptor_type):
    """Check that sequences which have known issues are flagged
    as such, and that deliberately invalid inputs are recognized."""
    # Pass dummy sequences with errors.
    aligner = PairedChainAnnotator(receptor_type = receptor_type)

    for invalid_seq in ["YaY", "YBW", "Y K", "Y-K"]:
        results = aligner.analyze_seq(invalid_seq)
        assert results[0][3].startswith("Sequence contains")
        assert results[1][3].startswith("Sequence contains")



def test_mab_analyze_seqs(load_mab_nmbr_test_data):
    """Test that analyze_seqs returns the same
    results as analyze_seq applied to the list.
    (analyze_seqs calls analyze_seq, so highly
    unlikely it does anything different...but just
    in case)."""
    seqs, _ = load_mab_nmbr_test_data
    sc_aligner = SingleChainAnnotator(scheme="imgt")
    pc_aligner = PairedChainAnnotator(scheme="imgt",
                receptor_type="mab")
    alignments = [sc_aligner.analyze_seq(seq) for seq in seqs]

    heavy_chains = [(seq, a) for a, seq in zip(alignments, seqs)
            if a[2] in ["H", "A", "G"]]
    light_chains = [(seq, a) for a, seq in zip(alignments, seqs)
            if a[2] not in ["H", "A", "G"]]

    # Join the sequences with a very arbitrary linker.
    merged_seqs = [l[0] + "YYYSSSGGG" + h[0] for
            h, l in zip(heavy_chains, light_chains)]

    heavy_test, light_test = pc_aligner.analyze_seqs(merged_seqs)
    assert len(heavy_test)==len(merged_seqs)
    assert len(light_test)==len(merged_seqs)

    for i, seq in enumerate(merged_seqs):
        h, l = pc_aligner.analyze_seq(seq)
        assert h==heavy_test[i]
        assert l==light_test[i]


@pytest.mark.parametrize("receptor_type", ["tcr", "mab"])
def test_analyze_seqs(load_mab_nmbr_test_data,
        load_tcr_nmbr_test_data, receptor_type):
    """Test that analyze_seqs returns the same
    results as analyze_seq applied to the list.
    (analyze_seqs calls analyze_seq, so highly
    unlikely it does anything different...but just
    in case)."""
    if receptor_type == "mab":
        seqs, _ = load_mab_nmbr_test_data
        sc_aligner = SingleChainAnnotator(scheme="imgt")
    else:
        seqs, _ = load_tcr_nmbr_test_data
        sc_aligner = SingleChainAnnotator(["A", "B", "D", "G"],
                    scheme="imgt")

    pc_aligner = PairedChainAnnotator("imgt", receptor_type)
    alignments = [sc_aligner.analyze_seq(seq) for seq in seqs]

    heavy_chains = [(seq, a) for a, seq in zip(alignments, seqs)
            if a[2] in ["H", "A", "G"]]
    light_chains = [(seq, a) for a, seq in zip(alignments, seqs)
            if a[2] not in ["H", "A", "G"]]

    # Join the sequences with a very arbitrary linker.
    merged_seqs = [l[0] + "YYYSSSGGG" + h[0] for
            h, l in zip(heavy_chains, light_chains)]

    heavy_test, light_test = pc_aligner.analyze_seqs(merged_seqs)
    assert len(heavy_test)==len(merged_seqs)
    assert len(light_test)==len(merged_seqs)

    for i, seq in enumerate(merged_seqs):
        h, l = pc_aligner.analyze_seq(seq)
        assert h==heavy_test[i]
        assert l==light_test[i]




@pytest.mark.parametrize("receptor_type", ["mab", "tcr"])
def test_sequence_extraction(load_mab_nmbr_test_data,
                load_tcr_nmbr_test_data, receptor_type):
    """Take test heavy and light chains and number them using
    SingleChainAnnotator. Next, join pairs of chains with a random
    number of random amino acids and add a random number of
    random amino acids at the beginning and end. Compare the
    numbering assigned by SingleChainAnnotator to that assigned
    by PairedChainAnnotator (this of course assumes that
    SingleChainAnnotator is functioning correctly)."""
    if receptor_type == "mab":
        seqs, _ = load_mab_nmbr_test_data
        sc_aligner = SingleChainAnnotator(scheme="imgt")
    else:
        seqs, _ = load_tcr_nmbr_test_data
        sc_aligner = SingleChainAnnotator(["A", "B", "D", "G"],
                    scheme="imgt")

    m_aligner = PairedChainAnnotator("imgt", receptor_type)

    alignments = [sc_aligner.analyze_seq(seq) for seq in seqs]

    heavy_chains = [(seq, a) for a, seq in zip(alignments, seqs)
            if a[2] in ["H", "A", "G"]]
    light_chains = [(seq, a) for a, seq in zip(alignments, seqs)
            if a[2] not in ["H", "A", "G"]]

    random.seed(0)

    for hc, lc in zip(heavy_chains, light_chains):
        prefixes = [[SCCONST.aa_list[random.randint(0,19)]
            for i in range(random.randint(5,25))]
            for j in range(4)]

        # Try with both the heavy chain first and the light chain
        # first.
        merged_hc = "".join( ["".join(prefixes[0]), hc[0], "".join(prefixes[1]) ] )
        merged_lc = "".join( ["".join(prefixes[2]), lc[0], "".join(prefixes[3]) ] )
        merged_chains = [merged_hc + merged_lc, merged_lc + merged_hc]

        hc_align = sc_aligner.analyze_seq(merged_hc)
        lc_align = sc_aligner.analyze_seq(merged_lc)
        if hc_align[1] < 0.7 or lc_align[1] < 0.7:
            continue

        if hc_align[-1] != "" or lc_align[-1] != "":
            continue

        for i, lcpos in enumerate(lc_align[0]):
            if lcpos != "-":
                break
        for j, lcpos in reversed(list(enumerate(lc_align[0]))):
            if lcpos != "-":
                break

        trimmed_lc_seq = merged_lc[i:j+1]
        trimmed_lc_align = lc_align[0][i:j+1]
        if not trimmed_lc_align[0] == "1":
            continue

        for i, hcpos in enumerate(hc_align[0]):
            if hcpos != "-":
                break
        for j, hcpos in reversed(list(enumerate(hc_align[0]))):
            if hcpos != "-":
                break

        trimmed_hc_seq = merged_hc[i:j+1]
        trimmed_hc_align = hc_align[0][i:j+1]

        for merged_chain in merged_chains:
            mc_heavy, mc_light = m_aligner.analyze_seq(merged_chain)
            assert len(mc_heavy[0]) == len(merged_chain)
            assert len(mc_light[0]) == len(merged_chain)

            _, mchn, hstart, hend = m_aligner.trim_alignment(merged_chain, mc_heavy)
            _, mcln, lstart, lend = m_aligner.trim_alignment(merged_chain, mc_light)

            # Skip sequences where the first letter was missing from the heavy
            # chain. These can sometimes cause variable gapping depending
            # on how many prefix letters were added.
            if mchn[0] != "1" or mchn[1] == "3":
                if mchn != trimmed_hc_align:
                    continue

            assert mcln == trimmed_lc_align
            assert mchn == trimmed_hc_align
            assert merged_chain[lstart:lend] == trimmed_lc_seq
            assert merged_chain[hstart:hend] == trimmed_hc_seq


@pytest.mark.parametrize("receptor_type", ["mab", "tcr"])
def test_single_chain_behavior(load_mab_nmbr_test_data,
        load_tcr_nmbr_test_data, receptor_type):
    """Take known heavy / light chains only and run them 
    through PairedChainAnnotator. Compare the
    numbering assigned by SingleChainAnnotator to that assigned
    by PairedChainAnnotator (this of course assumes that
    SingleChainAnnotator is functioning correctly)."""
    if receptor_type == "mab":
        seqs, _ = load_mab_nmbr_test_data
        sc_aligner = SingleChainAnnotator(scheme="imgt")
    else:
        seqs, _ = load_tcr_nmbr_test_data
        sc_aligner = SingleChainAnnotator(["A", "B", "D", "G"],
                    scheme="imgt")

    m_aligner = PairedChainAnnotator("imgt", receptor_type)

    # Skip sequences with only one cysteine, which are hard
    # to analyze.
    seqs = [s for s in seqs if len([l for l in s if l == "C"])
            >= 2]

    alignments = [sc_aligner.analyze_seq(seq) for seq in seqs]

    heavy_chains = [(seq, a) for a, seq in zip(alignments, seqs)
            if a[2] in ["H", "A", "G"]]
    light_chains = [(seq, a) for a, seq in zip(alignments, seqs)
            if a[2] not in ["H", "A", "G"]]

    random.seed(0)

    for (seq, alignment) in heavy_chains:
        pc_result = m_aligner.analyze_seq(seq)[0]
        assert pc_result == alignment

    for (seq, alignment) in light_chains:
        pc_result = m_aligner.analyze_seq(seq)[1]
        assert pc_result == alignment
