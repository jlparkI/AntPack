"""Tests basic functionality for the LiabilitySearchTool class."""
import pytest
from antpack import LiabilitySearchTool, SingleChainAnnotator



def test_liability_tool(generic_test_seq):
    """Check that sequences with liabilities are flagged as such
    while other sequences are passed."""
    search_tool = LiabilitySearchTool()
    annotation_tool = SingleChainAnnotator()

    alignment = annotation_tool.analyze_seq(generic_test_seq)
    result = search_tool.analyze_seq(generic_test_seq, alignment,
            "imgt", "imgt")
    assert len(result)==3
    assert result[0][1] == "Methionine / Tryptophan oxidation (severity: medium)"
    assert result[1][1] == "Deamidation (severity: low)"
    assert result[2][1] == "Deamidation (severity: high)"

    corrected_seq = list(generic_test_seq)
    corrected_seq[102:104] = "AL"
    corrected_seq = "".join(corrected_seq)
    alignment = annotation_tool.analyze_seq(corrected_seq)
    assert len(search_tool.analyze_seq(corrected_seq,
        alignment, "imgt", "imgt")) == 2

    corrected_seq = list(corrected_seq)
    corrected_seq[31] = "A"
    corrected_seq = "".join(corrected_seq)
    alignment = annotation_tool.analyze_seq(corrected_seq)
    assert len(search_tool.analyze_seq(corrected_seq,
        alignment, "imgt", "imgt")) == 1

    corrected_seq = list(corrected_seq)
    corrected_seq[56:58] = "AL"
    corrected_seq = "".join(corrected_seq)
    alignment = annotation_tool.analyze_seq(corrected_seq)
    assert len(search_tool.analyze_seq(corrected_seq,
        alignment, "imgt", "imgt")) == 0

    corrected_seq = list(corrected_seq)
    liability_seq = corrected_seq.copy()
    liability_seq[3] = "C"
    liability_seq = "".join(liability_seq)
    alignment = annotation_tool.analyze_seq(liability_seq)
    result = search_tool.analyze_seq(liability_seq,
            alignment, "imgt", "imgt")
    assert result[0][1] == "Unusual cysteine (severity: high)"

    liability_seq = corrected_seq.copy()
    liability_seq[-5:-2] = "NAS"
    liability_seq = "".join(liability_seq)
    alignment = annotation_tool.analyze_seq(liability_seq)
    result = search_tool.analyze_seq(liability_seq,
            alignment, "imgt", "imgt")
    assert result[0][1] == "N-glycosylation (severity: high)"

    liability_seq = corrected_seq.copy()
    liability_seq[3:5] = "M"
    liability_seq = "".join(liability_seq)
    alignment = annotation_tool.analyze_seq(liability_seq)
    result = search_tool.analyze_seq(liability_seq,
            alignment, "imgt", "imgt")
    assert len(result) == 0


@pytest.fixture
def generic_test_seq():
    """Get a generic heavy chain that can be used for
    checking for liabilities we introduce."""
    return "VHLQQSGAELMKPGASVKISCKASGYTFITYWIEWVKQRPGHGLEWIGDILPGSGSTNYNENFKGKATFTADSSSNTAYMQLSSLTSEDSAVYYCARSGYYGNSGFAYWGQGTLVTVSA"
