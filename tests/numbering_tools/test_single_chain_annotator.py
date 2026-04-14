"""Tests basic functionality for the SingleChainAnnotator class."""
import random
import pytest
from antpack import SingleChainAnnotator

from ..conftest import (standard_aa_list,
                    load_mab_nmbr_test_data, load_tcr_nmbr_test_data)


def test_sca_error_checking():
    """Ensure the single chain annotator returns errors
    when unacceptable parameters are supplied."""
    aligner = SingleChainAnnotator(chains=["H", "K", "L"])
    for invalid_seq in ["YaY", "YBW", "Y K", "Y-K", "yAy"]:
        results = aligner.analyze_seqs([invalid_seq])
        assert results[0][3].startswith("Sequence contains invalid")

    # Repeat these tests with a SingleChainAnnotator for TCRs.
    aligner = SingleChainAnnotator(chains=["A", "B", "D", "G"])
    for invalid_seq in ["YaY", "YBW", "Y K", "Y-K", "yAy"]:
        results = aligner.analyze_seqs([invalid_seq])
        assert results[0][3].startswith("Sequence contains invalid")

    # It should be impossible to create an annotator for both TCRs and
    # mAbs -- let's check.
    with pytest.raises(RuntimeError):
        aligner = SingleChainAnnotator(chains=["A", "B", "D", "H"])
    with pytest.raises(RuntimeError):
        aligner = SingleChainAnnotator(chains=["H", "K", "A"])

    # The only scheme allowed for TCRs should be IMGT.
    with pytest.raises(RuntimeError):
        aligner = SingleChainAnnotator(chains=["A", "B", "D", "G"], scheme="martin")



@pytest.mark.parametrize("scheme", ["martin", "imgt", "kabat", "aho"])
def test_mab_chain_recognition(known_chain_mab_seqs, scheme):
    """Ensure that single chain annotator can correctly recognize the
    input chain when supplied with something that could be L or H,
    and ensure it can correctly detect sequences with large deletions
    that remove one or more conserved residues."""
    known_h, known_l, known_k = known_chain_mab_seqs
    aligner = SingleChainAnnotator(chains=["H", "K", "L"],
                scheme = scheme)
    results = aligner.analyze_seqs([known_k, known_l, known_h])
    assert results[0][2] == "K"
    assert results[1][2] == "L"
    assert results[2][2] == "H"

    assert aligner.analyze_seq(known_k)[2] == "K"
    assert aligner.analyze_seq(known_l)[2] == "L"
    assert aligner.analyze_seq(known_h)[2] == "H"

    bad_chain = known_h[:100]
    assert aligner.analyze_seqs([bad_chain])[0][3].startswith("Unexpected")
    assert aligner.analyze_seq(bad_chain)[3].startswith("Unexpected")


def test_tcr_chain_recognition(known_chain_tcr_seqs):
    """Ensure that single chain annotator can correctly recognize the
    input chain when supplied with something that could be A, B, D, G."""
    known_a, known_b, known_d, known_g = known_chain_tcr_seqs
    aligner = SingleChainAnnotator(chains=["A", "B", "D", "G"])
    results = aligner.analyze_seqs([known_a, known_b, known_d, known_g])
    assert results[0][2] == "A"
    assert results[1][2] == "B"
    assert results[2][2] == "D"
    assert results[3][2] == "G"



@pytest.mark.parametrize("scheme", ["martin", "imgt", "kabat", "aho"])
def test_mab_performance(load_mab_nmbr_test_data, scheme):
    """Run a batch of test data (approximately 1600 sequences from the
    PDB) to ensure that numbering is consistent with numbering generated
    by another tool. There will occasionally be small differences in
    cases where there are multiple possible acceptable alignments,
    but in general we expect the numbering to be the same the vast
    majority of the time."""
    seqs, nmbr = load_mab_nmbr_test_data
    numbering = nmbr[scheme]

    aligner = SingleChainAnnotator(chains=["H", "K", "L"],
        scheme=scheme)
    aligner_results = aligner.analyze_seqs(seqs)
    total_comparisons, num_correct = compare_results(aligner_results,
                        numbering, seqs)
    print(f"{scheme}: Total comparisons: {total_comparisons}. "
            f"Num matching: {num_correct}.", flush=True)
    assert num_correct / total_comparisons > 0.97



def test_msa_construction(load_mab_nmbr_test_data):
    """Check that sequences can be correctly formed into an MSA."""
    # The last one produced is IMGT. Use this to test MSA construction.
    seqs, _ = load_mab_nmbr_test_data
    aligner = SingleChainAnnotator(chains=["H", "K", "L"],
        scheme="imgt")
    aligner_results = aligner.analyze_seqs(seqs)

    hnumbering, lnumbering = [], []
    hseqs, lseqs = [], []

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
            assert len(msa_seq) == len(position_set)


def test_tcr_performance(load_tcr_nmbr_test_data):
    """Run a batch of test data (approximately 525 sequences from the
    PDB) to ensure that numbering is consistent with numbering generated
    by another tool. There will occasionally be small differences in
    cases where there are multiple possible acceptable alignments,
    but in general we expect the numbering to be the same the vast
    majority of the time. Also check that the sequences can be correctly
    formed into an MSA. This test is for TCRs specifically, mabs are
    tested separately since they use a different alignment workflow."""
    seqs, imgt_num = load_tcr_nmbr_test_data

    aligner = SingleChainAnnotator(chains=["A", "B", "D", "G"])
    aligner_results = aligner.analyze_seqs(seqs)
    total_comparisons, num_correct = compare_results(aligner_results,
                                imgt_num, seqs)
    print(f"*TCRS, IMGT: Total comparisons: {total_comparisons}. "
            f"Num matching: {num_correct}.", flush=True)
    assert num_correct / total_comparisons > 0.97

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
            assert len(msa_seq) == len(position_set)



def test_extended_length_handling(load_mab_nmbr_test_data,
            standard_aa_list):
    """Check situations where the sequence is much longer than
    a typical variable chain. For now this test is only for
    mAbs."""
    seqs, _ = load_mab_nmbr_test_data
    std_aas = standard_aa_list[:20]
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

        muddled_seq = seq + "YYYGGGGGSSSS" + "".join([random.choice(std_aas) for i in range(250)])
        alt_alignment = aligner.analyze_seq(muddled_seq)
        assert len(alt_alignment[0]) == len(muddled_seq)


        trimmed_alt = aligner.trim_alignment(muddled_seq, alt_alignment)
        trimmed_gt = aligner.trim_alignment(seq, alignment)
        assert trimmed_alt[0] == trimmed_gt[0]



def test_x_handling(load_mab_nmbr_test_data):
    """Check situations where one or more letters has
    been replaced with X. This is specifically for mAbs."""
    seqs, _ = load_mab_nmbr_test_data
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

        assert alignment[0] == pre_x_alignment[0]



def compare_results(results, comparator_numbering, seqs):
    """Compares the numbering generated by AntPack with a comparator,
    and returns the number correct vs total comparisons. Also contains
    some (commented-out) code for writing non-matching results to
    a temporary file for closer inspection."""
    total_comparisons, num_correct = 0, 0
    #outhandle = open(f"temp.csv", "w+", encoding="utf-8")
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


@pytest.fixture
def known_chain_mab_seqs():
    """Returns sequences of known chain type in the order
    H, K, L."""
    known_k = ("DIVMTQSPSSLTVTAGEKVTMSCKSSQSLLSSGNQKNYLTWYQQIPGQPPKLLIYWASTR"
                "ESGVPDRFTGSGSGTDFTLTINSVQAEDLAVYYCQNDYTYPLTFGAGTKLELKRTV")
    known_l = ("QSALTQPASVSGSPGQSITISCTGTTSDVGTYNFVSWYQQHPGKAPKAIIFDVTNRPSGI"
                "SNRFSGSKFGNTASLTISGLQAEDEADYYCAAYTVASTLLFGGGTKVTVLRQP")
    known_h = ("VHLQQSGAELMKPGASVKISCKASGYTFITYWIEWVKQRPGHGLEWIGDILPGSGSTNYN"
                "ENFKGKATFTADSSSNTAYMQLSSLTSEDSAVYYCARSGYYGNSGFAYWGQGTLVTVSA")
    return known_h, known_l, known_k


@pytest.fixture
def known_chain_tcr_seqs():
    """Returns sequences of known chain type in the order
    A,B,D,G."""
    known_a = ("MQQVRQSPQSLTVWEGETAILNCSYENSAFDYFPWYQQFPGEGPALLIAIRSVSDKK"
            "EDGRFTIFFNKREKKLSLHITDSQPGDSATYFCAASATGANTGKLTFGHGTILRVHP")
    known_b = ("DAGVIQSPRHEVTEMGQEVTLRCKPISGHNSLFWYRQTMMRGLELLIYFNNNVPIDD"
            "SGMPEDRFSAKMPNASFSTLKIQPSEPRDSAVYFCASTWGRASTDTQYFGPGTRLTVL")
    known_d = ("AQKVTQAQSSVSMPVRKAVTLNCLYETSWWSYYIFWYKQLPSKEMIFLIRQGSDE"
            "QNAKSGRYSVNFKKAAKSVALTISALQLEDSAKYFCALGDPGGLNTDKLIFGKGTRVTVEP")
    known_g = ("SSNLEGGTKSVTRPTRSSAEITCDLTVINAFYIHWYLHQEGKAPQRLLYYDVSNSKDVLE"
            "SGLSPGKYYTHTPRRWSWILILRNLIENDSGVYYCATWDRGNPKTHYYKKLFGSGTTLVVTD")
    return known_a, known_b, known_d, known_g
