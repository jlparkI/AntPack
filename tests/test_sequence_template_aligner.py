"""Tests SequenceTemplateAligners to ensure that
new sequences can be added to existing msas in
a straightforward way, and that the aligner used
for database search is functioning correctly."""
import os
import random
import gzip
import unittest
from Bio import SeqIO
from antpack import SingleChainAnnotator as SCA
from antpack.antpack_cpp_ext import (SequenceTemplateAligner)


class TestSequenceTemplateAligner(unittest.TestCase):


    def test_error_checking(self):
        """Ensure that absurd inputs generate expected
        errors."""
        _, position_codes, _ = load_antibody_test_data()
        with self.assertRaises(RuntimeError):
            sta = SequenceTemplateAligner(position_codes,
                    "L", "")

        with self.assertRaises(RuntimeError):
            sta = SequenceTemplateAligner(position_codes,
                    "Z", "imgt")

        with self.assertRaises(RuntimeError):
            sta = SequenceTemplateAligner(position_codes,
                    "H", "imgt", "gollum")

        with self.assertRaises(RuntimeError):
            sta = SequenceTemplateAligner(position_codes + ['Z'],
                    "H", "imgt", "")


    def test_alignment_and_slicing(self):
        """Check to make sure the template aligner correctly performs
        alignment and slicing operations."""
        sc = SCA()

        for chain_type in ["H", "L"]:
            seqs, position_codes, region_labels = load_antibody_test_data(chain_type)
            unaligned_seqs = [s.replace('-', '') for s in seqs]
            unaligned_seq_nmbr = sc.analyze_seqs(unaligned_seqs)

            sta = SequenceTemplateAligner(position_codes,
                        chain_type, "imgt", "")
            realigned_seqs = [sta.align_sequence(s, u[0], False) for u, s in zip(
                unaligned_seq_nmbr, unaligned_seqs)]
            for r, s in zip(realigned_seqs, seqs):
                self.assertTrue(r==s)

            for region in ["all", "fmwk", "cdr2"]:
                realigned_seqs = [sta.align_and_slice_sequence(s, u[0], region, False)
                        for u, s in zip(unaligned_seq_nmbr, unaligned_seqs)]
                if region == "all":
                    sliced_seqs = seqs
                else:
                    sliced_seqs = [''.join([l for token, l in zip(region_labels, s)
                        if token.startswith(region)])
                            for s in seqs]
                for r, s in zip(realigned_seqs, sliced_seqs):
                    self.assertTrue(r==s)

                if region == "all":
                    expected_region_size = len(seqs[0])
                else:
                    expected_region_size = len([token for token in region_labels
                        if token.startswith(region)])
                self.assertTrue(expected_region_size==sta.get_region_size(region))

            self.assertTrue(sta.get_template_length() ==
                    len(position_codes))



    def test_forward_reverse_numbering(self):
        """Check to make sure that the forward and reverse numbering
        returned by the template aligner makes sense."""
        sc = SCA()

        for chain_type in ["H", "L"]:
            seqs, position_codes, _ = load_antibody_test_data(chain_type)
            unaligned_seqs = [s.replace('-', '') for s in seqs]
            unaligned_seq_nmbr = sc.analyze_seqs(unaligned_seqs)

            pcode_dict = {p:i for i, p in enumerate(position_codes)}
            sta = SequenceTemplateAligner(position_codes,
                        chain_type, "imgt", "")

            for unaligned_seq, nmbr in zip(
                    unaligned_seqs, unaligned_seq_nmbr):
                gt_fwd_numbering = [pcode_dict[p] for p in nmbr[0]]
                self.assertTrue(gt_fwd_numbering ==
                        sta.retrieve_alignment_forward_numbering(
                            unaligned_seq, nmbr[0], False))

                reverse_dict = {p:i for i, p in enumerate(nmbr[0])}
                gt_rev_numbering = []
                for p in position_codes:
                    if p in reverse_dict:
                        gt_rev_numbering.append(reverse_dict[p])
                    else:
                        gt_rev_numbering.append(-1)
                self.assertTrue(gt_rev_numbering ==
                        sta.retrieve_alignment_back_numbering(
                            unaligned_seq, nmbr[0], False))



def load_antibody_test_data(chain_type="L"):
    """Loads some saved test data specific for antibodies."""
    start_dir = os.path.abspath(os.path.dirname(__file__))

    data_path = os.path.join(start_dir, "test_data",
            "addtnl_test_data.fasta.gz")
    with gzip.open(data_path, "rt") as fhandle:
        seqrecs = list(SeqIO.parse(fhandle, "fasta"))

    sc = SCA()
    if chain_type == "L":
        light_seqs = [str(seqrec.seq) for seqrec in seqrecs if "light" in seqrec.name]
        light_alignments = sc.analyze_seqs(light_seqs)
        lpos_codes, lmsa = sc.build_msa(light_seqs, light_alignments)
        region_labels = sc.assign_cdr_labels(lpos_codes, chain_type)
        return lmsa, lpos_codes, region_labels


    heavy_seqs = [str(seqrec.seq) for seqrec in seqrecs if "heavy" in seqrec.name]
    heavy_alignments = sc.analyze_seqs(heavy_seqs)
    hpos_codes, hmsa = sc.build_msa(heavy_seqs, heavy_alignments)
    region_labels = sc.assign_cdr_labels(hpos_codes, chain_type)
    return hmsa, hpos_codes, region_labels



if __name__ == "__main__":
    unittest.main()
