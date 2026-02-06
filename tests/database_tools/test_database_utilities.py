"""Test database utility functions that the end user or
Python functions may need to access."""
from antpack import (number_imgt_imgt_cdr,
        SingleChainAnnotator, VJGeneTool)
from antpack.antpack_cpp_ext import vjgene_parser
from .database_utilities import get_vgene_code




def test_mab_renumbering(load_mab_nmbr_test_data):
    """Run a batch of test data (approximately 1600 sequences from the
    PDB) to ensure that numbering isolated cdrs is consistent with
    numbering the full sequence and extracting the cdr."""
    seqs, numberings = load_mab_nmbr_test_data
    numbering = numberings["imgt"]

    aligner = SingleChainAnnotator(chains=["H", "K", "L"],
        scheme="imgt")
    for seq, nmbr in zip(seqs, numbering):
        # Chain here is irrelevant since we are using IMGT.
        labels = aligner.assign_cdr_labels(nmbr, "H", "imgt")
        for i in range(1, 4):
            extracted_cdr, extracted_tokens = [], []
            for s, c, n in zip(seq, labels, nmbr):
                if c == f"cdr{i}":
                    extracted_cdr.append(s)
                    extracted_tokens.append(n)

            assigned_tokens = number_imgt_imgt_cdr(
                    ''.join(extracted_cdr), i)
            assert assigned_tokens == extracted_tokens



def test_tcr_renumbering(load_tcr_nmbr_test_data):
    """Run a batch of test data (approximately 1600 sequences from the
    PDB) to ensure that numbering isolated cdrs is consistent with
    numbering the full sequence and extracting the cdr."""
    seqs, numbering = load_tcr_nmbr_test_data

    aligner = SingleChainAnnotator(chains=["H", "K", "L"],
        scheme="imgt")
    for seq, nmbr in zip(seqs, numbering):
        # Chain here is irrelevant since we are using IMGT.
        labels = aligner.assign_cdr_labels(nmbr, "B", "imgt")
        for i in range(1, 4):
            extracted_cdr, extracted_tokens = [], []
            for s, c, n in zip(seq, labels, nmbr):
                if c == f"cdr{i}":
                    extracted_cdr.append(s)
                    extracted_tokens.append(n)

            assigned_tokens = number_imgt_imgt_cdr(
                    ''.join(extracted_cdr), i)
            assert assigned_tokens == extracted_tokens




def test_vjgene_parsing():
    """Test that vjgene parsing returns correct and
    expected data for a long list of v/j genes."""
    vjt = VJGeneTool()
    seq_lists = vjt.get_seq_lists()[1]

    for species, name_list in seq_lists.items():
        for name in name_list:
            res = vjgene_parser(name, species.split("_")[0])
            gt = get_vgene_code(name, species.split("_")[0])
            assert res == gt
            reassembled_gene = reassemble_gene(res,
                        name[2], name[3], name[:2])
            assert reassembled_gene == name



def reassemble_gene(gene_extracts, chain_type,
    gene_type, locus_type):
    """Reassembles the vgene from the extracted portions.
    Used for testing only, to confirm that the extracted
    information correctly maps back to the original gene."""
    output_gene = locus_type + chain_type + gene_type + \
            gene_extracts[1]
    if not gene_extracts[1].endswith("S") and \
            len(gene_extracts[2]) > 0:
        output_gene += '-'

    if len(gene_extracts[2]) > 0:
        output_gene += gene_extracts[2]

    if len(gene_extracts[3]) > 0:
        output_gene += "*" + gene_extracts[3]

    return output_gene
