from Bio import SeqIO
from antpack import PairedChainAnnotator


with open("input_fasta_data.fasta", "r") as fh:
    seqrecs = [s for s in SeqIO.parse(fh, "fasta")]
    seqs = [str(s.seq) for s in seqrecs]

pca = PairedChainAnnotator()

heavy_annot, light_annot = [], []
for seq in seqs:
    h_an, l_an = pca.analyze_seq(seq)
    heavy_annot.append(h_an)
    light_annot.append(l_an)

woah = pca.build_msa(seqs, heavy_annot)
