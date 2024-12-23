seq = 'QMQLQESGPGLVKPSETLSLTCSVSGASISDSYWSWIRRSPGKGLEWIGYVHKSGDTNYSPSLKSRVNLSLDTSKNQVSLSLVAATAADSGKYYCARTLHGRRIYGIVAFNEWFTYFYMDVWGNGTQVTVSS'

from antpack import VJGeneTool, SingleChainAnnotator

sch = "kabat"
#imgt_tool = VJGeneTool(scheme="imgt")
#imgt_aligner = SingleChainAnnotator(scheme="imgt")

vj = VJGeneTool(scheme=sch)
alternate_aligner = SingleChainAnnotator(scheme=sch)

annotation2 = alternate_aligner.analyze_seq(seq)

pred_vgene2, pred_jgene2, pidv2, pidj2 = vj.assign_vj_genes(
                            annotation2, seq, "human", "identity")
import pdb
pdb.set_trace()

#trualign = imgt_aligner.analyze_seq(seq)
#tru_v, tru_j, tru_pidv, tru_pidj = imgt_tool.assign_vj_genes(trualign,
#        seq, "human", "identity")
