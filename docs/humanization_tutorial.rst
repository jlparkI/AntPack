Humanizing sequences with AntPack
===============================================

AntPack can also suggest mutations to "humanize" a sequence or make
it more human. Humanization always involves a tradeoff between
improving the humanness score as much as possible and preserving
as much of the original sequence as possible (to reduce the risk
of loss of affinity). You can of course manually humanize a sequence
by 1) scoring it using AntPack, 2) retrieving the closest cluster
(see the "Generating new human sequences" page), 3) determining
which regions of the sequence are least human, 4) mutating these
and 5) rescoring the sequence. This can be however a little tedious.
AntPack offers an easy way to automatically
choose mutations that lie at a selected location along that tradeoff
curve. 

To automatically humanize a sequence, start by creating a scoring
tool:::

  from antpack import SequenceScoringTool
  scoring_tool = SequenceScoringTool()


Next, feed it the sequence you'd like to humanize using the following
function:

.. autoclass:: antpack.SequenceScoringTool
   :members: suggest_mutations


The "knob" you can turn here is ``s_thresh``. A value <= 1 means that
AntPack will basically do a straight *in silico* CDR graft. This is
not usually desirable. A value > 1 means that AntPack will 1) score
the straight *in silico* CDR graft of your sequence, which will generally
achieve a very high humanness score (the "optimal score"), then 2) revert
as many positions as possible to the original sequence without causing
the score to go over s_thresh * optimal score. Larger values of s_thresh
cause more of the original sequence to be preserved at the expense of
a smaller improvement in humanness.

To illustrate what this tradeoff looks like, here are the results of using
AntPack to humanize 25 antibodies from Prihoda et al. These antibodies were
originally humanized using other techniques (the "Experimental" column) or
using a more traditional straight CDR graft. We can see how using AntPack / SAM
with a straight in silico graft, with s_thresh set to 1.10 (110 percent)
and with s_thresh set to 1.25 (125 percent) offers a tradeoff between
preservation of the original sequence and increased humanness score.
