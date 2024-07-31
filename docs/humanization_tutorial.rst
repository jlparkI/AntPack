Humanizing sequences with AntPack
===============================================

AntPack can also suggest mutations to "humanize" a sequence or make
it more human. You can manually humanize a sequence
by 1) scoring it using AntPack, 2) retrieving the closest clusters
(see the "Generating new human sequences" page), 3) determining
which regions of the sequence are least human, 4) mutating these
and 5) rescoring the sequence. This can be however a little tedious.
AntPack offers an easy way to automatically
choose mutations that lie at a selected location along that tradeoff
curve.

Regardless of which approach you take, we suggest you carefully review
suggested mutations and re-score the altered sequence to be sure that
the changes that have been made are compatible with your objectives.

To automatically humanize a sequence, start by creating a humanization
tool:::

  from antpack import HumanizationTool
  h_tool = HumanizationTool()


Next, feed it the sequence you'd like to humanize using the following
function:

.. autoclass:: antpack.HumanizationTool
   :members: suggest_mutations


The "knob" you can turn here is ``s_thresh``. A value <= 1 means that
AntPack will basically do a straight *in silico* CDR graft. This is
straightforward but may lose affinity. A value > 1 means that AntPack will
1) score the straight *in silico* CDR graft of your sequence, which will
generally achieve a very high humanness score (the "optimal score"), then
2) revert as many positions as possible to the original sequence without
causing the score to go over s_thresh * optimal score. Larger values of s_thresh
cause more of the original sequence to be preserved at the expense of
a smaller improvement in humanness.

Kabat-defined CDRs are excluded from humanization. If you want to exclude
additional key positions, you can do this by passing a list of IMGT-numbered
positions you would like to be excluded from consideration.
