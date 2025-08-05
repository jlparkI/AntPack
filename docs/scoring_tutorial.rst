Scoring sequences for humanness in AntPack
===============================================

Here's how to create a tool to use for scoring sequences:::

  from antpack import SequenceScoringTool
  scoring_tool = SequenceScoringTool(offer_classifier_option = False,
                    normalization = "none")


Note that the scoring tool can work with both heavy and light chains and
will by default figure out what type of chain you're giving it
for any sequence you provide.

The scoring tool can normalize scores; this is mostly useful if you
want to combine heavy and light chain scores into a single antibody
score, or if you're trying to compare scores assigned to specific
regions (frameworks or CDRs). There are two normalization options.
Setting ``normalization='training_set_adjust'`` subtracts the median
training set score for heavy and light chains from scores assigned to
those chains. This is useful for combining the heavy and light chain
scores for an antibody if you haven't masked out any regions or used
any masking options. If masking regions, ``normalization='normalize'``
is probably most useful; this divides the score by the number of non-
masked residues. ``'none'`` is the default.

The scoring tool can run in classifier mode if ``offer_classifier_option``
is True. (See the Background section for more on how this works and why
we don't generally recommend it.) Setting ``offer_classifier_option``
to False makes the tool a little more lightweight, since it doesn't
need to load some additional model parameters.

Scoring a sequence
--------------------

To score sequences with various options, use ``score_seqs``.
If you want more detailed information about why a sequence
If you want a mask for all positions EXCEPT a standard IMGT-defined
framework or cdr region, use ``get_standard_mask``; the resulting
mask can then be passed to ``score_seqs``. This will enable you
to just get the score for a specific region.

Finally and most importantly, you can use ``get_closest_clusters``
to find out which cluster(s) are closest to your sequence then
call ``calc_per_aa_probs`` to determine the probability of each
AA in your sequence given those clusters. This enables easy identification
of low-probability residues and motifs that could pose developability
issues and is illustrated in more detail in the examples.

See below for details.

.. autoclass:: antpack.SequenceScoringTool
   :members: score_seqs, get_diagnostic_info, get_standard_mask, get_standard_positions, convert_sequence_to_array, retrieve_cluster


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

To automatically humanize a sequence, start by creating a SequenceScoringTool
as shown above. (This is a significant change in AntPack v0.4; note that prior
versions had a separate HumanizationTool). Next, feed it the sequence you'd like
to humanize using the following function:

.. autoclass:: antpack.SequenceScoringTool
   :members: suggest_humanizing_mutations


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
