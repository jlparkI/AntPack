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
may have been scored the way it was, use ``get_diagnostic_info``.
If you want a mask for all positions EXCEPT a standard IMGT-defined
framework or cdr region, use ``get_standard_mask``; the resulting
mask can then be passed to ``score_seqs``. This will enable you
to just get the score for a specific region.

See below for details.

.. autoclass:: antpack.SequenceScoringTool
   :members: score_seqs, get_diagnostic_info, get_standard_mask, get_standard_positions, convert_sequence_to_array, retrieve_cluster
