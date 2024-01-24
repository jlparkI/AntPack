Scoring sequences for humanness in AntPack
===============================================

Here's how to create a tool to use for scoring sequences:::

  from antpack import SequenceScoringTool
  scoring_tool = SequenceScoringTool(adjusted_scores = True,
                offer_classifier_option = False)


Note that the scoring tool can work with both heavy and light chains and
will by default figure out what type of chain you're giving it
for any sequence you provide.

The scoring tool adjusts scores by default by subtracting the median
training set score for heavy and light chains from scores assigned to
those chains. If you don't like this, you can turn this off by
setting ``adjusted_scores`` to False.

The scoring tool can run in classifier mode if ``offer_classifier_option``
is True. (See the Background section for more on how this works and why
we don't generally recommend it.) Setting ``offer_classifier_option``
to False makes the tool a little more lightweight, since it doesn't
need to load some additional model parameters.

Scoring a sequence
--------------------

To score a single sequence and retain lots of flexibility, use ``score_seq``.

.. autoclass:: antpack.SequenceScoringTool
   :members: score_seq

If scoring lots of sequences, ``score_seq`` is slow. Instead, use
``batch_score_seqs``.

.. autoclass:: antpack.SequenceScoringTool
   :members: batch_score_seqs
