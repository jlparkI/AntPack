Finding possible sequence liabilities
======================================

Some motifs are known for propensity to oxidation, N-glycosylation,
isomerization, hydrolysis, or deamidation under certain conditions. These
kinds of liabilities can be problematic for manufacturing and stability.
AntPack contains a simple tool you can use to search for motifs that *may*
pose some risk of one or more of these reactions. We use the list of
known possible problematic motifs from `Satlawa et al. <https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1011881>`_,
but excluding low risk motifs (e.g. [STK]N which is low-risk for deamidation).

Note that this type of search is *prone* to false-positives. A motif that
can in principle be N-glycosylated or undergo pH-dependent hydrolysis will
not *always* undergo these reactions; it depends on the context. You can think
of this as casting a 'wide net' to find possible liabilities which can then
be narrowed down to find those which should be removed.

.. autoclass:: antpack.LiabilitySearchTool
   :special-members: __init__
   :members: analyze_seq
