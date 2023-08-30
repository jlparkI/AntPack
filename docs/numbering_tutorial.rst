Antibody numbering in AntPack
===============================================

There are two key tools for antibody numbering in
AntPack: the ``SingleChainAnnotator`` and ``MultiChainAnnotator``.
The former is most useful if you think the input sequence is
a single-domain sequence, e.g. a heavy chain. It can be used
to align a sequence containing both a heavy and light chain,
but will only return numbering for one or the other. Similarly,
``MultiChainAnnotator`` numbers sequences that likely contain
both heavy and light chains. It can be used on a single domain
sequence without any issues.

**NOTE: MultiChainAnnotator is not yet implemented / not available
in v0.0.1 -- but will be in the next version**.

Here's SingleChainAnnotator:

.. autoclass:: antpack.SingleChainAnnotator
   :special-members: __init__
   :members: analyze_online_seqs, analyze_fasta


Notice that ``SingleChainAnnotator`` does not have a multithreading
option. That's because the *best* way to do multithreading / multiprocessing
for antibody numbering for a large number of sequences is to split the
antibody sequences up into batches, and you can do this easily by using
Python ``multiprocessing``. In each process, create a ``SingleChainAnnotator``
and use that to number one batch of the sequences.

As noted above, in the first release only IMGT numbering is supported -- more
options coming soon.
