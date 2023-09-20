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
in v0.0.2 -- but will be in the next version**.

Here's how to use ``SingleChainAnnotator``:::

  from antpack import SingleChainAnnotator
  aligner = SingleChainAnnotator(species=["all"], chains=["H", "K", "L"], scheme = "imgt")


You can then use the ``analyze_online_seqs`` (for a list of sequences), ``analyze_fasta``
(a generator for analyzing sequences in a fasta file) or ``analyze_seq`` (for analyzing
a single sequence). Details on these methods are below.

.. autoclass:: antpack.SingleChainAnnotator
   :special-members: __init__
   :members: analyze_online_seqs, analyze_fasta, analyze_seq


Notice that ``SingleChainAnnotator`` does not have a multithreading
option. That's because the *best* way to do multithreading / multiprocessing
for antibody numbering for a large number of sequences is to split the
antibody sequences up into batches, and you can do this easily by using
Python ``multiprocessing``. In each process, create a ``SingleChainAnnotator``
and use that to number one batch of the sequences.

Currently IMGT, Martin ("modern Chothia") and Kabat are supported numbering schemes.
We don't support original Chothia because (weirdly) it's defined inconsistently in
the literature, and Martin is basically an improved version of Chothia. IMGT and Kabat
are both fairly popular. We may add an updated version of IMGT (Aho) in a future version.
