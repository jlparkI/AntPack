Antibody numbering in AntPack
===============================================

Currently, to number sequences in AntPack, you can use
the ``SingleChainAnnotator`` tool. This tool finds the
region of your input sequence that corresponds most
closely to one of a selection of chains and numbers it.
In future versions, we'll add a second tool for numbering
extracting heavy and light chain regions from sequences
that likely contain both.

Here's how to use ``SingleChainAnnotator``:::

  from antpack import SingleChainAnnotator
  aligner = SingleChainAnnotator(chains=["H", "K", "L"], scheme = "imgt",
                        compress_init_gaps = False)


If you don't know what type of chain you're working with, leave
``chains`` as default and SingleChainAnnotator will figure out the chain
type for each input sequence. If you DO know that all of your chains are either
heavy ("H") or light (["K", "L"]) you can set SingleChainAnnotator to only
look for that chain type, which will speed things up slightly. ``scheme`` can
be one of ``martin``, ``imgt`` or ``kabat``.

``compress_init_gaps`` rearranges gaps in the first five positions of the
numbered sequence when there are small n-terminal deletions. Other tools
like to have the gaps at the beginning of the numbering (i.e. if a gap is at 
IMGT position 3, move it to IMGT position 1). You can mimic this
behavior by setting ``compress_init_gaps`` to True, or turn it off by
setting it to False, in which case the gaps will be left in what seem to be the
most sensible locations based on the alignment.

You can use the ``analyze_online_seqs`` (for a list of sequences), ``analyze_fasta``
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
We're interested in adding Aho and may do this soon.
