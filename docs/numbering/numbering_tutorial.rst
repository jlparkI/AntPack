Antibody numbering in AntPack
===============================================

In general use the ``SingleChainAnnotator`` tool:::

  from antpack import SingleChainAnnotator
  aligner = SingleChainAnnotator(chains=["H", "K", "L"], scheme="imgt",
                        compress_init_gaps = False)

If you don't know what type of chain you're working with, leave
``chains`` as default and SingleChainAnnotator will figure out the chain
type. If you DO know that all of your chains are either
heavy ["H"] or light ["K", "L"] set SingleChainAnnotator to only
look for that chain.

``compress_init_gaps`` rearranges gaps in the first five positions of the
numbered sequence when there are small n-terminal deletions. Some other tools
like to have the gaps at the beginning of the numbering. You can mimic this
behavior by setting ``compress_init_gaps`` to True. We recommend leaving this
as False (default).

.. autoclass:: antpack.SingleChainAnnotator
   :special-members: __init__
   :members: analyze_seqs, analyze_seq


Notice that ``SingleChainAnnotator`` does not have a multithreading
option. That's because the *best* way to do multithreading / multiprocessing
for antibody numbering for a large number of sequences is to split the
antibody sequences up into batches, and you can do this easily by using
Python ``multiprocessing``. In each process, create a ``SingleChainAnnotator``
and use that to number one batch of the sequences.

Currently IMGT, Martin ("modern Chothia") and Kabat are supported numbering schemes.
We're interested in adding Aho and may do this soon.


Numbering sequences that contain both heavy and light chain variable regions
=============================================================================

If convenient, you can extract the heavy and light chain variable regions from a sequence
that contains both using MultiChainAnnotator; they can then be
numbered as shown above. This tool just uses a SingleChainAnnotator to
find and extract chains from the input sequence.

.. autoclass:: antpack.MultiChainAnnotator
   :special-members: __init__
   :members: analyze_seq

