Antibody numbering in AntPack
===============================================

To start, use the ``SingleChainAnnotator`` tool:::

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
   :members: analyze_seqs, analyze_seq, sort_position_codes


Notice that ``SingleChainAnnotator`` does not have a multithreading
option. That's because the *best* way to do multithreading / multiprocessing
for antibody numbering for a large number of sequences is to split the
antibody sequences up into batches, and you can do this easily by using
Python ``multiprocessing``. In each process, create a ``SingleChainAnnotator``
and use that to number one batch of the sequences.

Currently IMGT, Martin ("modern Chothia") and Kabat are supported numbering schemes.
We're interested in adding Aho and may do this soon.

The ``sort_position_codes`` function is useful if you're trying to merge many sequences
into an MSA or a fixed length array. It's easy to use AntPack to loop over the sequences
and number each of them while keeping track of all the unique position codes you've
seen so far, but of course at the end of this process, the set of unique position codes
you've extracted won't be in the correct order. ``sort_position_codes`` can convert a
list of unique position codes to the correct ordering for that scheme. It's then
easy to create a dictionary mapping position codes to positions in a fixed length
array, e.g.:::

  position_dict = {k:i for i, k in enumerate(sorted_position_codes)}

and then use this when writing sequences to file or e.g. one-hot encoding them in
an array. For an example of how to do this, see the numbering example on the
main page.


Numbering sequences that contain both heavy and light chain variable regions
=============================================================================

If convenient, you can extract the heavy and light chain variable regions from a sequence
that contains both using MultiChainAnnotator; they can then be
numbered as shown above. This tool just uses a SingleChainAnnotator to
find and extract chains from the input sequence. You can do the
same thing yourself by setting up two SingleChainAnnotators, one of which
has chains ``['K', 'L']`` and therefore only recognizes light chains and 
the other of which has chains ``['H']`` and therefore only recognizes
heavy. Sometimes however the MultiChainAnnotator may be more convenient
to use.

.. autoclass:: antpack.MultiChainAnnotator
   :special-members: __init__
   :members: analyze_seq

