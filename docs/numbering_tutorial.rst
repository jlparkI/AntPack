Antibody numbering in AntPack
===============================================

To start, use the ``SingleChainAnnotator`` or ``PairedChainAnnotator`` tools.
You can use these to number sequences, "trim" the resulting alignment and
convert a list of sequences into a fixed-length multiple sequence alignment
or MSA:::

  from antpack import SingleChainAnnotator, PairedChainAnnotator
  aligner = SingleChainAnnotator(chains=["H", "K", "L"], scheme="imgt")
  paired_aligner = PairedChainAnnotator(scheme = "imgt")

Both tools have the same basic methods / functions available.

For single chains: if you don't know what type of chain you're working with, leave
``chains`` as default and SingleChainAnnotator will figure out the chain
type. If you DO know that all of your chains are either
heavy ["H"] or light ["K", "L"] set SingleChainAnnotator to only
look for that chain.

``PairedChainAnnotator`` is designed to work with sequences that contain both a
light and a heavy chain (in any order) but can also handle single chains. Keep in
mind that PairedChainAnnotator *will* try to find two chains in the input sequence,
so your clue that there is only one chain present will be a very low percent
identity and/or an error message for one of the two chains. It is a little slower than
``SingleChainAnnotator`` because it has to do some additional operations.

Some prior versions accepted an option called ``compress_init_gaps``. This option
is deprecated as of v0.3.6.

.. autoclass:: antpack.SingleChainAnnotator
   :special-members: __init__
   :members: analyze_seq, analyze_seqs, assign_cdr_labels, sort_position_codes, build_msa, trim_alignment

.. autoclass:: antpack.PairedChainAnnotator
   :special-members: __init__
   :members: analyze_seq, analyze_seqs, assign_cdr_labels, sort_position_codes, build_msa, trim_alignment

Notice that these tools do not have a multithreading
option. That's because the *best* way to do multithreading / multiprocessing
for antibody numbering for a large number of sequences is to split the
antibody sequences up into batches, and you can do this easily by using
Python ``multiprocessing``. In each process, create a ``SingleChainAnnotator``
and use that to number one batch of the sequences.

Aho, IMGT, Martin ("modern Chothia") and Kabat are supported numbering schemes.

Notice that it's easy to convert the output of either annotator into a fixed-length
MSA by using ``build_msa``. If you have a large set of sequences that's too large to
work with in memory, however, ``build_msa`` obviously will not help. In this case, you can
instead loop over the sequences, number each of them and keep track of all the unique
position codes you've seen (e.g. by adding them to a set). The ``sort_position_codes``
function can then sort the resulting list of unique position codes to the correct
ordering for that scheme. You can then create a dictionary mapping position codes
to positions in a fixed length array, e.g.:::

  position_dict = {k:i for i, k in enumerate(sorted_position_codes)}

and then use this when writing sequences to file or e.g. one-hot encoding them in
an array. For examples of how to do this, see the numbering example on the
main page. If "-" is present, ``sort_position_codes`` will always remove it when
sorting; thus, if the list you pass to this function contains '-', that character
will be removed before sorting.

``assign_cdr_labels`` is useful to figure out which portions of a numbered sequence are
cdr or framework.
