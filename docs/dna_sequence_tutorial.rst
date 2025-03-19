Handling DNA sequences
==========================

All of the existing antibody numbering schemes were designed
to work with amino acid sequences. Consequently, while it's
possible to align germline genes to DNA sequences, this is
slower and for *most* (though not all) workflows does not
offer a significant advantage.

AntPack offers some tools to facilitate translating DNA
sequences into amino acids using the ``DNASeqTranslator``.
If you know the reading frame and whether to use forward or
reverse complement for a given sequence, the ``DNASeqTranslator``
can translate it for you (which is fairly trivial). More importantly,
however, the ``DNASeqTranslator`` can quickly figure out which
reading frame is correct and whether to use forward or reverse
complement if you don't already know. It does this by checking
all three possible reading frames (and doing the same for the
reverse complement if you indicate that the reverse complement
may contain the sequence) for kmers which are common in
antibody / TCR sequences. This check can be done very quickly and
thus makes it easy to determine the correct reading frame
without doing any expensive alignments.

Currently DNASeqTranslator supports DNA sequences consisting of
uppercase A, C, T, G or N (any codon containing N is translated
to X, which is an allowed letter for the AntPack numbering tools).
It does not accept sequences containing gaps.

**IMPORTANT:** ``DNASeqTranslator`` may return amino acid sequences
which contain stop codons, which are not allowed inputs for the
AntPack numbering / humanization tools. It's a good idea then to
check the sequences that are returned for stop codons and if
found decide how you want to handle these (depending on the kind
of data you're working with and your application).

For more on how to use ``DNASeqTranslator``, see below.

.. autoclass:: antpack.DNASeqTranslator
   :special-members: __init__
   :members: translate_dna_known_rf, translate_dna_unknown_rf
