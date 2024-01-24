How antibody numbering in AntPack works
===============================================

Antibody numbering is essentially a form of multiple sequence alignment;
it inserts gaps into an input antibody heavy chain or light chain sequence
or into a template so that each AA in the input is assigned to a position
in an existing multiple sequence alignment. This process converts
antibody sequences into fixed-length vectors, so that subsequent
statistical analysis / machine learning is straightforward.
Additionally (and importantly), it assists in identifying CDRs (known
regions of high variability) and framework regions (known regions
of relatively lower variability).

Confusingly, there are a number of different numbering schemes for antibodies.
These are mostly very similar, but they tend to differ in how positions in / 
close to CDRs are assigned and in how CDRs are defined. The most popular
schemes include IMGT and Kabat. AntPack currently supports IMGT, Martin ("modern 
Chothia") and Kabat. For some background on how these
schemes differ and why you might prefer one or the other, we suggest reading
`this paper <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6198058/>`_.

Several different procedures have been suggested for antibody numbering.
The AbNum tool looks for specific subsequences typical in antibody heavy and
light chains and expands the numbering out from there. The ANARCI tool
uses HMMer to align the input sequences to a profile of germline sequences.
Finally, the AbRSA tool uses a custom global alignment procedure to
align the input sequence to a consensus sequence. (There are other
proprietary tools but we do not address these here.)

All of these tools have significant drawbacks. AbNum is only available as
a webserver, which means that sequences need to be uploaded, so that use
is very slow. ANARCI is available locally but is very slow for three reasons.
First, it aligns each input sequence against many different possible profiles,
even when those possible profiles are not relevant. You may know, for example,
that your input sequence is not a T-cell receptor. ANARCI will try to align it against
a variety of T-cell receptor profiles anyway. Second, ANARCI loads all input
sequences from a fasta file into memory before processing, then writes them
back to disk to query HMMer, which is highly inefficient if there are many
sequences. Third, ANARCI has to make a series of "tweaks" and corrections to
the HMMer alignment. If the input is a single domain sequence (e.g. a heavy
chain), ANARCI will align part of it using HMMer *a second time* to correct possible
issues with the first alignment.

AbRSA is faster than ANARCI or AbNum but is not open source. The available
AbRSA software has very few options available to users -- it will for
example print each result to the terminal and there does not seem (as of
this writing) to be any way to disable this, which slows down processing
dramatically. It's also harder to build into a Python-based data processing
pipeline, because it's a command-line tool with no Python wrapper.

In initial testing on about 1600 chains from PDB (more extensive benchmarking
available soon), AntPack is > 50x faster than ANARCI and > 25x faster
than AbRSA, taking 0.5 seconds to align the same number of sequences that
take > 35 seconds on ANARCI. If you know what chain type you're dealing with
(e.g. heavy or light), it's > 100x faster than ANARCI.

AntPack is written in Python-wrapped C++ and uses a very similar approach
to AbRSA with a custom global alignment. It uses position-specific scoring
to ensure that known highly conserved positions (e.g. the two cysteines)
are maintained and that gaps are inserted at desired places. The scoring
is constructed in such a way that most of the manual "tweaking" performed
by some other tools is unnecessary. It's also easy to build into a Python-based
workflow: create an e.g. ``SingleChainAnnotator`` class and use it to
annotate any sequences or fasta files you like.
