How antibody / TCR numbering in AntPack works
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
schemes include IMGT and Kabat. AntPack currently supports Aho, IMGT, Martin ("modern 
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
is very slow. ANARCI is available locally but is very slow. AbRSA is faster
than ANARCI or AbNum but is not open source and is still relatively slow.

In initial testing on about 1600 chains from PDB (see the links at the
github repo for more), AntPack is > 100x faster than ANARCI and > 25x faster
than AbRSA, taking < 0.2 seconds to align the same number of sequences that
take > 35 seconds on ANARCI.

TCRs are numbered using a different algorithm than the one used for antibodies
which is somewhat slower, although it is still faster by orders of magnitude
than ANARCI. There is plenty of room to further optimize the TCR alignment
algorithm and we may be able to improve this considerably in future versions.

``SingleChainAnnotator`` and ``PairedChainAnnotator`` contain functionality
for numbering sequences that may contain a single chain (heavy or light) or
up to two paired chains (heavy and light). They also contain functionality for
easily converting a list of numbered sequences into a fixed-length MSA
which is convenient for many kinds of analysis.
