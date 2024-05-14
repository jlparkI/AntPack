Assigning V and J genes
========================

AntPack can find the most similar V and J genes by looking
for the V and J amino acid sequences from a species of
interest that have the highest percent identity. Note that
in general this is not as reliable as using nucleotide
sequences, and that there are some other tools (e.g. IGOR)
that can generate a probability for each possible
recombination scenario. For a detailed analysis, especially
for repertoire data, these tools may be better. If you
just want to get the most similar V and J genes to your amino
acid sequence, AntPack should work fine.

To do this, use the VJGeneTool. The tool can tell you the
name of the most similar V and J genes for an input sequence.
You can either provide numbering (if you've already numbered
the sequence) or ask VJGeneTool to number it for you. You
can also retrieve the sequences of the assigned V and J
genes or the amino acid sequences of all genes in the same
V or J family (e.g. IGHV1).

Finally, note that AntPack uses its own internal database, which
is periodically updated based on the IMGT VQuest DB. Some
sequences found in IMGT VQuest are excluded (pseudogenes,
partial sequences, sequences that do not have 'F' in the functionality
section etc.), so not all V and J genes in the IMGT db are in
AntPack. You can find out when AntPack's db was last updated
using the VJGeneTool (see below).

Also see the example for an example of how to use these capabilities.

.. autoclass:: antpack.VJGeneTool
   :special-members: __init__
   :members: get_vj_gene_sequence, get_vj_gene_family, assign_numbered_sequence, assign_sequence, retrieve_db_dates

