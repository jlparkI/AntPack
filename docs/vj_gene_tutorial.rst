Assigning V and J genes
========================

AntPack can find the most similar V and J genes by looking
for the V and J amino acid sequences from a species of
interest that have the highest percent identity. Note that
there are some other tools that can generate a
probability for each possible recombination scenario. For
a detailed analysis, especially for repertoire data, these
tools may be better. If you just want to get the most similar
V and J genes to your amino acid sequence, AntPack can do
this easily.

To do this, use the VJGeneTool. The tool can tell you the
name of the most similar V and J genes for an input sequence.
It can determine similarity using either percent identity or
e-value (using the assigned numbering as the alignment).
You'll need to number the sequence first which could be done
with another tool (but is most easily done using AntPack).
You can next if desired retrieve the sequence of the assigned
vj genes.

Many tools try to assign a *single* V-gene and J-gene. In general
this is not necessarily correct -- it is not uncommon to find
situations where more than one V-gene and J-gene have the same
percent identity or very similar e-values. There are also
some germline genes that have different DNA sequences but the
same AA sequence. In these cases, AntPack
returns a list of the v-genes and j-genes that achieved the same
(or essentially equivalent) score, delimited or separated by the
character "_".

AntPack can also use either the "imgt" database or the "ogrdb"
database. In general we prefer "imgt". "ogrdb" assigns different
names to sequences with the same amino acid sequence, so in some
cases there is more than one name for the same gene; when this
happens, these are all reported separated by spaces, i.e.
"gene1 gene2". Still, the OGRDB database is available if useful
for your particular project. Where IMGT is concerned, note that some
sequences found in IMGT VQuest are excluded (pseudogenes,
partial sequences, sequences that do not have 'F' in the functionality
section etc.), so not all V and J genes in the IMGT db are in
AntPack. You can find out when AntPack's db was last updated
using the VJGeneTool (see below).

Also see the example, which illustrates how to use these capabilities.

.. autoclass:: antpack.VJGeneTool
   :special-members: __init__
   :members: get_vj_gene_sequence, assign_vj_genes, retrieve_db_dates

