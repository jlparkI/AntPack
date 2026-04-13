Tools for search and clustering
===============================================

Starting in v0.5, AntPack has tools for quickly
building antibody databases from millions to hundreds
of millions (or even billions) of sequences. AntPack
constructs a special indexing that enables it to quickly
search these databases for sequences that are similar
to a query sequence (they have CDRs that differ by up
to some percent identity, may have the same vgene or
vgene family etc.)

Search speeds are faster than
MMSeqs or KA-Search by orders of magnitude; AntPack
can run thousands of searches against a hundred million
sequence database in a few seconds. In fact, search
speed is fast enough that (depending on hardware)
AntPack is able to run single linkage clustering
(aka clonotyping) on 10-40 million sequence databases
in a few hours. Single linkage clustering has :math:O(N^2) .,
so this procedure is unlikely to scale to hundred million
or billion sequence datasets (at least not very
efficiently)! (Upcoming versions of AntPack will include
even more scalable alternatives for clustering large
datasets).

If you want to cluster a small dataset (e.g. low thousands),
AntPack has some very flexible tools to quickly construct a
distance matrix; this in turn can be used by a wide variety
of clustering algorithms in scikit-learn. Finally, if you
want to cluster framework regions (which are less variable
than CDRs), AntPack also includes a mixture-model based tool
that provides :math:O(N) . clustering. It's not as high-resolution
as single-linkage but its excellent scalability makes it
useful for some specific situations.

Overall, we recommend:

* If interested in creating a large searchable database,
  use the database construction and search tools (see
  "Database search and clustering tutorial" on the main
  page).

* If interested in clustering by similarity in CDR regions
  and shared vgene assignment for datasets < 40 million
  sequences or so, use the database construction and
  search tools (see "Database search and clustering tutorial"
  on the main page).

* If clustering a small (low thousands or less) dataset,
  use AntPack's distance matrix construction tools, then choose
  an algorithm from scikit-learn or scipy. See 
  "Clustering for small datasets and frameworks" on the main
  page.

* If you need a highly scalable (but possible lower-resolution)
  clustering, or if you're interested in clustering less
  variable framework regions for one reason or another,
  consider using a mixture model. See
  "Clustering for small datasets and frameworks" on the main page.
