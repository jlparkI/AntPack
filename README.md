# AntPack

AntPack is a Python package / toolkit for antibody numbering, data processing,
statistical inference and machine learning for antibody sequences. For usage,
see [the docs](https://antpack.readthedocs.io/en/latest/index.html).
It is currently in active development, so we are adding new features
and making more improvements periodically.

## What's new in v0.3.8.5

v0.3.8 adds numbering (using the IMGT scheme) and VJ
gene assignment for TCRs. For both TCRs and mAbs, if
you don't know which species to check for germline
gene assignment, you can now check all of them by
passing "unknown" for species. For mAbs, germline
gene assignment for rabbits and alpacas (in addition
to humans and mice) is now supported. It also now supports
cross-scheme CDR assignment, i.e. you can number using
one scheme but assign CDR labels using the CDR definitions
from another (e.g. number using IMGT but assign CDR
labels and search for liabilities using Kabat CDR definitions).

There are various other minor improvements to the Python API, the CLI
and the GUI, which is useful for quickly viewing
a few sequences and comparing them to their assigned
VJ genes.

v0.3.8.5 contains a couple of minor bug fixes and error messages
that are more clear when a highly conserved residue is missing
(the message now more clearly indicates which group of highly
conserved residues was altered).

## Installation

The only required dependency is numpy. If you want to run the GUI,
however, there are two additional dependencies you'll need to install
(in future versions these dependencies may no longer be needed):
```
pip install pyside6 qt_material
```

If you don't plan to use the GUI, you don't need those dependencies.
Either way, to install antpack run:
```
pip install antpack
```

AntPack is distributed as a wheel precompiled for most platforms and CPython >= 3.8,
so installation should be very straightforward. A source distribution is also available
(C++17) in case there is any need to compile from source.

## Capabilities


#### Antibody numbering

Numbering antibody sequences is an important precursor for many statistical inference /
machine learning applications. AntPack is orders of magnitude faster for numbering
antibody sequences than existing tools in the literature (e.g. ANARCI, AbRSA),
while providing >= reliability. AntPack also provides tools for merging a list
of numbered sequences into an MSA, for easy extraction of specific CDRs and
framework regions, and for TCR numbering.


#### V / J genes

Identifying the most similar human/mouse/llama/rabbit V / J gene sequences is useful
for a variety of purposes. AntPack provides tools for determining which human
V and J gene sequences are most similar to the variable region chain provided
as input.


#### Humanness and developability

Minimizing the risk of immunogenicity is important for selecting clinical
candidates. In AntPack v0.1.0, we introduce a simple, fully interpretable
generative model for human heavy and light chains that outperforms all
comparators in the literature on a large held-out test set for distinguishing
human sequences from those of other species. This scoring tool can be used
to score sequences for humanness, suggest modifications to make them more
human, identify liabilities, and generate highly human sequences that contain
selected motifs.


#### Finding developability liabilities

Some sequence motifs are known to be associated with developability issues -- certain
motifs are known, for example, to be prone to N-glycosylation or deamidation. AntPack
provides a tool for finding these "liability" motifs in an input sequence. Note that
that identifying liabilities through finding motifs in this way is known to be prone
to false positives (an N-glycosylation motif, for example, will not always be glycosylated).
Still, these kinds of alerts can be useful for making yourself aware of potential
developability issues.


### Licensing

AntPack is licensed under a GPL license, which means that you are free to
use it for your own data analysis in any manner you like irrespective of
whether you work on academic research or industrial R&D / QC.

If you are writing software intended for sale & distribution, however,
any software you create that uses AntPack must also be open-source, must
use the GPL license and must acknowledge AntPack appropriately. If you
are interested in using AntPack in a closed-source software product,
please contact us to obtain a version of AntPack licensed for this use.


### Citing this work

If using AntPack in research intended for publication, please cite
either the preprint:

Jonathan Parkinson and Wei Wang. 2024. For antibody sequence generative modeling,
mixture models may be all you need.
[bioRxiv: https://doi.org/10.1101/2024.01.27.577555](https://www.biorxiv.org/content/10.1101/2024.01.27.577555v1)

or the final paper in [Bioinformatics.](https://academic.oup.com/bioinformatics/article/40/5/btae278/7656770)

#### Acknowledgements

The authors & maintainers gratefully acknowledge the contributions of
Japanese artist Yusuke Kamiyamane for the creation of the Fugue icon set.
