# AntPack

![Static Badge](https://img.shields.io/badge/coverage-80%25-green)

AntPack is a Python package / toolkit for numbering / alignment, VJ gene
assignment, developability prediction, database construction,
clustering, search, statistical inference and machine learning for
antibody sequences and TCRs.
It is currently in active development, so we are adding new features
and making more improvements periodically.

It is intended for use on Linux or Windows platforms.

## Usage

For usage for v0.3.9 and later, see the [docs at this link.](https://antpackdocumentationlatest.pages.dev/)

[You can find the docs here for older versions.](https://antpackdocumentation.pages.dev/)


## What's new in v0.5

v0.5 contains tools for constructing searchable antibody & TCR
databases. Once you've built a local database you can search it
for sequences similar to your query (percent identity above some
threshold, same vgene or same jgene etc.) at speeds hundreds of
times faster than what other tools (e.g. MMSeqs, KA-Search) can
achieve. AntPack now provides tools for clonotyping a database using
single-linkage clustering, which scales reasonably well to the 10
million to 50 million sequence range (depending on hardware). We
hope to introduce even more scalable tools for single linkage
clustering in future versions. There is also a mixture model
tool that provides lower-resolution clustering for even larger
datasets, and highly efficient tools for clustering small datasets.

The API for the SequenceScoringTool and humanization has been updated
and streamlined so if you are scoring sequences for humanness or
humanizing them please see the updated section of the docs.

## Licensing

Versions v0.3.9 and afterwards are licensed for academic and noncommercial
use only; you must first obtain a [free license key](https://antpackdocumentationlatest.pages.dev/installation)
to use. You do not have to re-setup the key when upgrading AntPack in
an existing venv or conda environment.

Versions prior to 0.3.9 were made available under the GPL, which means they
can be used it for your own data analysis in any manner you wish
whether you work in academia or industry,
but any software built using AntPack and intended for sale or distribution 
must also be open-source under the GPL license. Therefore if
you are currently using e.g. v0.3.8.6.2 for data analysis,
you should feel free to continue to do so under the
terms of your existing license. Also note that v0.3.8.6.2 (the most
recent open source version) is still
available on PyPi and can be installed using:
```
pip install antpack==0.3.8.6.2
```
Please note that v0.3.8.6.2 contains a bug fix from v0.3.8.5 and
should be preferred to v0.3.8.5.


## Installation

Starting with v0.3.9, AntPack is only available for noncommercial
academic use. To use AntPack versions
v0.3.9 or later, you must first obtain a [free license key](https://antpackdocumentationlatest.pages.dev/installation).
To install AntPack, run:
```
pip install antpack
AntPack-setup
```

then paste in the license key and your email address when prompted.

The only required dependency is numpy. If you want to run the GUI,
there are two additional dependencies you'll need to install:
```
pip install pyside6 qt_material
```

If you don't plan to use the GUI, you don't need those dependencies.

AntPack is distributed as a wheel precompiled for Linux and Windows and CPython >= 3.9. Support for v0.3.8 is deprecated.

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


### Clustering and fast database search

![Static Badge](https://img.shields.io/badge/New!-grey) Starting in v0.4,
AntPack contains tools for quickly clustering small sequence datasets using
all or part of the antibody sequence. Starting in v0.5, AntPack further
includes tools to quickly build large antibody / TCR sequence databases
and to run lightning-fast searches for sequences similar to queries on
those databases (thousands of searches against a 100 million sequence
database in seconds). AntPack now provides some limited capabilities
for clustering / clonotyping of large datasets. The search tool is
fast enough to run single linkage clustering on 10 - 50 million
sequences. For still larger datasets, AntPack currently provides
a highly scalable mixture model. We hope to further improve the
scalability of single linkage clustering to reach the billion
sequence range in the near future.


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

The authors have also made use of the excellent [parallel-hashmap library by Gregory Popovitch](https://github.com/greg7mdp/parallel-hashmap) in the development of this toolkit.
