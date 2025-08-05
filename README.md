# AntPack

+ [![cov](https://<you>.github.io/<repo>/badges/coverage.svg)](https://github.com/jlparki/antpack/actions)

AntPack is a Python package / toolkit for antibody numbering, data processing,
statistical inference and machine learning for antibody sequences and TCRs.
It is currently in active development, so we are adding new features
and making more improvements periodically.

## Usage

For usage for v0.3.9 and later, see the [docs at this link.](https://antpackdocumentationlatest.pages.dev/)

[You can find the docs here for older versions.](https://antpackdocumentation.pages.dev/)

## What's new in v0.4

v0.4 contains a bug fix for VJ gene assignment for numbering schemes
*other* than IMGT. If perforiming VJ gene assignment with numbering
schemes *other* than IMGT, please prefer v0.4 to prior versions.

v0.4 also contains capabilities for clustering antibody datasets
ranging in size from a few dozen sequences to tens of millions
using a variety of algorithms (see below and see docs). These
will be significantly expanded in the next version with further
clustering and search capabilities and structure analysis-related
tools.

## Licensing

Versions v0.3.9 and afterwards are licensed for academic and noncommercial
use only; you must first obtain a [free license key](https://pwslicensekey.pythonanywhere.com/)
to use. You do not have to re-setup the key when upgrading AntPack in
an existing venv or conda environment.

Versions prior to 0.3.9 were made available under the GPL, which means they
can be used it for your own data analysis in any manner you wish
whether you work in academia or industry,
but any software built using AntPack and intended for sale or distribution 
must also be open-source under the GPL license. Therefore if
you are currently using e.g. v0.3.8.5 for data analysis,
you should feel free to continue to do so under the
terms of your existing license. Also note that v0.3.8.5 is still
available on PyPi and can be installed using:
```
pip install antpack==0.3.8.5
```

## Installation

Starting with v0.3.9, AntPack is only available for noncommercial
academic use. To use AntPack versions
v0.3.9 or later, you must first obtain a [free license key](https://pwslicensekey.pythonanywhere.com/).
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

AntPack is distributed as a wheel precompiled for most platforms and CPython >= 3.8.

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

### Clustering

Starting in v0.4, AntPack contains tools for quickly clustering large sequence
datasets using all or part of the antibody sequence.


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
