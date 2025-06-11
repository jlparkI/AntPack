# AntPack

AntPack is a Python package / toolkit for antibody numbering, data processing,
statistical inference and machine learning for antibody sequences and TCRs.
It is currently in active development, so we are adding new features
and making more improvements periodically.

## Usage

For usage for v0.3.9 and later, see the [docs at this link.](https://antpackdocumentationlatest.pages.dev/)

[You can find the docs here for older versions.](https://antpackdocumentation.pages.dev/)

## Major updates in v0.3.9

As we continue to add features to AntPack, we are finding some
features that we will add in upcoming versions useful for our own
discovery efforts. Consequently we have reluctantly decided starting
with v0.3.9 to license AntPack for noncommercial use only. This was
not an easy decision and not one that we took lightly.

Prior versions were made available under the GPL, which means they
can be used it for your own data analysis in any manner you wish
whether you work in academia or industry,
but any software built using AntPack and intended for sale or distribution 
must also be open-source under the GPL license. Therefore if
you are currently using e.g. v0.3.8 for data analysis,
you should feel free to continue to do so under the
terms of your existing license.

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
however, there are two additional dependencies you'll need to install
(in future versions these dependencies may no longer be needed):
```
pip install pyside6 qt_material
```

If you don't plan to use the GUI, you don't need those dependencies.

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
