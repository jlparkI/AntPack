# AntPack

AntPack is a toolkit for antibody numbering, data processing, statistical inference and
machine learning for antibody sequences. For usage,
see [the docs](https://antpack.readthedocs.io/en/latest/index.html).
It is currently in active development, so we are adding new features
and making more improvements periodically. We will try to avoid breaking
changes but nonetheless recommend checking the docs after you install
a new version to be sure the component you are using is unaffected.

### Are there specific features / changes you'd like to see in AntPack?

If so, feel free to take our short [anonymous user survey!](https://www.surveymonkey.com/r/FWQJKZS)
There are quite a few features we plan to add to AntPack --
this survey will help us decide which we should prioritize /
which should come first.


## What's new in v0.3.7

v0.3.7 contains some minor improvements to the API,
specifically for `build_msa`, `assign_vj_genes` and
the liability search tool. It also adds support for
DNA sequences in the form of a tool that determines
the correct reading frame and forward / reverse complement
for an input DNA sequence and translates it to AAs
(this is faster than aligning to / numbering the
DNA sequence). We've added support for llama germline
gene assignment (VHH only). Finally, we've added an
experimental GUI which we plan to expand considerably
in upcoming versions. You can quickly launch the GUI from
the terminal by typing:
```
AntPack-GUI
```

whereas to run the command line interface use:
```
AntPack-CLI
```

For now, you can use the GUI to add and align / compare some arbitrary
number of paired and/or single chain sequences. You
can also retrieve the VJ gene assignments for an input sequence
and see an alignment of your input sequence with those genes.

The command line interface is useful for some fairly standard
analyses of moderate-large datasets, while the more powerful
Python API is useful for building your own workflows and pipelines.
The GUI (at least right now) is mostly useful for a fast analysis of a few
sequences that you want to examine quickly without needing to
start a Jupyter notebook and write code.

## Coming soon

We're planning to continue improving AntPack's numbering capabilities
and fixing "edge cases" so that we keep improving speed and quality.
We're also looking at adding support for TCR sequences, at least for
numbering. We plan to extensively expand the GUI. We're also planning to link
AntPack to a searchable online database of human antibody sequences
and add tools for fast sequence clustering and stability prediction --
coming soon.

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
of numbered sequences into an MSA and for easy extraction of specific CDRs and
framework regions.


#### V / J genes

Identifying the most similar human/mouse/llama V / J gene sequences is useful
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
