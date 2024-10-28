# AntPack

AntPack is a toolkit for antibody numbering, data processing, statistical inference and
machine learning for antibody sequences. For usage,
see [the docs](https://antpack.readthedocs.io/en/latest/index.html).
It is currently in active development, so we are adding new features
and making more improvements periodically. We will try to avoid breaking
changes but nonetheless recommend checking the docs after you install
a new version to be sure your code is unaffected.


## What's new in v0.3.6

v0.3.6 incorporates a couple of bug fixes and a for now very basic
graphical user interface (GUI) (we're planning to add more features
to this soon). The GUI can be launched with a single command
from the terminal / Windows command line. AntPack can now therefore be used
either from the command line, from Python or from the GUI. The GUI
is convenient if you just want to take a quick look at a small
number of sequences. The command line tool is ideal if you want to
run a fairly simple standard analysis. The Python API is best if
you want to build your own pipeline, need to do a more custom analysis
and/or are comfortable writing scripts / working in Jupyter Notebook.

## Upcoming in planned versions

We're planning to make some extensive additions to the GUI, to refactor
the numbering code for improved efficiency, and to add a ChainAnnotator
to the Python API. Currently AntPack has SingleChainAnnotator if you
know your sequence contains one variable region or PairedChainAnnotator
if you know it contains heavy and light. ChainAnnotator will extract
however many variable regions are in your sequence with no advance
knowledge needed about number of regions present.

We're also planning to link AntPack to a searchable online database of
human antibody sequences -- coming soon.


## Installation

```
pip install antpack
```

AntPack is distributed as a wheel precompiled for most platforms and CPython >= 3.8,
so installation should be very straightforward. A source distribution is also available
(C++17) in case there is any need to compile from source.

The only required dependency is numpy. If you want to use the GUI, however, you will
need to install two additional optional dependencies (optional because they are required
only for the GUI and not for other AntPack functionality).

## Capabilities


#### Antibody numbering

Numbering antibody sequences is an important precursor for many statistical inference /
machine learning applications. AntPack is orders of magnitude faster for numbering
antibody sequences than existing tools in the literature (e.g. ANARCI, AbRSA),
while providing >= reliability. AntPack also provides tools for merging a list
of numbered sequences into an MSA and for easy extraction of specific CDRs and
framework regions.


#### V / J genes

Identifying the most similar human V / J gene sequences is useful for a variety of
purposes. AntPack provides tools for determining which human V and J gene sequences
are most similar to the variable region chain provided as input.


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


#### Licensing

AntPack is licensed under a GPL license, which means that you are free to
use it for your own data analysis irrespective of whether you work on
academic research or industrial R&D / QC.

If you are writing software intended for sale / distribution, however,
any software you create that uses AntPack must also be open-source, must
use the GPL license and must acknowledge AntPack appropriately. If you
are interested in using AntPack in a closed-source software product,
please contact us to obtain a version of AntPack licensed for this use.


#### Citing this work

If using AntPack in research intended for publication, please cite
either the preprint:

Jonathan Parkinson and Wei Wang. 2024. For antibody sequence generative modeling,
mixture models may be all you need.
[bioRxiv: https://doi.org/10.1101/2024.01.27.577555](https://www.biorxiv.org/content/10.1101/2024.01.27.577555v1)

or the final paper in [Bioinformatics.](https://academic.oup.com/bioinformatics/article/40/5/btae278/7656770)
