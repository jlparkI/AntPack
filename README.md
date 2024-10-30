# AntPack

AntPack is a toolkit for antibody numbering, data processing, statistical inference and
machine learning for antibody sequences. For usage,
see [the docs](https://antpack.readthedocs.io/en/latest/index.html).
It is currently in active development, so we are adding new features
and making more improvements periodically. We will try to avoid breaking
changes but nonetheless recommend checking the docs after you install
a new version to be sure the component you are using is unaffected.


## What's new in v0.3.6

v0.3.6 incorporates a couple of bug fixes and some improvements to
the command line tool. It's also about 2x faster than the previous
version. On an Intel-i7-13700K, the time taken to number 3500
sequences for different tools was as follows:

| Tool     | Time (s)        |
| -------- | --------------- |
| ANARCI   | 45    +/- 1     |
| AntPack  | 0.18  +/- 0.01  |
| AbRSA    | 13.02 +-/ 0.01  |

Previous versions of AntPack achieved about 0.35 seconds. Finally,
in previous versions you had to know whether the sequences you
were using were single or paired. AntPack still has the
SingleChainAnnotator and PairedChainAnnotator *but* PairedChainAnnotator
has been made more flexible and better at dealing with
single chains, so if your input data contains both paired and single
chains, it's ok to default to using PairedChainAnnotator. (SingleChain
will be much faster, so if you know you're dealing with single chains
only, always use SingleChainAnnotator instead.) You can also use
SingleChainAnnotator on chains that may be single or paired; if you
set it to look for a specific chain type (e.g. heavy), it will
extract that chain from paired sequences.

## Coming soon

We're adding a GUI that can be launched with a single command
from the terminal / Windows command line so that AntPack can be
used either as a command line tool, a Python library or a GUI,
depending on your needs and what you're comfortable with.

We're also planning to link AntPack to a searchable online database of
human antibody sequences and add tools for fast sequence clustering
and stability prediction -- coming soon.

## Installation

```
pip install antpack
```

AntPack is distributed as a wheel precompiled for most platforms and CPython >= 3.8,
so installation should be very straightforward. A source distribution is also available
(C++17) in case there is any need to compile from source.

The only required dependency is numpy. In future, to use the
GUI there will be a couple of additional optional dependencies
(although this is not added yet).

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
