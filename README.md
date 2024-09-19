# AntPack

AntPack is a toolkit for antibody numbering, data processing, statistical inference and
machine learning for antibody sequences. For usage,
see [the docs](https://antpack.readthedocs.io/en/latest/index.html).


## What's new in v0.3.5

In v0.3.5, we've added the Aho numbering scheme. We've made AntPack numbering
about 30% faster (it was already the fastest numbering tool, but there's no
such thing as too fast!) We've also improved VJ 
gene assignment. AntPack can now assign v- and j-genes using either percent
identity or e-value, and if multiple V- or J-genes have essentially the
same similarity to the query sequence, it will return a list of them
rather than trying to assign a single V or J gene (which doesn't make much
sense in that case).

The API has been slightly redesigned. This may cause some breaking changes,
mainly for VJ gene assignment. Our apologies for any inconvenience -- the
API should be stable from here on out. Finally, for numbering and VJ gene
assignment, AntPack now accepts sequences which contain the letter 'X'
(although we suggest using this only if necessary -- an alignment with
all AAs will generally be better quality).

Finally, AntPack now has a command line tool for running a quick standard
analysis on small-medium size datasets. For more customized analyses and/or
larger datasets the Python API can be used as before.

Python3.8 or later is now required.

Starting in v0.4 we're planning to add a GUI -- currently in active
development -- which makes it easy to load and view sequences if you only
have a few. We're also planning to add additional functionality for antibody
analysis -- coming soon!


## Installation

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
