# AntPack

AntPack is a toolkit for antibody numbering, data processing, statistical inference and
machine learning for antibody sequences. For usage,
see [the docs](https://antpack.readthedocs.io/en/latest/index.html).


## What's new in v0.3

v0.3 provides some API redesign for greater ease of use and convenience. The
MultiChainAnnotator has been replaced with the PairedChainAnnotator, which
is designed specifically for paired light and heavy chains. Both this tool
and SingleChainAnnotator now offer a `build_msa` command that makes it easy
to combine many numbered sequences into an fixed-length array for further
analysis / writing to file.

Coming soon: we're planning to add a GUI front-end you can use to review
sequences numbered by AntPack, humanize sequences and determine how different
mutations will impact predicted sequence properties. For more on this project
as it matures, [see this repo](https://github.com/jlparkI/RESP_proteins).


# Installation

```
pip install antpack
```

AntPack is distributed as a wheel precompiled for most platforms and CPython >= 3.7,
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


#### Citing this work

If using AntPack in research intended for publication, please cite
either the preprint:

Jonathan Parkinson and Wei Wang. 2024. For antibody sequence generative modeling,
mixture models may be all you need.
[bioRxiv: https://doi.org/10.1101/2024.01.27.577555](https://www.biorxiv.org/content/10.1101/2024.01.27.577555v1)

or the final paper in [Bioinformatics.](https://academic.oup.com/bioinformatics/article/40/5/btae278/7656770)
