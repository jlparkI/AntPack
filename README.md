# AntPack

AntPack is a toolkit for antibody numbering, data processing, statistical inference and
machine learning for antibody sequences. For usage,
see [the docs](https://antpack.readthedocs.io/en/latest/index.html).


## What's new in v0.3

In v0.3.5, we've added the Aho numbering scheme. We've also improved the
accuracy of VJ gene assignment -- previously AntPack used an algorithm
very similar to ANARCI for VJ gene assignment which is suboptimal. Starting
in v0.3.5 we are using an improved algorithm which offers greater accuracy.
The API has been slightly redesigned, with VJ gene assignment merged with
the numbering tools for greater ease of use. We've also added some multithreading
to further speed up sequence numbering.

Finally, starting in v0.3.5 we're providing a GUI -- currently in active
development -- which makes it easy to load and view sequences if you only
have a few. You can align / number sequences, see how and where they differ,
view predicted properties of each sequence, and see which regions are least
human / should be changed in order to humanize a given sequence. The GUI
is handy for qualitative analysis of a few sequences, while
the Python package is better for building your own tools and workflows.
We're planning to add more functionality to this (e.g. sequence clustering)
and to AntPack generally -- coming soon!


## :hammer_and_wrench: Installation

```
pip install antpack
```

AntPack is distributed as a wheel precompiled for most platforms and CPython >= 3.7,
so installation should be very straightforward. A source distribution is also available
(C++17) in case there is any need to compile from source.

Starting with v0.3.5, AntPack can also be used / accessed as a GUI tool. If you wish
to make use of this tool, you'll have to follow a separate installation procedure;
for more on this, see [the docs](https://antpack.readthedocs.io/en/latest/index.html).

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
