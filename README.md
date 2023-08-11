# AntPack

AntPack is a toolkit for processing, annotating and making inferences
about antibody sequences. Currently it's in VERY early development,
so it may look like a mess for a while. The following are some forward-
looking statements about what we hope to incorporate into AntPack.


### Antibody numbering

The situation surrounding antibody numbering is moderately confusing since
over 7 different schemes have been described in the literature. To make
matters worse, of the tools developed for antibody sequence numbering
(AbRSA, AbNum, DIGIT, ANARCI, IMGT DomainGapAlign) most are closed-source
and are made available as webservers -- hence, they are hard to incorporate
or build into a pipeline. The open-source options have some room for improvement.
For example, the ANARCI tool is a pure Python script that if supplied
with a fasta file of input sequences will

1) load all the sequences into memory at once
2) write them all back to disk in a new location
3) Run this through HMMER, bouncing it off of a db containing all available species
and chain types, even if only one species or chain type is of interest;
4) read the HMMER output from disk back into memory
5) Make a series of corrections to the HMMER alignment
6) For single-domain sequences, rewrite all of them to disk again;
7) Run the single-domain sequences past all HMMER dbs again and read the results
from disk back into memory;
8) Use this second round of results to fix issues with the first round

and then do the numbering in pure Python. To be fair, the ANARCI tool was (like the
other tools) most probably a set of quick and dirty scripts intended for light
use; the efficiency doesn't matter much if you're just processing a few 100s to 1000s of
sequences. However, for larger volumes and more complicated projects, a more
sophisticated and efficient tool would be helpful, and that's one thing we'll
introduce here.


### Mapping antibody sequence space

We've built a set of human-interpretable statistical models of antibody sequence space --
more on this project soon!


### Human-ness and developability

Immunogenicity is a frequent cause of failure for antibody drug candidates in clinical trials.
Avoiding such costly failures is highly desirable. As part of this package, we'll introduce
tools to help predict the likelihood of immunogenicity or other serious developability
problems (coming soon!)
