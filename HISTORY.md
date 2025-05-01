### Version 0.3.8.5
Minor updates to the LiabilitySearchTool. Updated error messages
from numbering to distinguish between situations where one of
the three most highly conserved residues is altered (the two
cysteines and trytophan) and situations where the FWGxG motif
nearly always (though not ALWAYS) present near the end of the
sequence has been altered. Fixed a bug in the VJGeneTool that
caused it to search human and mouse only when all species were
specified (and should have been searched).

### Version 0.3.8
Updated PairedChainAnnotator, SingleChainAnnotator and VJGeneTool
to accept and work with TCRs. Updated the command line interface
to require species as an option, to allow for writing output to
a fasta file, to allow TCRs as inputs, and to allow processing
either in memory (suitable for smaller datasets, e.g. < 100K seqs)
or on disk (a little slower but preferable for large datasets).
Updated the VJGeneTool to allow the user to supply "unknown" for
the species, in which case all species are checked. Removed support
for the OGRDB database in VJGeneTool since we have found IMGT to
be preferable. Updated the assign_cdr_labels function to accept
a scheme as input, defaulting to "" which indicates use same scheme
as used to create class; if this scheme is DIFFERENT from the one
used to create class, the set of CDR definitions used to assign
cdr labels is different from the original numbering scheme (e.g.
you can now number with IMGT but use Kabat CDR definitions to
assign CDR labels). Updated the LiabilitySearchTool to allow user
to specify a cdr assignment scheme (e.g. Kabat) that is different
from the numbering scheme (e.g. IMGT).

### Version 0.3.7
Add sensible defaults to build_msa and assign_vj_genes.
Added support for DNA sequences by adding tools to determine
the correct reading frame and forward / reverse complement
for a DNA sequence and translate it to AA. Rewrote the API
for the liability search tool to simplify it. Added an
experimental basic GUI for comparing an arbitrary number
of sequences and for comparing a given sequence with the
closest assigned VJ genes. Added support for germline / VJ
gene assignment for llamas (heavy chain only).

### Version 0.3.6.1
Bug fixes (the rare substr error). Added option to build msa to
add unobserved positions. Rewrote algorithm for paired chain annotator
to improve support for single chains supplied by accident (or
intentionally). Added prefiltering for alignment which improves speed
by roughly 2x. Updated CLI to include unobserved positions by default
and to include percent identity. Updated `assign_cdr_labels` to
always use the same numbering scheme as the Annotator object was
assigned upon creation (avoids possible confusion).

### Version 0.3.5.1
Minor changes to docstrings in 0.3.5.

### Version 0.3.5
Revised and updated API for PairedChainAnnotator, SingleChainAnnotator.
Added bug fix for unusual alignment errors. Updated numbering tools to
accept letter X in input. Converted PairedChainAnnotator,
SingleChainAnnotator, VJGeneTool to C++. Updated vj gene assignment
to use either e-value or percent identity and either IMGT or OGRDB,
and to report multiple matches if multiple equivalent matches are
found. Switched C++ wrapping to Nanobind. Rewrote core alignment code
to produce a 30% improvement in numbering speed / efficiency.
Added Aho numbering scheme. Added code to convert between numbering
schemes when assigning VJ so that VJ can be assigned using any
scheme not just IMGT. Added command line interface for easy
use on fasta files from the command line.

### Version 0.2.7
Added cibuildwheel support to build cross-platform automatically.

### Version 0.2.6
Updated the VJ gene tool to return BOTH the assigned V and J genes and the
percent identity.

### Version 0.2.5
Added the sort position codes function to SingleChainAnnotator, which enables
the user to sort a list of position codes (this is useful and convenient when
merging many lists of numbers). Set SingleChainAnnotator up to return a list
indicating whether each position is cdr or framework if the user so requests.

### Version 0.2.0
Added MultiChainAnnotator, LiabilitySearchTool and VJGeneTool. This
considerably expands the capabilities of the library to cover most basic
antibody amino acid sequence analysis tasks. Removed BioPython as a dependency
(it is not necessary to include a fasta file generator in SingleChainAnnotator
since any users using Biopython can easily add this themselves.)

### Version 0.1.5
Minor updates to alignment; fixed a bug relevant to old versions
of GCC that prevents compilation with old versions of GCC; removed
the scipy dependency; fixed the double-normalization issue.

### Version 0.1.2
Added scoring with a user-supplied mask, so that the
user can filter out regions they would like to not include
in the score.

### Version 0.1.1
Fixed bug that caused gapped scoring (scoring where all
gaps are ignored) to yield masked scores only.

### Version 0.1.0
Added sequence scoring functionality for humanness.
Updated all numbering tools to improve behavior in
the presence of N-terminal deletions.

### Version 0.0.2
Added Kabat, Martin schemes.

### Version 0.0.1
Initial release.
