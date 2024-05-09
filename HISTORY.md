### Version 0.0.1
Initial release.

### Version 0.0.2
Added Kabat, Martin schemes.

### Version 0.1.0
Added sequence scoring functionality for humanness.
Updated all numbering tools to improve behavior in
the presence of N-terminal deletions.

### Version 0.1.1
Fixed bug that caused gapped scoring (scoring where all
gaps are ignored) to yield masked scores only.

### Version 0.1.2
Added scoring with a user-supplied mask, so that the
user can filter out regions they would like to not include
in the score. Also default to -O3 when compiling.

### Version 0.1.5
Minor updates to alignment; fixed a bug relevant to old versions
of GCC that prevents compilation with old versions of GCC; removed
the scipy dependency; fixed the double-normalization issue.

### Version 0.2.0
Added MultiChainAnnotator, LiabilitySearchTool and VDJAssignmentTool. This
considerably expands the capabilities of the library to cover most basic
antibody amino acid sequence analysis tasks. Removed BioPython as a dependency
(it is not necessary to include a fasta file generator in SingleChainAnnotator
since any users using Biopython can easily add this themselves.)
