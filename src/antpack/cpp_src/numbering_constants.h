#ifndef IG_NUMBERING_CONSTANTS_HEADER_H
#define IG_NUMBERING_CONSTANTS_HEADER_H


/*************************************************
 * There are numerous "magic numbers" involved in antibody
 * numbering -- the start and end of the CDRs in each numbering
 * scheme, for example -- which have to be defined somewhere.
 * This header ensures they are all in one convenient place.
 *************************************************/



// The minimum number of amino acids in a sequence to try to align it.
// Less than this and it will be immediately rejected. This is fairly
// arbitrary, we didn't put much thought into the selection of 25 --
// a typical chain is > 100 AAs, so anything MUCH less than that is
// clearly a fragment that probably can't be reliably numbered.
#define MINIMUM_SEQUENCE_LENGTH 25


// The amount by which to offset from the start of a c-terminal region
// when using tools that look for typical J-gene sequences.
#define CTERMINAL_STANDARD_OFFSET 9


// The IMGT numbering system will always have (at least) 128 positions.
// Technically light chains also have 128, but due to another weird quirk
// position 128 is never used for light chains.
#define NUM_HEAVY_IMGT_POSITIONS 128
#define NUM_LIGHT_IMGT_POSITIONS 127

// Expected number of positions for Martin, Kabat.
#define NUM_HEAVY_MARTIN_KABAT_POSITIONS 113
#define NUM_LIGHT_MARTIN_KABAT_POSITIONS 107

// Expected number of positions for Aho.
#define NUM_HEAVY_AHO_POSITIONS 149
#define NUM_LIGHT_AHO_POSITIONS 149


// These are "magic number" positions in the IMGT framework at
// which "forwards-backwards" insertion numbering must be applied.
// This is a nuisance, but is out of our control -- the IMGT #ing
// system has this quirk built-in... Note that because IMGT numbers
// from 1, these positions are the actual IMGT position - 1.
#define CDR1_INSERTION_PT 32
#define CDR2_INSERTION_PT 60
#define CDR3_INSERTION_PT 110

// Highly conserved positions in the IMGT scheme. These are the IMGT #s - 1.
#define HIGHLY_CONSERVED_IMGT_1 22
#define HIGHLY_CONSERVED_IMGT_2 40
#define HIGHLY_CONSERVED_IMGT_3 103
#define HIGHLY_CONSERVED_IMGT_4 117
#define HIGHLY_CONSERVED_IMGT_5 118
#define HIGHLY_CONSERVED_IMGT_6 120

// Highly conserved positions in the AHO scheme. These are the AHO #s - 1.
#define HIGHLY_CONSERVED_AHO_1 22
#define HIGHLY_CONSERVED_AHO_2 42
#define HIGHLY_CONSERVED_AHO_3 105
#define HIGHLY_CONSERVED_AHO_4 138
#define HIGHLY_CONSERVED_AHO_5 139
#define HIGHLY_CONSERVED_AHO_6 141


// Highly conserved positions in the Martin / Kabat schemes for heavy chains.
// These are the #s - 1.
#define HIGHLY_CONSERVED_KABAT_HEAVY_1 21
#define HIGHLY_CONSERVED_KABAT_HEAVY_2 35
#define HIGHLY_CONSERVED_KABAT_HEAVY_3 91
#define HIGHLY_CONSERVED_KABAT_HEAVY_4 102
#define HIGHLY_CONSERVED_KABAT_HEAVY_5 103
#define HIGHLY_CONSERVED_KABAT_HEAVY_6 105

// Highly conserved positions in the Martin / Kabat schemes for light chains.
// These are the #s - 1.
#define HIGHLY_CONSERVED_KABAT_LIGHT_1 22
#define HIGHLY_CONSERVED_KABAT_LIGHT_2 34
#define HIGHLY_CONSERVED_KABAT_LIGHT_3 87
#define HIGHLY_CONSERVED_KABAT_LIGHT_4 97
#define HIGHLY_CONSERVED_KABAT_LIGHT_5 98
#define HIGHLY_CONSERVED_KABAT_LIGHT_6 100


// The start and end of CDRs in IMGT (ignoring insertions).
#define IMGT_CDR_BREAKPOINT_1 27
#define IMGT_CDR_BREAKPOINT_2 39
#define IMGT_CDR_BREAKPOINT_3 56
#define IMGT_CDR_BREAKPOINT_4 66
#define IMGT_CDR_BREAKPOINT_5 105
#define IMGT_CDR_BREAKPOINT_6 118

// The start and end of CDRs in AHO (ignoring insertions).
#define AHO_CDR_BREAKPOINT_1 25
#define AHO_CDR_BREAKPOINT_2 41
#define AHO_CDR_BREAKPOINT_3 58
#define AHO_CDR_BREAKPOINT_4 78
#define AHO_CDR_BREAKPOINT_5 109
#define AHO_CDR_BREAKPOINT_6 138


// The start and end of CDRs in Kabat (ignoring insertions).
#define KABAT_HEAVY_CDR_BREAKPOINT_1 31
#define KABAT_HEAVY_CDR_BREAKPOINT_2 36
#define KABAT_HEAVY_CDR_BREAKPOINT_3 50
#define KABAT_HEAVY_CDR_BREAKPOINT_4 66
#define KABAT_HEAVY_CDR_BREAKPOINT_5 95
#define KABAT_HEAVY_CDR_BREAKPOINT_6 103

#define KABAT_LIGHT_CDR_BREAKPOINT_1 24
#define KABAT_LIGHT_CDR_BREAKPOINT_2 35
#define KABAT_LIGHT_CDR_BREAKPOINT_3 50
#define KABAT_LIGHT_CDR_BREAKPOINT_4 57
#define KABAT_LIGHT_CDR_BREAKPOINT_5 89
#define KABAT_LIGHT_CDR_BREAKPOINT_6 98

// The start and end of CDRs in Martin (ignoring insertions).
#define MARTIN_HEAVY_CDR_BREAKPOINT_1 26
#define MARTIN_HEAVY_CDR_BREAKPOINT_2 33
#define MARTIN_HEAVY_CDR_BREAKPOINT_3 52
#define MARTIN_HEAVY_CDR_BREAKPOINT_4 57
#define MARTIN_HEAVY_CDR_BREAKPOINT_5 95
#define MARTIN_HEAVY_CDR_BREAKPOINT_6 103

#define MARTIN_LIGHT_CDR_BREAKPOINT_1 26
#define MARTIN_LIGHT_CDR_BREAKPOINT_2 33
#define MARTIN_LIGHT_CDR_BREAKPOINT_3 50
#define MARTIN_LIGHT_CDR_BREAKPOINT_4 53
#define MARTIN_LIGHT_CDR_BREAKPOINT_5 91
#define MARTIN_LIGHT_CDR_BREAKPOINT_6 97



#endif
