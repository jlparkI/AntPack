#ifndef IG_NUMBERING_CONSTANTS_HEADER_H
#define IG_NUMBERING_CONSTANTS_HEADER_H


/*************************************************
 * There are numerous "magic numbers" involved in antibody
 * numbering -- the start and end of the CDRs in each numbering
 * scheme, for example -- which have to be defined somewhere.
 * This header ensures they are all in one convenient place.
 *************************************************/






// The IMGT numbering system will always have (at least) 128 positions.
// Technically light chains also have 128, but due to another weird quirk
// position 128 is never used for light chains.
#define NUM_HEAVY_IMGT_POSITIONS 128
#define NUM_LIGHT_IMGT_POSITIONS 127

// Expected number of positions for Martin, Kabat.
#define NUM_HEAVY_MARTIN_KABAT_POSITIONS 113
#define NUM_LIGHT_MARTIN_KABAT_POSITIONS 107


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
// These are the IMGTs - 1.
#define IMGT_CDR_BREAKPOINT_1 26
#define IMGT_CDR_BREAKPOINT_2 38
#define IMGT_CDR_BREAKPOINT_3 55
#define IMGT_CDR_BREAKPOINT_4 65
#define IMGT_CDR_BREAKPOINT_5 104
#define IMGT_CDR_BREAKPOINT_6 117


// The start and end of CDRs in Kabat (ignoring insertions).
// These are the #s - 1.
#define KABAT_HEAVY_CDR_BREAKPOINT_1 30
#define KABAT_HEAVY_CDR_BREAKPOINT_2 35
#define KABAT_HEAVY_CDR_BREAKPOINT_3 49
#define KABAT_HEAVY_CDR_BREAKPOINT_4 65
#define KABAT_HEAVY_CDR_BREAKPOINT_5 94
#define KABAT_HEAVY_CDR_BREAKPOINT_6 102

#define KABAT_LIGHT_CDR_BREAKPOINT_1 23
#define KABAT_LIGHT_CDR_BREAKPOINT_2 34
#define KABAT_LIGHT_CDR_BREAKPOINT_3 49
#define KABAT_LIGHT_CDR_BREAKPOINT_4 56
#define KABAT_LIGHT_CDR_BREAKPOINT_5 88
#define KABAT_LIGHT_CDR_BREAKPOINT_6 97

// The start and end of CDRs in Martin (ignoring insertions).
// These are the #s - 1.
#define MARTIN_HEAVY_CDR_BREAKPOINT_1 25
#define MARTIN_HEAVY_CDR_BREAKPOINT_2 32
#define MARTIN_HEAVY_CDR_BREAKPOINT_3 51
#define MARTIN_HEAVY_CDR_BREAKPOINT_4 56
#define MARTIN_HEAVY_CDR_BREAKPOINT_5 94
#define MARTIN_HEAVY_CDR_BREAKPOINT_6 102

#define MARTIN_LIGHT_CDR_BREAKPOINT_1 25
#define MARTIN_LIGHT_CDR_BREAKPOINT_2 32
#define MARTIN_LIGHT_CDR_BREAKPOINT_3 49
#define MARTIN_LIGHT_CDR_BREAKPOINT_4 52
#define MARTIN_LIGHT_CDR_BREAKPOINT_5 90
#define MARTIN_LIGHT_CDR_BREAKPOINT_6 96



#endif
