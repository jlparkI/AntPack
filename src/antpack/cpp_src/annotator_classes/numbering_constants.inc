#ifndef IG_NUMBERING_CONSTANTS_INCLUDE_H
#define IG_NUMBERING_CONSTANTS_INCLUDE_H


#include <array>

/*************************************************
 * There are numerous "magic numbers" involved in antibody
 * numbering -- the start and end of the CDRs in each numbering
 * scheme, for example -- which have to be defined somewhere.
 * This header ensures they are all in one convenient place.
 *************************************************/

namespace NumberingTools {

// Codes for the pathways that can link a score
// to the best-scoring parent.
static constexpr int LEFT_TRANSFER = 1;
static constexpr int UP_TRANSFER = 2;
static constexpr int DIAGONAL_TRANSFER = 0;

// The columns of the score matrix that are accessed for gap penalties.
static constexpr int QUERY_GAP_COLUMN = 21;
static constexpr int TEMPLATE_GAP_COLUMN = 22;

// We universally use 1 for valid result, 0 for error.
static constexpr int VALID_SEQUENCE = 1;
static constexpr int INVALID_SEQUENCE = 0;


// The expected number of AAs.
static constexpr int EXPECTED_NUM_AAS_FOR_ALIGNERS = 23;


// The minimum number of amino acids in a sequence to try to align it.
// Less than this and it will be immediately rejected. This is fairly
// arbitrary, we didn't put much thought into the selection of 25 --
// a typical chain is > 100 AAs, so anything MUCH less than that is
// clearly a fragment that probably can't be reliably numbered.
static constexpr int MINIMUM_SEQUENCE_LENGTH = 25;


// The IMGT numbering system will always have (at least) 128 positions.
// Technically light chains also have 128, but due to another weird quirk
// position 128 is never used for light chains.
static constexpr int NUM_HEAVY_IMGT_POSITIONS = 128;
static constexpr int NUM_LIGHT_IMGT_POSITIONS = 127;

// Expected number of positions for Martin, Kabat.
static constexpr int NUM_HEAVY_MARTIN_KABAT_POSITIONS = 113;
static constexpr int NUM_LIGHT_MARTIN_KABAT_POSITIONS = 107;

// Expected number of positions for Aho.
static constexpr int NUM_HEAVY_AHO_POSITIONS = 149;
static constexpr int NUM_LIGHT_AHO_POSITIONS = 149;


// These are "magic number" positions in the IMGT framework at
// which "forwards-backwards" insertion numbering must be applied.
// This is a nuisance, but is out of our control -- the IMGT #ing
// system has this quirk built-in... Note that because IMGT numbers
// from 1, these positions are the actual IMGT position - 1.
static constexpr int CDR1_INSERTION_PT = 32;
static constexpr int CDR2_INSERTION_PT = 60;
static constexpr int CDR3_INSERTION_PT = 110;

// Highly conserved positions in the IMGT scheme. These are the IMGT #s - 1.
static constexpr int HIGHLY_CONSERVED_IMGT_1 = 22;
static constexpr int HIGHLY_CONSERVED_IMGT_2 = 40;
static constexpr int HIGHLY_CONSERVED_IMGT_3 = 103;
static constexpr int HIGHLY_CONSERVED_IMGT_4 = 117;
static constexpr int HIGHLY_CONSERVED_IMGT_5 = 118;
static constexpr int HIGHLY_CONSERVED_IMGT_6 = 120;

// Highly conserved positions in the AHO scheme. These are the AHO #s - 1.
static constexpr int HIGHLY_CONSERVED_AHO_1 = 22;
static constexpr int HIGHLY_CONSERVED_AHO_2 = 42;
static constexpr int HIGHLY_CONSERVED_AHO_3 = 105;
static constexpr int HIGHLY_CONSERVED_AHO_4 = 138;
static constexpr int HIGHLY_CONSERVED_AHO_5 = 139;
static constexpr int HIGHLY_CONSERVED_AHO_6 = 141;


// Highly conserved positions in the Martin / Kabat schemes for heavy chains.
// These are the #s - 1.
static constexpr int HIGHLY_CONSERVED_KABAT_HEAVY_1 = 21;
static constexpr int HIGHLY_CONSERVED_KABAT_HEAVY_2 = 35;
static constexpr int HIGHLY_CONSERVED_KABAT_HEAVY_3 = 91;
static constexpr int HIGHLY_CONSERVED_KABAT_HEAVY_4 = 102;
static constexpr int HIGHLY_CONSERVED_KABAT_HEAVY_5 = 103;
static constexpr int HIGHLY_CONSERVED_KABAT_HEAVY_6 = 105;

// Highly conserved positions in the Martin / Kabat schemes for light chains.
// These are the #s - 1.
static constexpr int HIGHLY_CONSERVED_KABAT_LIGHT_1 = 22;
static constexpr int HIGHLY_CONSERVED_KABAT_LIGHT_2 = 34;
static constexpr int HIGHLY_CONSERVED_KABAT_LIGHT_3 = 87;
static constexpr int HIGHLY_CONSERVED_KABAT_LIGHT_4 = 97;
static constexpr int HIGHLY_CONSERVED_KABAT_LIGHT_5 = 98;
static constexpr int HIGHLY_CONSERVED_KABAT_LIGHT_6 = 100;


// Gap open vs extend penalties for c and n terminals.
static constexpr int NTERMINAL_QUERY_GAP_PENALTY = -1;
static constexpr int CTERMINAL_QUERY_GAP_PENALTY = -1;
static constexpr int CTERMINAL_TEMPLATE_GAP_PENALTY = -1;
static constexpr int MAX_NTERMINAL_GAP_PENALTY = -5;

}  // namespace NumberingTools

#endif
