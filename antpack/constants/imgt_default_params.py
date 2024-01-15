"""Defines important defaults for the IMGT numbering scheme."""
import copy

#The expected length of the heavy and light chains
NUM_HEAVY = 128
NUM_LIGHT = 127

############################The following are alignment parameters that
#were selected to try to ensure a good alignment. Making large changes
#to these WILL change the results, possibly substantially.

DEFAULT_C_TERMINAL_QUERY_GAP_PENALTY = -1
DEFAULT_TERMINAL_TEMPLATE_GAP_PENALTY = -1

#These are positions that are essentially mandatory. A sequence that deviates
#from one of these is an alignment issue or has a very large deletion. Note that
#IMGT defines position 89 as highly conserved in the sense that a hydrophobic
#amino acid is always observed here, but "hydrophobic amino acid" is not defined
#and a variety of (hydrophobic) aas are observed at that position. Hence we
#handle position 89 separately below.
heavy_conserved_positions = {23:"C", 41:"W", 104:"C", 118:"W", 119:"G", 121:"G"}
light_conserved_positions = {23:"C", 41:"W", 104:"C", 118:"F", 119:"G", 121:"G"}


#These are positions where either A) a blank in the query sequence is very common or
#B) insertions are very common. We want special template and query gap penalties for
#these positions. One or two of these are chain dependent. The first value is
#the query gap column, the second is the template gap column.
heavy_special_positions = {10:[-1,-12], 33:[-1.0,-1.0], 61:[-1.0,-1.0],
        73:[-1,-11], 111:[-1.0,-1.0]}

light_special_positions = {10:[-1,-12],
        33:[-1.0,-1.0], 61:[-1.0,-1.0], 73:[-1,-11], 81:[-1.0,-25],
        82:[-1.0,-25], 111:[-1.0,-1.0]}

#IMGT-defined CDRs. We have a specially defined gap penalty for each to ensure they are filled
#in the correct order. For IMGT, light and heavy cdrs have the same definition.
heavy_cdrs = {27:-5.1, 28:-4.6, 29:-4.1, 30:-3.6, 31:-3.1, 32:-2.6, 33:-2.5, 34:-3, 35:-3.5,
        36:-4, 37:-4.5, 38:-5,
        # CDR 2
        56:-5.1, 57:-4.6, 58:-4.1, 59:-3.6, 60:-3.1, 61:-3, 62:-3.5, 63:-4, 64:-4.5, 65:-5,
        #CDR 3
        105:-5.1, 106:-4.6, 107:-4.1, 108:-3.6, 109:-3.1, 110:-2.6, 111:-1,
        112:-2.5, 113:-3, 114:-3.5, 115:-4, 116:-4.5, 117:-5}
light_cdrs = copy.copy(heavy_cdrs)
