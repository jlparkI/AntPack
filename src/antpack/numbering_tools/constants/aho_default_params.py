"""Defines important defaults for the AHO numbering scheme."""
import copy

#The expected length of the heavy and light chains
NUM_HEAVY = 149
NUM_LIGHT = 149

############################The following are alignment parameters that
#were selected to try to ensure a good alignment. Making large changes
#to these WILL change the results, possibly substantially.

DEFAULT_C_TERMINAL_QUERY_GAP_PENALTY = -1
DEFAULT_TERMINAL_TEMPLATE_GAP_PENALTY = -1

#These are positions that are essentially mandatory. A sequence that deviates
#from one of these is an alignment issue or has a very large deletion. Note that
#AHO defines position 90 as highly conserved in the sense that a hydrophobic
#amino acid is always observed here, but "hydrophobic amino acid" is not defined
#and a variety of (hydrophobic) aas are observed at that position.
heavy_conserved_positions = {23:"C", 43:"W", 106:"C", 139:"W", 140:"G", 142:"G"}
light_conserved_positions = {23:"C", 43:"W", 106:"C", 139:"F", 140:"G", 142:"G"}


#These are positions where either A) a blank in the query sequence is very common or
#B) insertions are very common. We want special template and query gap penalties for
#these positions. One or two of these are chain dependent. The first value is
#the query gap column, the second is the template gap column.
heavy_special_positions = {8:[-1,-12], 27:[-1,-12], 28:[-1,-12],
        36:[-1.0,-1.0], 63:[-1.0,-1.0],
        74:[-1.,-12.], 75:[-1.,-12.],
        85:[-1.,-12.], 86:[-1.,-12.],
        123:[-1.0,-1.0]}

light_special_positions = {8:[-1,-12], 27:[-1,-12], 28:[-1,-12],
        36:[-1.0,-1.0], 63:[-1.0,-1.0],
        74:[-1.,-12.], 75:[-1.,-12.],
        85:[-1.,-12.], 86:[-1.,-12.],
        123:[-1.0,-1.0]}

#Aho-defined CDRs. We have a specially defined gap penalty for each to ensure they are filled
#in the correct order. For Aho, light and heavy cdrs have the same definition. Notice that
#positions here do not exactly correspond to the Aho CDR definition, but rather to the region
#where insertions / deletions are likely.
heavy_cdrs = {29:-5.1, 30:-4.6, 31:-4.1, 32:-3.6, 33:-3.1, 34:-2.6, 35:-2.5, 36:-3, 37:-3.5,
        38:-4, 39:-4.5, 40:-5,
        # CDR 2
        56:-5.1, 57:-4.6, 58:-4.1, 59:-3.6, 60:-3.1, 61:-3, 62:-3.5, 63:-4, 64:-4.5, 65:-5,
        #CDR 3
        105:-5.1, 106:-4.6, 107:-4.1, 108:-3.6, 109:-3.1, 110:-2.6, 111:-1,
        112:-2.5, 113:-3, 114:-3.5, 115:-4, 116:-4.5, 117:-5}
light_cdrs = copy.copy(heavy_cdrs)
