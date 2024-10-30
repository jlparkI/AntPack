"""Defines important defaults for the AHO numbering scheme."""
import copy

#The expected length of the heavy and light chains
NUM_HEAVY = 149
NUM_LIGHT = 149

############################The following are alignment parameters that
#were selected to try to ensure a good alignment. Making large changes
#to these WILL change the results, possibly substantially.

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
#the query gap column, the second is the template gap column. Aho unlike other
#schemes has a separate special positions dictionary for kappa.
heavy_special_positions = {8:[-1,-12], 27:[-3.7,-12], 28:[-1.0,-12],
        36:[-1.1,-1.1], 63:[-1.1,-1.1],
        74:[-10.,-12.], 75:[-10.,-12.],
        85:[-1.5,-12.], 86:[-1.,-12.],
        87:[-8., -12.], 123:[-1.1,-1.1]}

light_special_positions = {8:[-1,-12], 27:[-3.8,-12.],
        28:[-0.2,-12.], 29:[-3.9,-12.],
        36:[-1.1,-1.1], 63:[-1.1,-1.1],
        74:[-10.,-12.], 75:[-10.,-12.],
        85:[-1.5,-12.], 86:[-1.1,-12.],
        87:[-8., -12.], 123:[-1.1,-1.1]}

kappa_special_positions = {8:[-1,-12], 27:[-0.9,-12.], 28:[-0.2,-12.],
        36:[-1.1,-1.1], 63:[-1.1,-1.1],
        74:[-10.,-12.], 75:[-10.,-12.],
        85:[-1.5,-12.], 86:[-1.,-12.],
        87:[-8., -12.], 123:[-1.05,-1.05]}

#Aho-defined CDRs. We have a specially defined gap penalty for each to ensure they are filled
#in the correct order. For Aho, light and heavy cdrs have the same definition. Notice that
#positions here do not exactly correspond to the Aho CDR definition, but rather to the region
#where insertions / deletions are likely.
heavy_cdrs = {26:-5.5, 30:-6., 31:-5.0, 32:-4.5, 33:-4.0, 34:-3.5, 35:-3.0, 36:-2.5, 37:-3.1,
        38:-3.6, 39:-4.1, 40:-4.6,
        # CDR 2
        58:-5., 59:-4.5, 60:-4., 61:-3.5, 62:-3., 63:-2.6,
        64: -3.1, 65:-3.6, 66:-4.1, 67:-4.6, 68:-5.1,
        #CDR 3
        109:-5.1, 110:-4.8, 111:-4.5, 112:-4.2, 113:-3.9, 114:-3.6, 115:-3.3, 116:-3.0,
        117:-2.7, 118:-2.4, 119:-2.1, 120:-1.8, 121:-1.5, 122:-1.2, 123:-1, 124:-1.1,
        125:-1.4, 126:-1.7, 127:-2.0, 128:-2.3, 129:-2.6, 130:-2.9, 131:-3.2, 132:-3.5,
        133:-3.8, 134:-4.1, 135:-4.4, 136:-4.7, 137:-5 }
light_cdrs = copy.copy(heavy_cdrs)
