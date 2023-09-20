"""Defines important defaults for the Kabat numbering scheme."""

#The expected length of the heavy and light chains
NUM_HEAVY = 113
NUM_LIGHT = 107

############################The following are alignment parameters that
#were selected to try to ensure a good alignment. Making large changes
#to these WILL change the results, possibly substantially.

DEFAULT_QUERY_GAP_PENALTY = -25
DEFAULT_TEMPLATE_GAP_PENALTY = -25

#These are positions that are essentially mandatory. A sequence that deviates
#from one of these is an alignment issue or has a very large deletion.
heavy_conserved_positions = {22:"C", 36:"W", 92:"C", 103:"W", 104:"G", 106:"G"}
light_conserved_positions = {22:"C", 35:"W", 88:"C", 98:"F", 99:"G", 101:"G"}

#The penalty for inserting a gap at a highly conserved position.
HIGHLY_CONSERVED_GAP_PENALTY = -65
#The bonus for complying at a highly conserved position.
HIGHLY_CONSERVED_BONUS = 60


#These are positions where either A) a blank in the query sequence is very common or
#B) insertions are very common. We want special template and query gap penalties for
#these positions. One or two of these are chain dependent. The first value is
#the query gap column, the second is the template gap column.
heavy_special_positions = {82:[-11.0,-1.0], 6:[-11.0,-1.0],
        35:[-1.0,-1.0], 52:[-1.0,-1.0], 100:[-1.0,-1.0]}

light_special_positions = {10:[-1.0,-11.0], 66:[-11.0,-1.0],
        27:[-1.0,-1.0], 52:[-1.0,-1.0], 81:[-5.0,-1.0]}

# These are Kabat defined CDRs. For Kabat, unlike IMGT, heavy and light CDRs are different.
heavy_cdr_idx = list(range(31,36)) + list(range(50,66)) + list(range(95,103))
light_cdr_idx = list(range(24,35)) + list(range(50,57)) + list(range(89,98))
heavy_cdrs = {k:-11 for k in heavy_cdr_idx}
light_cdrs = {k:-11 for k in light_cdr_idx}
#heavy_cdrs = {31:-5.1, 32:-4.6, 33:-4.1, 34:-3.6, 35:-1,
#        # CDR 2
#        50:-5.1, 51:-4.6, 52:-1, 53:-11, 54:-11, 55:-11, 56:-11, 57:-11, 58:-11,
#                59:-11, 60:-11, 61:-11, 62:-11, 63:-11, 64:-11, 65:-11,
#        #CDR 3
#        95:-5.1, 96:-4.6, 97:-4.1, 98:-3.6, 99:-3.1, 100:-1,
#            101:-11, 102:-11}
#light_cdrs = {31:-5.1, 32:-4.6, 33:-4.1, 34:-3.6, 35:-1,
#        # CDR 2
#        50:-5.1, 51:-4.6, 52:-1, 53:-11, 54:-11, 55:-11, 56:-11, 57:-11, 58:-11,
#                59:-11, 60:-11, 61:-11, 62:-11, 63:-11, 64:-11, 65:-11,
#        #CDR 3
#        95:-5.1, 96:-4.6, 97:-4.1, 98:-3.6, 99:-3.1, 100:-1,
#            101:-11, 102:-11}
