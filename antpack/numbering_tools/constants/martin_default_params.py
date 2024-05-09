"""Defines important defaults for the Martin numbering scheme."""

#The expected length of the heavy and light chains
NUM_HEAVY = 113
NUM_LIGHT = 107

############################The following are alignment parameters that
#were selected to try to ensure a good alignment. Making large changes
#to these WILL change the results, possibly substantially.

DEFAULT_C_TERMINAL_QUERY_GAP_PENALTY = -1
DEFAULT_TERMINAL_TEMPLATE_GAP_PENALTY = -1

#These are positions that are essentially mandatory. A sequence that deviates
#from one of these is an alignment issue or has a very large deletion.
heavy_conserved_positions = {22:"C", 36:"W", 92:"C", 103:"W", 104:"G", 106:"G"}
light_conserved_positions = {23:"C", 35:"W", 88:"C", 98:"F", 99:"G", 101:"G"}

#These are positions where either A) a blank in the query sequence is very common or
#B) insertions are very common. We want special template and query gap penalties for
#these positions. One or two of these are chain dependent. The first value is
#the query gap column, the second is the template gap column.
heavy_special_positions = {8:[-25,-1.0], 31:[-1.0,-1.0],
        40:[-11.4,-25], 41:[-11.3,-25], 42:[-11.2,-25], 43:[-11.1,-25], 44:[-11,-25],
        52:[-1.0,-1.0], 72:[-25,-1.0], 100:[-1.0,-1.0]}

light_special_positions = {10:[-1,-11.0], 68:[-1.0,-1.0],
        #CDR 1
        30:[-11.1,-1.0],
        #CDR 2, 3
        52:[-1.0,-1.0], 95:[-1.0,-1.0]}

# These exclude special positions. For Martin, unlike IMGT, heavy and light CDRs are different.
heavy_cdrs = {26:-11.5, 27:-11.4, 28:-11.3, 29:-11.2, 30:-11.1,
        # CDR 2
        50:-11.2, 51:-11.1, 53:-11.3, 54:-11.4, 55:-11.5, 56:-11.6,
        #CDR 3
        95:-11.5, 96:-11.4, 97:-11.3, 98:-11.2, 99:-11.1,
            101:-11.6, 102:-11.7}
light_cdrs = {26:-11.5, 27:-11.4, 28:-11.3, 29:-11.2, 31:-11.0,
        # CDR 2
        50:-11.2, 51:-11.1, 53:-11, 54:-11.15, -55:-11.25, 56:-11.35,
        #CDR 3
        89:-11.5, 90:-11.4, 91:-11.3, 92:-11.2, 93:-11.1, 94:-11, 96:-11.05, 97:-11.15}
