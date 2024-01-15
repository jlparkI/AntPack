"""Contains default params shared by all schemes. Currently this comprises
only gap penalties for the N-terminal. This region is tricky because there
are frequenly large N-terminal deletions -- and also sometimes 
deletions located 10 AAs or so in from the N-terminal. We use
'staircase' gap penalties to encourage the aligner to handle this
sensibly. In any case where positions present in this gap penalty list
are also present in a scheme-specific constants file, the scheme-specific
file takes precedence."""

n_terminal_gap_positions = {}

for i in range(1,10):
    n_terminal_gap_positions[i] = [-1 - i, -25]

for i in range(10,18):
    n_terminal_gap_positions[i] = [-11 - (i-10) * 0.1, -25]


############################The following are alignment parameters that
#were selected to try to ensure a good alignment. Making large changes
#to these WILL change the results, possibly substantially.

DEFAULT_QUERY_GAP_PENALTY = -25
DEFAULT_TEMPLATE_GAP_PENALTY = -25

#The penalty for inserting a gap at a highly conserved position.
HIGHLY_CONSERVED_GAP_PENALTY = -65
#The bonus for complying at a highly conserved position.
HIGHLY_CONSERVED_BONUS = 65
