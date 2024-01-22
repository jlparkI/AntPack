"""Contains constants used across multiple procedures."""

#The number of positions in the input sequence.
SEQUENCE_LENGTH = 173


#The list of aas; the position in the list determines the
#number assigned for one-hot encoding.
aa_list = ["A", "C", "D", "E", "F", "G", "H", "I", "K",
        "L", "M", "N", "P", "Q", "R", "S", "T", "V",
        "W", "Y", "-"]

#The allowed amino acids.
allowed_aas = set(aa_list)

#The minimum percent identity allowed. This is pretty generous...most should
#be much higher.
MIN_PERCENT_IDENTITY = 0.8

#The number of positions required to be non-blank. Fairly generous.
#MIN_HEAVY_NOT_BLANK = 110
MIN_HEAVY_NOT_BLANK = 105
MIN_LIGHT_NOT_BLANK = 105


#The number of IMGT positions expected in each sequence.
NUM_IMGT_POSITIONS = 235

#The number of sequences to save in each encoding file.
CHUNK_SIZE = 10000

#The rescaling parameters for heavy and light chain scores
#when combining them. Based on training set.
HEAVY_MEDIAN_SCORE = -59.5617
LIGHT_MEDIAN_SCORE = -33.22445
HEAVY_IQR_SCORE = 42.022
LIGHT_IQR_SCORE = 24.1216
