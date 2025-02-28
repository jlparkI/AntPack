"""Contains utilities useful for compiling the raw sequence data
and encoding."""
from ..scoring_tools.scoring_constants import allowed_imgt_pos as ahip


def get_position_dict(chain_type="heavy"):
    """Loads the dictionary mapping IMGT numbers to positions."""
    if chain_type == "heavy":
        positions = ahip.heavy_allowed_positions
    elif chain_type == "light":
        positions = ahip.light_allowed_positions
    else:
        raise ValueError("Unrecognized chain type passed.")

    position_dict = {pos:i for (i, pos) in enumerate(positions)}
    return position_dict, len(position_dict)


def get_reverse_position_dict(chain_type="heavy"):
    """Sets up a dictionary mapping positions to IMGT numbers."""
    if chain_type == "heavy":
        positions = ahip.heavy_allowed_positions
    elif chain_type == "light":
        positions = ahip.light_allowed_positions
    else:
        raise ValueError("Unrecognized chain type passed.")

    return dict(enumerate(positions))
