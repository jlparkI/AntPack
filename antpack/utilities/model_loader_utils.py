"""Contains utilities useful for loading fitted models."""
import os
import pickle
import numpy as np
from ..categorical_mix import CategoricalMixture
from ..constants.scoring_constants import allowed_imgt_pos as ahip


def load_model(start_dir, chain_type="heavy", species="human"):
    """Loads the model for the specified species."""
    os.chdir(os.path.join(start_dir, "model_data"))
    fname = f"{species}_{chain_type}_model.pk"
    if fname not in os.listdir():
        raise ValueError("Final run not present in expected location.")

    if chain_type == "heavy":
        sequence_length = len(ahip.heavy_allowed_positions)
    elif chain_type == "light":
        sequence_length = len(ahip.light_allowed_positions)
    else:
        raise ValueError("Unrecognized chain type supplied.")

    model = CategoricalMixture(200, num_possible_items = 21,
                         sequence_length = sequence_length)
    with open(f"{species}_{chain_type}_model.pk", "rb") as fhandle:
        model_data = pickle.load(fhandle)

    model.n_components = model_data["mixweights"].shape[0]
    model.load_params(model_data["mu"], model_data["mixweights"])

    return model
