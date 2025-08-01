"""Contains utilities useful for loading fitted models."""
import os
import numpy as np
from ..clustering_tools.em_categorical_mixture import EMCategoricalMixture
from ..scoring_tools.scoring_constants import allowed_imgt_pos as ahip


def load_model(start_dir, chain_type="heavy",
        species="human", max_threads=2):
    """Loads the model for the specified species."""
    current_dir = os.getcwd()
    os.chdir(os.path.join(start_dir, "model_data"))

    if chain_type == "heavy":
        allowed_positions = ahip.heavy_allowed_positions
    elif chain_type == "light":
        allowed_positions = ahip.light_allowed_positions
    else:
        raise ValueError("Unrecognized chain type supplied.")

    try:
        mu = np.load(f"{species}_{chain_type}_mu.npy").astype(np.float64)
        mixweights = np.load(f"{species}_{chain_type}_mixweights.npy").astype(np.float64)
    except Exception as exc:
        raise ValueError("Final run not present in expected location.") from exc

    model = EMCategoricalMixture(mixweights.shape[0],
            numbering=allowed_positions,
            region="all", max_threads=max_threads)

    model.load_params(mu, mixweights)

    os.chdir(current_dir)
    return model
