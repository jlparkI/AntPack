"""Tests fitting procedures for the categorical mixture."""
import os
import numpy as np
import timeit
import yep
#from categorical_mix import CategoricalMixture
from antpack.clustering_tools import EMCategoricalMixture




def load_test_data():
    """Loads some saved test data."""
    current_dir = os.getcwd()
    start_dir = os.path.abspath(os.path.dirname(__file__))
    data_path = os.path.join(start_dir, "..", "..", "tests", "test_data",
            "mixture_model_test_data", "encoded_test_data.npy")
    xd = np.load(data_path)
    xd = np.vstack([xd for i in range(100)])
    xf = 100 * [os.path.abspath(data_path)]
    print(f"There are {xd.shape[0]} datapoints.")
    os.chdir(current_dir)
    return xd, xf


def build_default_model(nclusters = 100, seq_length = 408,
        verbose=False, n_threads=10):
    """Builds an initial default model for the test set."""
    cat_mix = EMCategoricalMixture(n_components=nclusters,
            sequence_length=seq_length, max_threads=n_threads,
            verbose=verbose)
    #cat_mix = CategoricalMixture(n_components=nclusters,
    #        sequence_length=seq_length)
    return cat_mix



def fit_test(xfiles, base_model):
    base_model.fit(xfiles, n_restarts=1, max_iter=150, tol=1e-3,
            prune_after_fitting=False)




if __name__ == "__main__":
    #xdata, xfiles = load_test_data()
    #base_model = build_default_model(verbose=False)

    #yep.start('SPEED_PROFILING.prof')
    #base_model.fit(xfiles, max_iter = 150, n_restarts=3,
    #        prune_after_fitting = False)
    #yep.stop()

    timetest_setup = f"""
import numpy as np
from __main__ import build_default_model, fit_test, load_test_data
base_model = build_default_model(verbose=True)
xdata, xfiles = load_test_data()"""
    n_tests = 5
    time_taken = timeit.timeit("fit_test(xdata, base_model)",
            setup=timetest_setup, number=n_tests)
    print(f"Took {1e6 * time_taken / n_tests} us")
