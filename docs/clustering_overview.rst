Clustering antibody sequences in AntPack
===============================================

For clustering small datasets (a few dozen to a few thousand
sequences), AntPack makes it easy to construct a distance matrix
using Hamming distance for any specified subregion of your
sequences (for example, the framework, the CDRs, or a specific
framework or CDR region) or the full sequence if needed. To do
this, use the ``build_distance_matrix`` method of the
``SingleChainAnnotator`` and ``PairedChainAnnotator`` tools
used for numbering:

.. autoclass:: antpack.SingleChainAnnotator
   :special-members: __init__
   :members: build_distance_matrix

Because distance matrices have n^2 scaling in size and construction
time with dataset size, this method is recommended only for datasets
up to a few thousand sequences in size. Once you've built a distance
matrix, you can supply it to a variety of Scipy and scikit-learn
functions and easily cluster and visualize it using a variety of methods.
See the Clustering examples on the main page of the docs to see how to
do this.

For large datasets, AntPack currently offers just one *highly* scalable
option (other options coming soon). The ``EMCategoricalMixture`` is
designed to use multithreading and run with datasets too large to load
to memory; it can easily cluster datasets ranging from a few hundred to
tens of millions in size. This model is a mixture model which assigns
a probability to each amino acid at each position in each cluster. As
such, it is a probabilistic model that can calculate the probability
that a new sequence could have come from its training data or the
probability of a specific mutation given its training set. As with
distance matrix construction, you can cluster using the whole
sequence or any selected subregion.

Choosing the number of clusters for the EMCategoricalMixture can be a
little tricky. Currently it provides [BIC](https://en.wikipedia.org/wiki/Bayesian_information_criterion)
and [AIC](https://en.wikipedia.org/wiki/Akaike_information_criterion) and you can
use these to select the number of clusters by finding the number that gives the
lowest BIC and AIC. (For an example of how to do this -- and how to use the
EMCategoricalMixture more generally -- see the Clustering examples on the
main page of the docs.) If you have a small dataset and are clustering
however AIC and BIC will tend to select a number of clusters that is too small
(e.g. 1!) Also note that it's ok to use more clusters than needed because
the EMCategoricalMixture will eliminate empty clusters during fitting, so
if the number of clusters is larger than needed it can kill off some unneeded
ones. In future versions we will likely add a procedure to auto-select the
number of clusters so that no manual selection is required.

``EMCategoricalMixture`` is designed to use multithreading to process multiple
batches of data in parallel during fitting; selecting a ``max_threads`` > 1
will automatically enable multithreading and is highly recommended if you're
fitting a large dataset. If you have a large dataset (e.g. millions of sequences)
in a fasta file or gzipped fasta file and don't want to load it to memory,
``EMCategoricalMixture`` can encode it to a set of temporary files in a location
you specify by calling ``encode_fasta_file``. You supply the list of these temporary
files to the ``fit`` function and ``EMCategoricalMixture`` will take care of the
rest -- just remember to delete the temporary files once you're done.

.. autoclass:: antpack.EMCategoricalMixture
   :special-members: __init__
   :members: get_model_parameters, load_params, get_model_specs, fit, BIC, AIC, predict, predict_proba, score, encode_fasta
