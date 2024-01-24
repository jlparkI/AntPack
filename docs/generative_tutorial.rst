Generating new antibody sequences
===============================================

The mixture model used by AntPack has over 1800 clusters
for heavy chains and 1300 for light. Each cluster is
itself a distribution that you can sample from. By
retrieving the cluster that's most similar to an
input sequence, you can generate new sequences that
contain specific features and use this approach to
build small synthetic libraries.

First, create a scoring tool:::

  from antpack import SequenceScoringTool
  scoring_tool = SequenceScoringTool()

Next, take your sequence of interest and find the closest
cluster(s) to it in the mixture models in the scoring tool.
For a large number of sequences, you can call ``batch_score_seqs``
(see the humanness evaluation page) with ``mode="assign"`` to
get the closest cluster number for each sequence, then
use ``retrieve_cluster`` to get the per-position, per-amino-acid
probabilities associated with those clusters. For a single
sequence, however, it's easiest to use ``get_closest_clusters``
as shown below.

.. autoclass:: antpack.SequenceScoringTool
   :members: get_closest_clusters, retrieve_cluster

Each cluster you retrieve is a distribution and you can easily
sample from it or use it to visualize which sections of your input
sequence are least human / most problematic. See the example
notebook on the main page to see how to do this.
