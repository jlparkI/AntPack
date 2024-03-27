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
You'll use ``mode='assign'`` to do this.
Then get the parameters for that cluster.

.. autoclass:: antpack.SequenceScoringTool
   :members: score_seqs, retrieve_cluster

Each cluster you retrieve is a distribution and you can easily
sample from it or use it to visualize which sections of your input
sequence are least human / most problematic. See the example
notebook on the main page to see how to do this.
