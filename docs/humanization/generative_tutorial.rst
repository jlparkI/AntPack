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
notebooks on the main page to see how to do this.

It is important to realize that sampling from a specified cluster
can generate low-probability sequences, i.e. sequences
that are unlikely to be human, in the same way that sampling from
a normal distribution can occasionally generate outliers. It may
be useful to score any sequences generated in this fashion using
the scoring tool to ensure they are highly human. Alternatively,
you may want to use the highest-probability amino acid at each
position in a selected cluster *except* for the CDRs and sample
from the distribution for those regions specifically.
