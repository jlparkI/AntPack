#Used for determining the version if running setup.
#Only change if building a new version.
__version__ = "0.3"

from .numbering_tools import SingleChainAnnotator, PairedChainAnnotator, PyCTermFinder
#from .numbering_tools.consensus_update_tools import build_consensus_files
from .vj_tools.vj_gene_assignment import VJGeneTool
from .scoring_tools import SequenceScoringTool, HumanizationTool, LiabilitySearchTool
from .antpack_cli import run_cli_interface
