# Only change if building a new version.
__version__ = "0.3.9"

import importlib.util

from .numbering_tools import SingleChainAnnotator, PairedChainAnnotator
from .numbering_tools.dna_sequence_translation import DNASeqTranslator
from .vj_tools.vj_gene_assignment import VJGeneTool
from .scoring_tools import SequenceScoringTool, HumanizationTool
from antpack.antpack_cpp_ext import LiabilitySearchTool
from .antpack_cli import run_cli_interface
from .antpack_license import run_license_key_setter

pyside = importlib.util.find_spec("PySide6")
if pyside is not None:
    from .antpack_gui.seq_viewer import run_seq_viewer
