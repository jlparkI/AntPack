# Only change if building a new version.
__version__ = "0.4"

import importlib.util

from .numbering_tools import SingleChainAnnotator, PairedChainAnnotator
from .vj_tools.vj_gene_assignment import VJGeneTool
from .scoring_tools import SequenceScoringTool
from antpack.antpack_cpp_ext import (LiabilitySearchTool,
                SequenceTemplateAligner)
from .clustering_tools import EMCategoricalMixture
from .antpack_license import run_license_key_setter
from .database_tools import (build_database_from_fasta,
        LocalDBSearchTool,
        build_database_from_csv,
        build_tcr_database_from_csv)

pyside = importlib.util.find_spec("PySide6")
if pyside is not None:
    from .antpack_gui.seq_viewer import run_seq_viewer
