"""Launches the application event loop; entry point."""
import os
from antpack.antpack_gui_ext import run_antpack_gui



def run_seq_viewer():
    """Launches the sequence viewer application."""
    fpath = os.path.dirname(os.path.abspath(__file__))
    icon_paths = {"app_icon":os.path.join(fpath, "..",
        "resources", "beaker.png")}
    consensus_filepath = os.path.join(fpath, "..",
        "numbering_tools", "consensus_data")
    run_antpack_gui(consensus_filepath,
        icon_paths)
