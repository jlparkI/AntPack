"""Launches the application event loop; entry point."""
import os
from PySide6 import QtWidgets
#from qt_material import apply_stylesheet
from .gui_src.main_window import MainWindow
from qt_material import apply_stylesheet


def run_seq_viewer():
    """Launches the sequence viewer application."""
    fpath = os.path.dirname(os.path.abspath(__file__))
    #icon_paths = {"app_icon":os.path.join(fpath, "..",
    #    "resources", "beaker.png")}
    #consensus_filepath = os.path.join(fpath, "..",
    #    "numbering_tools", "consensus_data")
    #run_antpack_gui(consensus_filepath,
    #    icon_paths)
    resp_app = QtWidgets.QApplication()
    apply_stylesheet(resp_app, theme='light_blue.xml')
    #sheet_to_style(resp_app, os.path.join(fpath, "main_stylesheet.css.template"))
    main_window = MainWindow()
    main_window.show()
    resp_app.exec()
