"""Launches the application event loop; entry point."""
from PySide6.QtWidgets import QApplication
import qdarktheme
from .gui_src.main_window import MainWindow


def run_seq_viewer():
    """Launches the sequence viewer application."""
    resp_app = QApplication()
    qdarktheme.setup_theme()
    main_window = MainWindow()
    main_window.show()
    resp_app.exec()
