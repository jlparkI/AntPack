"""Launches the application event loop; entry point."""
from PySide6 import QtWidgets
from qt_material import apply_stylesheet
from .gui_src.main_window import MainWindow



def run_seq_viewer():
    """Launches the sequence viewer application."""
    resp_app = QtWidgets.QApplication()
    apply_stylesheet(resp_app, theme='light_blue.xml')

    main_window = MainWindow()
    main_window.show()
    resp_app.exec()
