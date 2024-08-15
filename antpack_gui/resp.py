"""Launches the application event loop; entry point."""
from PySide6.QtWidgets import QApplication
import qdarktheme
from gui_src.main_window import MainWindow


if __name__ == "__main__":
    resp_app = QApplication()
    qdarktheme.setup_theme()
    main_window = MainWindow()
    main_window.show()
    resp_app.exec()
