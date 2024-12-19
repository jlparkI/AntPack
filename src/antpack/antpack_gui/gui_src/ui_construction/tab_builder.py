"""Contains some components of the UI construction (to avoid making
the main window object too cluttered)."""
from PySide6.QtWidgets import QTabWidget, QLineEdit, QMenuBar
from PySide6.QtWidgets import QWidget, QVBoxLayout, QHBoxLayout
from PySide6.QtWidgets import QTableWidget, QScrollArea
from PySide6.QtWidgets import QPushButton, QSizePolicy, QLabel
from PySide6.QtGui import QAction


def build_menubar(parent):
    menu_bar = QMenuBar(parent)
    file_menu = menu_bar.addMenu("File")

    num_scheme_select = file_menu.addMenu("Select numbering scheme")
    imgt_select = QAction("IMGT", parent)
    imgt_select.triggered.connect(parent.set_imgt_scheme)
    #imgt_select.setChecked(1)
    num_scheme_select.addAction(imgt_select)

    # For now only IMGT is supported, until we can update
    # VJ gene formatting.
    #kabat_select = QAction("Kabat", checkable=True)
    #kabat_select.triggered.connect(parent.set_kabat_scheme)
    #num_scheme_select.addAction(kabat_select)

    #martin_select = QAction("Martin", checkable=True)
    #martin_select.triggered.connect(parent.set_martin_scheme)
    #num_scheme_select.addAction(martin_select)

    quit_command = QAction("Quit", parent)
    quit_command.triggered.connect(parent.quit_event)
    quit_command.setShortcut('Ctrl+Q')
    file_menu.addAction(quit_command)
    return menu_bar
