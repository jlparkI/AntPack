"""Simple utility for displaying message boxes of a requested type."""
from PySide6.QtWidgets import QMessageBox


def display_message_box(title, message, parent):
    """Displays a message box."""
    mbox = QMessageBox(parent=parent)
    mbox.setWindowTitle(title)
    mbox.setText(message)
    mbox.exec()
