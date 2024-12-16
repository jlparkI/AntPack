"""Dialog box for selecting the sequence type for an input file."""
from PySide6.QtWidgets import QDialog, QRadioButton, QVBoxLayout, QLabel, QDialogButtonBox



class SequenceTypeSelector(QDialog):
    """Selects the sequence type (paired or not)."""

    def __init__(self, parent):
        super().__init__(parent)

        self.setWindowTitle("Choose sequence type")

        self.layout = QVBoxLayout()
        q_button = QDialogButtonBox.Ok

        self.button_box = QDialogButtonBox(q_button)
        self.button_box.accepted.connect(self.accept)

        self.paired = QRadioButton("Paired (sequences contain heavy and\nlight chains)", self)
        self.single = QRadioButton("Single (sequences contain heavy or\nlight chains only)", self)

        self.paired.setChecked(1)

        self.layout = QVBoxLayout()
        message = QLabel("Please choose the type of the sequences in "
                "your input data.")
        self.layout.addWidget(message)
        self.layout.addWidget(self.paired)
        self.layout.addWidget(self.single)
        self.layout.addWidget(self.button_box)
        self.setLayout(self.layout)

        self.exec()
