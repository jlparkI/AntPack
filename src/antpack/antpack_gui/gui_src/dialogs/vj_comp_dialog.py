"""Dialog box to find VJ genes for an input sequence."""
from PySide6.QtWidgets import QDialog, QRadioButton, QVBoxLayout
from PySide6.QtWidgets import QLabel, QDialogButtonBox, QLineEdit, QComboBox
from PySide6.QtGui import QAction



class VJComparisonDialog(QDialog):
    """Selects the sequence type (paired or not)."""

    def __init__(self, parent):
        super().__init__(parent)

        self.setWindowTitle("Enter a sequence to see VJ genes:")
        self.setMinimumWidth(500)

        buttons = (
                    QDialogButtonBox.Ok
                    | QDialogButtonBox.Cancel
                )

        self.layout = QVBoxLayout()

        self.button_box = QDialogButtonBox(buttons)
        self.button_box.accepted.connect(self.accept)
        self.button_box.rejected.connect(self.reject)

        self.layout = QVBoxLayout()

        line_edit_label = QLabel("Enter your sequence here:")
        self.line_edit = QLineEdit()
        self.layout.addWidget(line_edit_label)
        self.layout.addWidget(self.line_edit)

        paired_select_label = QLabel("Please choose the type of sequence.")
        paired = QRadioButton("Paired (heavy and\nlight chains)", self)
        single = QRadioButton("Single (heavy or\nlight chain only)", self)
        unknown = QRadioButton("I don't know, decide for me")
        paired.setChecked(1)
        self.paired_radio_buttons = [paired, single, unknown]

        self.layout.addWidget(paired_select_label)
        self.layout.addWidget(paired)
        self.layout.addWidget(single)
        self.layout.addWidget(unknown)
        self.layout.addWidget(self.button_box)

        pid_edit_label = QLabel("Percent identity threshold to accept a chain as valid\n"
                    "(should be between 0 and 1):")
        self.pid_edit = QLineEdit("0.7")
        self.layout.addWidget(pid_edit_label)
        self.layout.addWidget(self.pid_edit)

        species_label = QLabel("Select species for VJ genes:")
        self.species_menu = QComboBox()
        self.species_list = ["human", "mouse", "unknown"]
        for species in self.species_list:
            self.species_menu.addItem(species)

        self.layout.addWidget(species_label)
        self.layout.addWidget(self.species_menu)

        self.setLayout(self.layout)




    def get_options(self):
        """Enables caller to retrieve the options that the user
        selected once the dialog has run."""
        seq_text = self.line_edit.text()
        if len(seq_text) == 0:
            return None, None, None, None

        species = self.species_list[self.species_menu.currentIndex()]
        if species not in ("human", "mouse", "unknown"):
            return None, None, None, None

        pid_thresh = self.pid_edit.text()
        try:
            pid_thresh = float(pid_thresh)
        except:
            return None, None, None, None

        if pid_thresh < 0 or pid_thresh > 1:
            return None, None, None, None

        button_choice = "paired"
        for button_choice, button in zip(("paired", "single", "unknown"),
                self.paired_radio_buttons):
            if button.isChecked():
                break
        return seq_text, button_choice, pid_thresh, species
