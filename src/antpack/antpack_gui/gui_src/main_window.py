"""Describes the MainWindow class which forms the core of the application."""
import os
from PySide6 import QtWidgets
from PySide6.QtWidgets import QApplication, QMainWindow, QTabWidget, QWidget, QVBoxLayout, QMenuBar, QFrame
from PySide6.QtWidgets import QTableWidget, QMessageBox, QFileDialog, QTableWidgetItem, QScrollArea
from PySide6.QtWidgets import QPushButton, QHBoxLayout, QStackedWidget, QSizePolicy, QLabel, QComboBox
from PySide6.QtCore import QSize, Qt
from PySide6.QtGui import QIcon, QAction, QPixmap, QColor, QBrush

from antpack import SingleChainAnnotator, PairedChainAnnotator, SequenceScoringTool
from antpack import VJGeneTool, LiabilitySearchTool

from .data_processing.selected_seq_processing import process_for_vj_comparison
from .data_processing.selected_seq_processing import MultiSequenceData, VJComparisonData
from .custom_widgets.sidebar import SideMenuWidget
from .dialogs.add_sequence_dialog import AddSequenceDialog
from .dialogs.vj_comp_dialog import VJComparisonDialog


class MainWindow(QMainWindow):
    """Main window for the application. Stores application state,
    launches sub-windows and other processes."""

    def __init__(self):
        super().__init__()

        self.setWindowTitle("AntPack")
        fpath = os.path.dirname(os.path.abspath(__file__))
        self.setWindowIcon(QIcon(os.path.join(fpath, "..",
            "..", "resources", "beaker.png")))
        self.setMinimumSize(QSize(800,400))

        central_widget = QWidget()
        self.setCentralWidget(central_widget)

        self.selected_seqs = None

        self.scheme = "imgt"
        self.sc_annotator = SingleChainAnnotator(scheme="imgt")
        self.pc_annotator = PairedChainAnnotator(scheme="imgt")
        self.vj_tool = VJGeneTool(scheme="imgt")
        self.liability_tool = LiabilitySearchTool()
        self.scoring_tool = SequenceScoringTool()

        # A default pixmap we will use for pieces which have not been
        # implemented yet.
        default_pixmap = QPixmap(os.path.join(fpath, "..", "..", "resources",
            "not_available.jpg"))

        # Add the special sidebar.
        sidebar_layout = QVBoxLayout()
        sidebar_contents = [
                ("Sequence Review", os.path.join(fpath, "..", "..",
                    "resources", "pipette--plus.png")),
                ("Dataset Review", os.path.join(fpath, "..", "..",
                    "resources", "flask--plus.png")),
                ("Process New\nDataset", os.path.join(fpath, "..", "..",
                    "resources", "database--plus.png"))
                ]
        self.sidebar = SideMenuWidget(sidebar_contents)
        self.sidebar.setSizePolicy(QSizePolicy.Fixed, QSizePolicy.Expanding)
        self.sidebar.listview.selectionModel().selectionChanged.connect(self.tab_select)
        sidebar_layout.addWidget(self.sidebar)

        # Next, create the simple sequence viewer. This will have two tabs,
        # one of which is for sequence comparison and one of which is
        # for plotting.
        seq_view_tabs = QTabWidget()
        tab1 = QWidget()
        seq_view_tab1_layout = QVBoxLayout()
        tab1.setLayout(seq_view_tab1_layout)

        top_box_layout = QHBoxLayout()
        add_seq_button = QPushButton("Add new sequence")
        add_seq_button.clicked.connect(self.seqview_add_sequence)
        top_box_layout.addWidget(add_seq_button)

        vj_compare = QPushButton("Compare single sequence\nwith VJ genes")
        vj_compare.clicked.connect(self.seqview_compare_with_vj)
        top_box_layout.addWidget(vj_compare)

        clear_seqs = QPushButton("Clear sequences")
        clear_seqs.clicked.connect(self.clear_sequences)
        top_box_layout.addWidget(clear_seqs)
        top_box_layout.addStretch(1)

        seq_view_tab1_layout.addLayout(top_box_layout)
        seq_view_tabs.addTab(tab1, "Sequence View")

        heavy_light_tabs = QTabWidget()
        seq_view_tab1_layout.addWidget(heavy_light_tabs)
        seq_view_tab1_layout.addStretch(1)

        self.heavy_chain_view = QTableWidget()
        self.light_chain_view = QTableWidget()
        self.heavy_chain_view.setEditTriggers(QTableWidget.NoEditTriggers)
        self.light_chain_view.setEditTriggers(QTableWidget.NoEditTriggers)
        self.heavy_chain_view.setSizeAdjustPolicy(QtWidgets.QAbstractScrollArea.AdjustToContents)
        self.light_chain_view.setSizeAdjustPolicy(QtWidgets.QAbstractScrollArea.AdjustToContents)


        heavy_light_tabs.addTab(self.heavy_chain_view, "Heavy")
        heavy_light_tabs.addTab(self.light_chain_view, "Light")
        self.update_seq_comparison_tabs()

        self.stacked_widget = QStackedWidget(self)
        self.stacked_widget.addWidget(seq_view_tabs)
        qp1_label = QLabel()
        qp2_label = QLabel()
        qp1_label.setPixmap(default_pixmap)
        qp2_label.setPixmap(default_pixmap)
        self.stacked_widget.addWidget(qp1_label)
        self.stacked_widget.addWidget(qp2_label)

        # Add the main widgets in each tab to the appropriate tab...
        main_layout = QHBoxLayout()
        main_layout.addLayout(sidebar_layout)
        main_layout.addWidget(self.stacked_widget)
        central_widget.setLayout(main_layout)

        # Next add menu options and the menu bar.
        menu_bar = QMenuBar()
        file_menu = menu_bar.addMenu("File")

        num_scheme_select = file_menu.addMenu("Select numbering scheme")
        self.imgt_select = QAction("IMGT", checkable=True)
        self.imgt_select.triggered.connect(self.set_imgt_scheme)
        self.imgt_select.setChecked(1)
        num_scheme_select.addAction(self.imgt_select)

        self.kabat_select = QAction("Kabat", checkable=True)
        self.kabat_select.triggered.connect(self.set_kabat_scheme)
        num_scheme_select.addAction(self.kabat_select)
        self.martin_select = QAction("Martin", checkable=True)
        self.martin_select.triggered.connect(self.set_martin_scheme)
        num_scheme_select.addAction(self.martin_select)
        self.aho_select = QAction("Aho", checkable=True)
        self.aho_select.triggered.connect(self.set_aho_scheme)
        num_scheme_select.addAction(self.aho_select)

        quit_command = QAction("Quit", self)
        quit_command.triggered.connect(self.quit_event)
        quit_command.setShortcut('Ctrl+Q')
        file_menu.addAction(quit_command)

        self.setMenuBar(menu_bar)



    def quit_event(self):
        """Callback when user selects quit under file."""
        self.close()


    def tab_select(self):
        """Callback when user selects a specific tab."""
        selected_idx = [idx.row() for idx in self.sidebar.listview.selectedIndexes()]
        if len(selected_idx) == 1 and selected_idx[0] < 3:
            self.stacked_widget.setCurrentIndex(selected_idx[0])
        else:
            self.stacked_widget.setCurrentIndex(0)


    def seqview_add_sequence(self):
        """Open the add sequence dialog box and if user does not
        cancel, add a sequence. Either way, update the seq compare
        and the properties windows."""
        if self.selected_seqs is not None:
            if isinstance(self.selected_seqs, VJComparisonData):
                _ = QMessageBox.warning(self, "Sequence not processed",
                    "You are currently viewing a VJ gene sequence comparison. "
                    "You must clear this before using the Add Sequence "
                    "option.", QMessageBox.Ok)
                return

        dialog = AddSequenceDialog(self)
        if dialog.exec():
            seq, seq_type, pid_thresh = dialog.get_options()
            if seq is None:
                _ = QMessageBox.warning(self, "Sequence not processed",
                        "The sequence and/or options you "
                        "entered were not valid.", QMessageBox.Ok)
                self.update_seq_comparison_tabs()
                return

            if self.selected_seqs is None:
                self.selected_seqs = MultiSequenceData()
            err = ""
            seq_name = str(self.selected_seqs.get_num_light())
            try:
                err = self.selected_seqs.add_selected_sequence(seq, self.sc_annotator,
                        self.pc_annotator, seq_type, pid_thresh, seq_name)
            except RuntimeError as e:
                print(e)
                # AntPack actually should not generate exceptions for unrecognized
                # characters -- it will return an error message -- so exceptions
                # here are unexpected and very unlikely. Just in case, though,
                # handle with the default error message.

            if err is not None:
                _ = QMessageBox.warning(self, "Sequence not processed",
                        "The sequence you entered was rejected "
                    "by AntPack's numbering tool. Common reasons include: "
                    "lowercase letters, unrecognized amino acids (X and the "
                    "typical 20 AAs are allowed but not gaps) and "
                    "punctuation.", QMessageBox.Ok)

        self.update_seq_comparison_tabs()




    def seqview_compare_with_vj(self):
        """Open a dialog box for settings for comparing a
        sequence with the most similar VJ genes."""
        if self.selected_seqs is not None:
            _ = QMessageBox.warning(self, "Sequence not processed",
                    "Before finding VJ genes for a sequence and "
                    "comparing it with those genes, you must "
                    "first clear the sequnce viewer.", QMessageBox.Ok)
            return
        if self.scheme != "imgt":
            _ = QMessageBox.warning(self, "Sequence not processed",
                    "Alignments with germline genes in the GUI tool "
                    "are currently supported only using IMGT numbering.",
                    QMessageBox.Ok)
            return


        dialog = VJComparisonDialog(self)
        if dialog.exec():
            seq, seq_type, pid_thresh, species = dialog.get_options()
            if seq is None:
                _ = QMessageBox.warning(self, "Sequence not processed",
                        "The sequence and/or options you "
                        "entered were not valid.", QMessageBox.Ok)
                self.update_seq_comparison_tabs()
                return

            try:
                self.selected_seqs = process_for_vj_comparison(seq, seq_type,
                    self.sc_annotator, self.pc_annotator,
                    self.vj_tool, pid_thresh, species)
            except RuntimeError as e:
                print(e)
                # AntPack actually should not generate exceptions for unrecognized
                # characters -- it will return an error message -- so exceptions
                # here are unexpected and very unlikely. Just in case, though,
                # handle with the default error message.
                self.selected_seqs = None

            if self.selected_seqs is None:
                _ = QMessageBox.warning(self, "Sequence not processed",
                        "The sequence you entered was rejected "
                    "by AntPack's numbering tool. Common reasons include: "
                    "lowercase letters, unrecognized amino acids (X and the "
                    "typical 20 AAs are allowed but not gaps) and "
                    "punctuation.", QMessageBox.Ok)

        self.update_seq_comparison_tabs()




    def update_seq_comparison_tabs(self):
        """Updates the spreadsheet rows for comparing sequences
        in the sequence viewer."""
        self.heavy_chain_view.clear()
        self.light_chain_view.clear()

        if self.selected_seqs is None:
            self.heavy_chain_view.setColumnCount(30)
            self.heavy_chain_view.setRowCount(2)
            self.light_chain_view.setColumnCount(30)
            self.light_chain_view.setRowCount(2)
            self.heavy_chain_view.horizontalHeader().setStyleSheet("""
                            QHeaderView::section {padding-left: 20px; border: 0px;
                                    padding-right: 5px}""")
            self.heavy_chain_view.horizontalHeader().setStyleSheet("""
                            QHeaderView::section {padding-left: 20px; border: 0px;
                                    padding-right: 5px}""")
            self.heavy_chain_view.resizeColumnsToContents()
            self.light_chain_view.resizeColumnsToContents()
            return

        # Light chain

        if self.selected_seqs.get_num_light() == 0:
            self.light_chain_view.setColumnCount(30)
            self.light_chain_view.setRowCount(2)
        else:
            self.light_chain_view.setRowCount(self.selected_seqs.get_num_light())
            light_nmbr = self.selected_seqs.get_light_numbering()
            light_labels = self.sc_annotator.assign_cdr_labels(light_nmbr, "L")

            self.light_chain_view.setColumnCount(len(light_nmbr) + 3)
            self.light_chain_view.setHorizontalHeaderItem(0, QTableWidgetItem("Description"))
            self.light_chain_view.setHorizontalHeaderItem(1, QTableWidgetItem("Percent\nidentity"))
            self.light_chain_view.setHorizontalHeaderItem(2,
                    QTableWidgetItem("Error message\n(if any)"))
            for i, (nmbr, label) in enumerate(zip(light_nmbr, light_labels)):
                if label.startswith("cdr"):
                    widget = QTableWidgetItem(nmbr + "\n__")
                else:
                    widget = QTableWidgetItem(nmbr)
                self.light_chain_view.setHorizontalHeaderItem(i+3, widget)

            self.light_chain_view.horizontalHeader().setStyleSheet("""
                            QHeaderView::section {padding-left: 0px; border: 0px;
                                    padding-right: 0px}""")

            light_data = self.selected_seqs.get_light_data()
            for i in range(len(light_data[0])):
                self.light_chain_view.setItem(i, 0, QTableWidgetItem(light_data[3][i]))
                self.light_chain_view.setItem(i, 1, QTableWidgetItem(light_data[1][i]))
                self.light_chain_view.setItem(i, 2, QTableWidgetItem(light_data[2][i]))
                for j, letter in enumerate(light_data[0][i]):
                    self.light_chain_view.setItem(i, j+3, QTableWidgetItem(letter))

        # Heavy chain

        if self.selected_seqs.get_num_heavy() == 0:
            self.heavy_chain_view.setColumnCount(30)
            self.heavy_chain_view.setRowCount(2)
        else:
            self.heavy_chain_view.setRowCount(self.selected_seqs.get_num_heavy())
            heavy_nmbr = self.selected_seqs.get_heavy_numbering()
            heavy_labels = self.sc_annotator.assign_cdr_labels(heavy_nmbr, "H")

            self.heavy_chain_view.setColumnCount(len(heavy_nmbr) + 3)
            self.heavy_chain_view.setHorizontalHeaderItem(0, QTableWidgetItem("Description"))
            self.heavy_chain_view.setHorizontalHeaderItem(1, QTableWidgetItem("Percent\nidentity"))
            self.heavy_chain_view.setHorizontalHeaderItem(2,
                    QTableWidgetItem("Error message\n(if any)"))
            for i, (nmbr, label) in enumerate(zip(heavy_nmbr, heavy_labels)):
                if label.startswith("cdr"):
                    widget = QTableWidgetItem(nmbr + "\n__")
                else:
                    widget = QTableWidgetItem(nmbr)
                self.heavy_chain_view.setHorizontalHeaderItem(i+3, widget)

            self.heavy_chain_view.horizontalHeader().setStyleSheet("""
                            QHeaderView::section {padding-left: 0px; border: 0px;
                                    padding-right: 0px}""")

            heavy_data = self.selected_seqs.get_heavy_data()
            for i in range(len(heavy_data[0])):
                self.heavy_chain_view.setItem(i, 0, QTableWidgetItem(heavy_data[3][i]))
                self.heavy_chain_view.setItem(i, 1, QTableWidgetItem(heavy_data[1][i]))
                self.heavy_chain_view.setItem(i, 2, QTableWidgetItem(heavy_data[2][i]))
                for j, letter in enumerate(heavy_data[0][i]):
                    self.heavy_chain_view.setItem(i, j+3, QTableWidgetItem(letter))

        self.heavy_chain_view.resizeColumnsToContents()
        self.light_chain_view.resizeColumnsToContents()



    def clear_sequences(self):
        """Clears the currently entered sequences."""
        self.selected_seqs = None
        self.update_seq_comparison_tabs()


    def set_imgt_scheme(self):
        """Sets numbering scheme to IMGT."""
        self.scheme = "imgt"
        self.sc_annotator = SingleChainAnnotator(scheme=self.scheme)
        self.pc_annotator = PairedChainAnnotator(scheme=self.scheme)
        self.imgt_select.setChecked(1)
        self.martin_select.setChecked(0)
        self.kabat_select.setChecked(0)
        self.aho_select.setChecked(0)
        self.clear_sequences()


    def set_martin_scheme(self):
        """Sets numbering scheme to martin."""
        self.scheme = "martin"
        self.sc_annotator = SingleChainAnnotator(scheme=self.scheme)
        self.pc_annotator = PairedChainAnnotator(scheme=self.scheme)
        self.martin_select.setChecked(1)
        self.imgt_select.setChecked(0)
        self.kabat_select.setChecked(0)
        self.aho_select.setChecked(0)
        self.clear_sequences()


    def set_kabat_scheme(self):
        """Sets numbering scheme to Kabat."""
        self.scheme = "kabat"
        self.sc_annotator = SingleChainAnnotator(scheme=self.scheme)
        self.pc_annotator = PairedChainAnnotator(scheme=self.scheme)
        self.kabat_select.setChecked(1)
        self.martin_select.setChecked(0)
        self.imgt_select.setChecked(0)
        self.aho_select.setChecked(0)
        self.clear_sequences()


    def set_aho_scheme(self):
        """Sets numbering scheme to aho."""
        self.scheme = "aho"
        self.sc_annotator = SingleChainAnnotator(scheme=self.scheme)
        self.pc_annotator = PairedChainAnnotator(scheme=self.scheme)
        self.kabat_select.setChecked(0)
        self.martin_select.setChecked(0)
        self.imgt_select.setChecked(0)
        self.aho_select.setChecked(1)
        self.clear_sequences()
