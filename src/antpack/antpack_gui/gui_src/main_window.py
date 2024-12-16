"""Describes the MainWindow class which forms the core of the application."""
import os
import qdarktheme
from PySide6 import QtGui
from PySide6.QtWidgets import QApplication, QMainWindow, QMenuBar, QTabWidget, QWidget, QVBoxLayout, QTableWidget, QMessageBox, QFileDialog, QTableWidgetItem, QScrollArea
from PySide6.QtCore import QSize
from PySide6.QtGui import QIcon, QAction
import pyqtgraph as pgp

from .data_loaders.dataset import Dataset
from .utils.message_box import display_message_box
from .dialogs.sequence_type import SequenceTypeSelector



class MainWindow(QMainWindow):
    """Main window for the application. Stores application state,
    launches sub-windows and other processes."""

    def __init__(self):
        super().__init__()

        self.theme = "dark"
        self.setWindowTitle("RESP / AntPack")
        fpath = os.path.dirname(os.path.abspath(__file__))
        self.setWindowIcon(QIcon(os.path.join(fpath, "..", "resources", "mAb.png")))
        self.setMinimumSize(QSize(600,400))


        # Any loaded data is stored as a Dataset object; the AntPack library
        # is used extensively in construction of the Dataset.
        self.dataset = None

        # Always store the user selected numbering scheme.
        self.scheme = "imgt"

        # There are several main tabs the user can access.
        self.tab_widget = QTabWidget()
        self.setCentralWidget(self.tab_widget)

        # Key elements of these tabs are the raw data table and the property
        # plot.
        self.property_graph = pgp.PlotWidget()
        self.property_graph.setMouseEnabled(False, False)
        self.update_property_plot()

        self.raw_data_table = QTableWidget()
        self.raw_data_table.setEditTriggers(QTableWidget.NoEditTriggers)
        self.reset_raw_data_table()
        selection = self.raw_data_table.selectionModel()
        selection.selectionChanged.connect(self.update_table_selection)

        self.raw_tab = QWidget()
        self.property_tab = QWidget()
        self.alignment_tab = QWidget()
        #self.cluster_tab = QWidget()

        # Add the main widgets in each tab to the appropriate tab...
        raw_tab_layout = QVBoxLayout(self.raw_tab)
        property_tab_layout = QVBoxLayout(self.property_tab)

        self.scroll = QScrollArea()
        raw_tab_layout.addWidget(self.raw_data_table)
        property_tab_layout.addWidget(self.scroll)
        self.scroll.setWidget(self.property_graph)
        self.scroll.setWidgetResizable(True)



        # ...now add the tabs to the tab widget.
        self.tab_widget.addTab(self.raw_tab, "Raw Data")
        self.tab_widget.addTab(self.property_tab, "Property Viewer")
        self.tab_widget.addTab(self.alignment_tab, "Alignment Viewer")

        # Next add menu options and the menu bar.
        theme_command = QAction("Switch color scheme", self)
        theme_command.triggered.connect(self.theme_event)
        quit_command = QAction("Quit", self)
        quit_command.triggered.connect(self.quit_event)
        quit_command.setShortcut('Ctrl+Q')

        menu_bar = QMenuBar()
        file_menu = menu_bar.addMenu("File")
        file_menu.addAction(theme_command)
        file_menu.addAction(quit_command)


        fasta_load = QAction("Load FASTA file", self)
        fasta_load.triggered.connect(self.fasta_load_event)
        csv_load = QAction("Load csv file", self)
        csv_load.triggered.connect(self.csv_load_event)

        data_menu = menu_bar.addMenu("Data")

        num_scheme_select = data_menu.addMenu("Select numbering scheme")
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

        data_menu.addAction(fasta_load)
        data_menu.addAction(csv_load)

        self.setMenuBar(menu_bar)


        # Next add tools specific to the property tab.
        self.selected_property = "Hydrophobicity"




    def quit_event(self):
        """Callback when user selects quit under file."""
        self.close()


    def fasta_load_event(self):
        """Prompts the user for a filepath to a fasta file,
        asks some questions about what kinds of sequences it
        contains, and loads the data."""
        self.reset_raw_data_table()
        self.dataset = None
        seqtype_sel = SequenceTypeSelector(self)

        if seqtype_sel.paired.isChecked():
            dtype = "paired"
        elif seqtype_sel.single.isChecked():
            dtype = "single"

        new_dataset = Dataset(dtype, self.scheme)
        file_name = QFileDialog.getOpenFileName(self,
                    "Open FASTA", filter="FASTA Files (*.fa *.fasta)")[0]

        err_message = new_dataset.load_fasta_file(file_name)

        if err_message is not None:
            display_message_box("Error", err_message, self)
            return

        self.dataset = new_dataset
        self.update_data_table()



    def csv_load_event(self):
        """Prompts the user for a filepath to a csv file,
        asks some questions about what kinds of sequences it
        contains, and loads the data."""
        self.reset_raw_data_table()
        self.dataset = None
        pass


    def reset_raw_data_table(self):
        """Resets the raw data table when loading data or
        initializing the application."""
        self.raw_data_table.clearContents()
        self.raw_data_table.setColumnCount(7)
        self.raw_data_table.setRowCount(20)
        self.raw_data_table.setHorizontalHeaderLabels(["ID", "Metadata",
            "Heavy chain\nsequence", "Light chain\nsequence",
            "Errors/\nissues", "Input file", "Cluster"])
        # Reset the property plot as well.
        self.update_property_plot()


    def update_data_table(self):
        """Updates the data table with the data from a newly-loaded dataset.
        Assumes that reset_raw_data_table was called recently so that
        num columns and column headers are already correct."""
        if self.dataset is None:
            return
        if self.dataset.get_num_seqs() == 0:
            return
        self.raw_data_table.setRowCount(self.dataset.get_num_seqs())

        for i, seq_data in enumerate(self.dataset.get_seq_data()):
            self.raw_data_table.setItem(i, 0, QTableWidgetItem(str(i)))
            self.raw_data_table.setItem(i, 1, QTableWidgetItem(seq_data.get_metadata()))
            self.raw_data_table.setItem(i, 2, QTableWidgetItem(seq_data.get_heavy_chain()))
            self.raw_data_table.setItem(i, 3, QTableWidgetItem(seq_data.get_light_chain()))
            self.raw_data_table.setItem(i, 4, QTableWidgetItem(seq_data.get_errors()))
            self.raw_data_table.setItem(i, 5, QTableWidgetItem(self.dataset.get_input_file()))



    def update_table_selection(self, selected, deselected):
        """Updates the properties tab if the selected sequence has changed."""
        if self.dataset is None:
            return
        selected_rows = self.raw_data_table.selectionModel().selectedRows()
        if len(selected_rows) != 1:
            self.update_property_plot()
            return

        self.update_property_plot(selected_rows[0].row())


    def update_property_plot(self, seqnum = None):
        """Updates the property plot tab. If seqnum or dataset is None,
        it resets it to a blank plot. If seqnum is not None, it pulls the
        corresponding sequence data from the dataset and plots it on
        the property tab."""
        self.property_graph.clear()
        self.property_graph.setDefaultPadding(0)

        #color = self.palette().color(QtGui.QPalette.Window)
        #self.property_graph.setBackground(color)
        if self.theme == "light":
            self.property_graph.setBackground("w")
        else:
            self.property_graph.setBackground("k")
            #pgp.setConfigOption("foreground", "gray")
        x = None

        if self.dataset is None or seqnum is None:
            return
        elif seqnum >= self.dataset.get_num_seqs():
            return

        idx, x, y, metadata = self.dataset.get_property(seqnum, self.selected_property)
        if len(x) == 0:
            return

        self.property_graph.setFixedWidth(len(x) * 30)

        bar_graph = pgp.BarGraphItem(x=idx, y1=y, width=1)
        self.property_graph.addItem(bar_graph)
        self.property_graph.setLabel("left", self.selected_property)
        self.property_graph.setLabel("bottom", metadata + "\n")
        self.property_graph.getPlotItem().getAxis('bottom').setTicks([[(idn, nmbr)
                for (idn, nmbr) in zip(idx, x) ]])
        #self.property_graph.getPlotItem().getAxis('left').setWidth(20)
        #self.property_graph.getPlotItem().getAxis('bottom').setHeight(20)



    def theme_event(self):
        """Switches between light and dark theme."""
        if self.theme == "dark":
            qdarktheme.setup_theme("light")
            self.theme = "light"
        else:
            qdarktheme.setup_theme("dark")
            self.theme = "dark"

        self.update_property_plot()


    def set_imgt_scheme(self):
        """Sets numbering scheme to IMGT."""
        self.scheme = "imgt"
        self.imgt_select.setChecked(1)
        self.martin_select.setChecked(0)
        self.kabat_select.setChecked(0)

    def set_martin_scheme(self):
        """Sets numbering scheme to martin."""
        self.scheme = "martin"
        self.martin_select.setChecked(1)
        self.imgt_select.setChecked(0)
        self.kabat_select.setChecked(0)

    def set_kabat_scheme(self):
        """Sets numbering scheme to Kabat."""
        self.scheme = "kabat"
        self.kabat_select.setChecked(1)
        self.martin_select.setChecked(0)
        self.imgt_select.setChecked(0)
