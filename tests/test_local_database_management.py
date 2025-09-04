"""Test database management & search.
Use sqlite for testing and prototyping."""
import os
import sqlite3
import unittest
from antpack import build_database_from_fasta, SingleChainAnnotator
from antpack.utilities import read_fasta
import numpy as np



class TestLocalDBManagement(unittest.TestCase):



    def test_local_db_search(self):
        """Check that searches yield expected results."""




if __name__ == "__main__":
    unittest.main()
