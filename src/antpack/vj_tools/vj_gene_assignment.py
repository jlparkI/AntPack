"""This module contains tools needed to assign the closest VJ germline
VJ gene (highest percent identity), using amino acid information
only. Bear in mind that tools which offer a probabilistic assigment
(e.g. IGOR) may be more informative."""
import os
import numpy as np
from ..antpack_license import get_license_key_info
from ..utilities.vj_utilities import load_vj_gene_consensus_db
from antpack.antpack_cpp_ext import VJMatchCounter



class VJGeneTool(VJMatchCounter):
    """Contains functionality needed to find the closest
    VJ genes for a given amino acid sequence and retrieve the
    sequence for those VJ genes. The date of the database
    used for these assignments is stored and can be retrieved
    if needed."""

    def __init__(self, scheme = "imgt"):
        """Class constructor.

        Args:
            scheme (str): One of 'aho, 'imgt', 'kabat' or 'martin'.
                Determines the numbering scheme that will be used
                when assigning vj genes based on alignments.
        """
        license_key, user_email = get_license_key_info()

        project_path = os.path.abspath(os.path.dirname(__file__))
        db_path = os.path.join(project_path, "consensus_data")
        vj_names, vj_seqs, retrieved_dates = load_vj_gene_consensus_db(os.getcwd(),
                db_path, "imgt")

        blosum_matrix = np.load(os.path.join(project_path, "..",
            "numbering_tools", "consensus_data", "mabs",
            "blosum_matrix.npy")).astype(np.float64)
        self.retrieved_dates = retrieved_dates
        ig_aligner_consensus_path = os.path.join(project_path, "..", "numbering_tools",
                "consensus_data")
        super().__init__(vj_names, vj_seqs, blosum_matrix,
                scheme, ig_aligner_consensus_path,
                license_key, user_email)



    def retrieve_db_dates(self):
        """Returns the dates when each VJ gene database
        used for this assignment was last updated by downloading
        from IMGT or OGRDB."""
        return self.retrieved_dates
