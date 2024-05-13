"""This module contains tools needed to assign the closest VJ germline
VJ gene (highest percent identity), using amino acid information
only. Bear in mind that nucleotide assignments are likely more
reliable, and that tools which offer a probabilistic assigment
(e.g. IGOR) may be more informative."""
import os
import gzip


class VJGeneTool:
    """Contains functionality needed to find the closest
    VJ genes for a given amino acid sequence and retrieve the
    sequence for those VJ genes. The date of the database
    used for these assignments is stored and can be retrieved
    if needed."""

    def __init__(self):
        """Class constructor."""
        project_path = os.path.abspath(os.path.dirname(__file__))
        current_dir = os.getcwd()

        self.retrieved_dates = {}

        vj_genes = self._consensus_db_load(current_dir, project_path)


    def _consensus_db_load(self, current_dir, project_path):
        """Loads the current vj gene databases and returns them
        as a dictionary."""
        vj_genes = {}

        try:
            os.chdir(os.path.join(project_path, "consensus_data"))
            db_files = [f for f in os.listdir() if f.endswith(".fa.gz")]
            for db_file in db_files:
                species = db_file.split("_")[0]
                receptor = db_file.split("_")[1]

                if species not in vj_genes:
                    vj_genes[species] = {}
                    self.retrieved_dates[species] = {}

                if receptor in vj_genes[species]:
                    raise RuntimeError("Duplicate db files found for "
                            f"{species}, {receptor}.")

                vj_genes[species][receptor] = {}
                self.retrieved_dates[species][receptor] = db_file.split(".fa.gz")[0].split("_")[-1]

                # We avoid using Biopython's SeqIO (since it introduces an additional
                # unnecessary dependency). Since we wrote the db files, we know
                # how they are formatted -- one aa seq for each gene name -- and
                # can use a very simple procedure here which is not applicable to all
                # fasta files.
                with gzip.open(db_file, "rt", encoding="utf-8") as fhandle:
                    for line in fhandle:
                        if line.startswith(">"):
                            description = line.strip()[1:]
                        else:
                            aa_seq = line.strip()
                            if aa_seq in vj_genes[species][receptor]:
                                raise RuntimeError("Duplicate vj genes encountered.")
                            vj_genes[species][receptor][line.strip()] = description

        except Exception as exc:
            os.chdir(current_dir)
            raise RuntimeError("The consensus data for the package either has been deleted or "
                    "moved or was never properly installed.") from exc

        os.chdir(current_dir)
        return vj_genes



    def retrieve_db_dates(self):
        """Returns the date when the VJ gene database
        used for this assignment was retrieved."""
        return self.retrieved_dates
