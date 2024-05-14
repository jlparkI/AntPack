"""This module contains tools needed to assign the closest VJ germline
VJ gene (highest percent identity), using amino acid information
only. Bear in mind that nucleotide assignments are likely more
reliable, and that tools which offer a probabilistic assigment
(e.g. IGOR) may be more informative."""
import os
import gzip
from ..numbering_tools.single_chain_annotator import SingleChainAnnotator
from antpack_cpp_ext import validate_sequence, VJMatchCounter



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
        self.default_aligner = SingleChainAnnotator()

        self.vj_gene_matchups = self._consensus_db_load(current_dir, project_path)
        self.std_positions = {str(i+1) for i in range(128)}



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

                self.retrieved_dates[species][receptor] = db_file.split(".fa.gz")[0].split("_")[-1]

                # We avoid using Biopython's SeqIO (since it introduces an additional
                # unnecessary dependency). Since we wrote the db files, we know
                # how they are formatted -- one aa seq for each gene name -- and
                # can use a very simple procedure here which is not applicable to all
                # fasta files.
                with gzip.open(db_file, "rt", encoding="utf-8") as fhandle:
                    sequences, names = [], []

                    for line in fhandle:
                        if line.startswith(">"):
                            description = line.strip()[1:]
                        else:
                            aa_seq = line.strip()
                            sequences.append(aa_seq)
                            names.append(description)

                vj_genes[species][receptor] = VJMatchCounter(sequences, names)

        except Exception as exc:
            os.chdir(current_dir)
            raise RuntimeError("The consensus data for the package either has been deleted or "
                    "moved or was never properly installed.") from exc

        os.chdir(current_dir)
        return vj_genes


    def get_vj_gene_sequence(self, vj_gene_name, species = "human"):
        """Retrieves the amino acid sequence of a specified V or J
        gene, if it is in AntPack's current database. You can use
        assign_sequence or assign_numbered_sequence to get VJ gene
        names and this function to see what the VJ sequences are (if
        needed).

        Args:
            vj_gene_name (str): A valid V or J gene name, as generated
                by for example assign_sequence.
            species (str): One of 'human', 'mouse'.

        Returns:
            sequence (str): The amino acid sequence of the V or J gene
                that was requested. If that V or J gene name does not
                match anything in AntPack's internal database, None
                is returned.
        """
        if species not in self.vj_gene_matchups:
            return None
        if vj_gene_name[:4] not in self.vj_gene_matchups[species]:
            return None

        match_tool = self.vj_gene_matchups[species][vj_gene_name[:4]]
        match_result = match_tool.findVJSequenceByName(vj_gene_name)
        if match_result == "":
            return None
        return match_result



    def get_vj_gene_family(self, vj_family, species = "human"):
        """Retrieves the amino acid sequences and names of all V or
        J genes in a specified family (e.g. IGHV1) that are
        currently in AntPack's database.

        Args:
            vj_family (str): A valid V or J gene family (e.g. IGHV1).
            species (str): One of 'human', 'mouse'.

        Returns:
            sequences (list): A list of the v or j gene sequences in
                the specified family. The list is empty if the family
                supplied was not found.
            names (list): A list of the v or j gene names in this family.
                The list is empty if the family was not found.
        """
        if species not in self.vj_gene_matchups:
            return [], []
        if vj_family[:4] not in self.vj_gene_matchups[species]:
            return [], []

        seq_list, name_list = self.vj_gene_matchups[species][vj_family[:4]].getSeqLists()
        idx = [i for i, name in enumerate(name_list) if name.startswith(vj_family)]
        if len(idx) == 0:
            return [], []
        sequences, names = [seq_list[i] for i in idx], [name_list[i] for i in idx]
        return sequences, names




    def assign_numbered_sequence(self, sequence, numbering,
            chain, species = "human"):
        """Assigns V and J genes for a sequence which has already been
        numbered, either by AntPack or some other tool.

        Args:
            sequence (str): A sequence containing the usual 20 amino acids -- no gaps.
            numbering (list): The numbering has to be in the format returned by the
                SingleChainAnnotator tool, which is a list where each element is a
                character that is either "-" or an allowed IMGT number.
            chain (str): Must be one of "H", "K" or "L". This information is returned
                by SingleChainAnnotator if you use SingleChainAnnotator to number the
                sequence.
            species (str): Currently must be one of 'human', 'mouse'. More options will
                be added in future.

        Returns:
            v_gene (str): The closest V-gene name, as measured by sequence identity.
            j_gene (str): The closest J-gene name, as measured by sequence identity.

        Raises:
            RuntimeError: A RuntimeError is raised if a sequence containing invalid amino
                acids is passed.
        """
        if chain not in ["H", "K", "L"]:
            raise RuntimeError("Unexpected chain passed to vj gene matcher.")

        fmt_seq = self._prep_sequence(sequence, numbering)

        v_matcher = self.vj_gene_matchups[species][f"IG{chain}V"]
        j_matcher = self.vj_gene_matchups[species][f"IG{chain}J"]

        v_gene, _ = v_matcher.vjMatch(fmt_seq)
        j_gene, _ = j_matcher.vjMatch(fmt_seq)
        return v_gene, j_gene


    def assign_sequence(self, sequence, species = "human"):
        """Aligns a sequence using SingleChainAnnotator then assigns
        it to V and J gene names. Use this function if your sequence
        is not yet numbered and if it is a single chain (to extract
        multiple chains from a multi-chain sequence, you can use
        MultiChainAnnotator).

        Args:
            sequence (str): A sequence containing the usual 20 amino acids, no gaps.
            species (str): One of 'human', 'mouse'.

        Returns:
            v_gene (str): The closest V-gene name, as measured by sequence identity.
                If there is an error in alignment None is returned.
            j_gene (str): The closest J-gene name, as measured by sequence identity.
                If there is an error in alignment None is returned.

        Raises:
            RuntimeError: A RuntimeError is raised if a sequence containing invalid amino
                acids is passed or if the alignment suggests the sequence may be
                problematic.
        """
        numbering, p_ident, chain, err = self.default_aligner.analyze_seq(sequence)
        if p_ident < 0.8 or err != "":
            return None, None

        fmt_seq = self._prep_sequence(sequence, numbering)

        v_matcher = self.vj_gene_matchups[species][f"IG{chain}V"]
        j_matcher = self.vj_gene_matchups[species][f"IG{chain}J"]

        v_gene, _ = v_matcher.vjMatch(fmt_seq)
        j_gene, _ = j_matcher.vjMatch(fmt_seq)
        return v_gene, j_gene


    def _prep_sequence(self, sequence, numbering):
        """Extracts the standard IMGT 128 positions from a numbered
        sequence and represents others as blanks.

        Args:
            sequence (str): A sequence containing the standard 20 AAs.
            numbering (list): An Antpack-style formatted numbering
                assignment.

        Returns:
            fmt_seq (str): A formatted length 128 sequence ready for
                input to the Cpp extension.
        """
        if not validate_sequence(sequence):
            raise RuntimeError("Sequence contains invalid aas.")

        if len(numbering) != len(sequence):
            raise RuntimeError("Sequence length and numbering length must "
                    "match.")

        fmt_seq = ['-' for i in range(128)]
        for letter, nb_assign in zip(sequence, numbering):
            if nb_assign in self.std_positions:
                fmt_seq[int(nb_assign) - 1] = letter

        return "".join(fmt_seq)



    def retrieve_db_dates(self):
        """Returns the dates when each VJ gene database
        used for this assignment was last updated by downloading
        from IMGT."""
        return self.retrieved_dates
