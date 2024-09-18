"""Contains some simple utilities for loading VJ gene databases."""
import os
import gzip



def load_vj_gene_consensus_db(current_dir, database_path, database = "imgt"):
    """Loads the current vj gene databases and returns them
    as two dictionaries mapping species name and receptor type
    to a list of gene names and to a list of gene sequences.

    Args:
        current_dir (str): The current user directory.
        project_path (str): The location of the consensus databases."""
    vj_names, vj_seqs, retrieved_dates = {}, {}, {}

    try:
        os.chdir(database_path)
        db_files = [f for f in os.listdir() if f.endswith(".fa.gz") and
                    f.split("_")[1] == database]
        for db_file in db_files:
            species = db_file.split("_")[0]
            receptor = db_file.split("_")[2]
            if species not in retrieved_dates:
                retrieved_dates[species] = {}

            retrieved_dates[species][receptor] = db_file.split(".fa.gz")[0].split("_")[-1]

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


            key = species + "_" + receptor
            if key in vj_names:
                raise RuntimeError("Duplicate db files found for "
                            f"{species}, {receptor}.")

            vj_names[key] = names
            vj_seqs[key] = sequences

    except Exception as exc:
        os.chdir(current_dir)
        raise RuntimeError("The consensus data for the package either has been deleted or "
                    "moved or was never properly installed.") from exc

    os.chdir(current_dir)
    return vj_names, vj_seqs, retrieved_dates
