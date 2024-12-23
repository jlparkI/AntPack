"""Tests basic functionality for the antpack cli. Note that this
relies on the functionality of the SingleChainAnnotator,
PairedChainAnnotator and VJGeneTool classes, so if those tests
do not pass, this one should not either."""
import os
import subprocess
import random
import gzip
import unittest
import numpy as np
from antpack import PairedChainAnnotator, SingleChainAnnotator, VJGeneTool



class TestAntPackCLI(unittest.TestCase):




    def test_paired_chain_csv(self):
        """Take test heavy and light chains and pair them,
        then use the CLI to do basic analysis written to a
        csv file."""
        merged_seqs, heavy_pc_annotations, light_pc_annotations, pc_aligner = get_paired_seqs()

        heavy_code_map = create_manual_mapping(heavy_pc_annotations, pc_aligner)
        light_code_map = create_manual_mapping(light_pc_annotations, pc_aligner)

        # Create a temporary fasta input file.
        with open("input_fasta_data.fasta", "w+", encoding="utf-8") as fhandle:
            for i, merged_seq in enumerate(merged_seqs):
                fhandle.write(f">this is sequence # {i}\n{merged_seq}\n")


        _ = subprocess.run(["AntPack-CLI", "input_fasta_data.fasta",
                "output_data", "imgt", "--paired"])

        self.assertTrue("output_data_heavy.csv" in os.listdir())
        self.assertTrue("output_data_light.csv" in os.listdir())

        with open("output_data_heavy.csv", "r",
                encoding="utf-8") as fhandle:
            _ = fhandle.readline()
            for i, line in enumerate(fhandle):
                elements = line.strip().split(",")
                self.assertTrue(elements[0] == f"this is sequence # {i}")
                self.assertTrue("".join(elements[2:-1])==map_sequence(merged_seqs[i],
                    heavy_pc_annotations[i], heavy_code_map))
                self.assertTrue(elements[-1]==heavy_pc_annotations[i][-1])

            self.assertTrue((i+1)==len(heavy_pc_annotations))

        with open("output_data_light.csv", "r",
                encoding="utf-8") as fhandle:
            temp = fhandle.readline()
            for i, line in enumerate(fhandle):
                elements = line.strip().split(",")
                self.assertTrue(elements[0] == f"this is sequence # {i}")
                self.assertTrue("".join(elements[2:-1])==map_sequence(merged_seqs[i],
                    light_pc_annotations[i], light_code_map))
                self.assertTrue(elements[-1]==light_pc_annotations[i][-1])

            self.assertTrue((i+1)==len(light_pc_annotations))

        os.remove("input_fasta_data.fasta")
        os.remove("output_data_heavy.csv")
        os.remove("output_data_light.csv")




    def test_paired_chain_vj(self):
        """Take test heavy and light chains and pair them,
        then use the CLI to do basic analysis written to a
        csv file with human VJ gene assignment."""
        merged_seqs, heavy_pc_annotations, light_pc_annotations, pc_aligner = get_paired_seqs()
        vj_tool = VJGeneTool(database="imgt", scheme="imgt")

        heavy_code_map = create_manual_mapping(heavy_pc_annotations, pc_aligner)
        light_code_map = create_manual_mapping(light_pc_annotations, pc_aligner)

        # Create a temporary fasta input file.
        with open("input_fasta_data.fasta", "w+", encoding="utf-8") as fhandle:
            for i, merged_seq in enumerate(merged_seqs):
                fhandle.write(f">this is sequence # {i}\n{merged_seq}\n")


        _ = subprocess.run(["AntPack-CLI", "input_fasta_data.fasta", "output_data",
            "imgt", "--paired", "--vj", "human", "identity"])

        self.assertTrue("output_data_heavy.csv" in os.listdir())
        self.assertTrue("output_data_light.csv" in os.listdir())

        with open("output_data_heavy.csv", "r",
                encoding="utf-8") as fhandle:
            _ = fhandle.readline()

            for i, line in enumerate(fhandle):
                elements = line.strip().split(",")
                self.assertTrue(elements[0] == f"this is sequence # {i}")
                self.assertTrue(elements[2] == "human")
                self.assertTrue(elements[3] == "identity")
                v_gene, j_gene, vident, jident = vj_tool.assign_vj_genes(heavy_pc_annotations[i],
                        merged_seqs[i], "human", "identity")
                self.assertTrue(elements[4]==v_gene)
                self.assertTrue(np.allclose(float(elements[5]), vident))
                self.assertTrue(elements[6]==j_gene)
                self.assertTrue(np.allclose(float(elements[7]), jident))
                self.assertTrue("".join(elements[8:-1])==map_sequence(merged_seqs[i],
                    heavy_pc_annotations[i], heavy_code_map))
                self.assertTrue(elements[-1]==heavy_pc_annotations[i][-1])

            self.assertTrue((i+1)==len(heavy_pc_annotations))

        with open("output_data_light.csv", "r",
                encoding="utf-8") as fhandle:
            _ = fhandle.readline()
            for i, line in enumerate(fhandle):
                elements = line.strip().split(",")
                self.assertTrue(elements[0] == f"this is sequence # {i}")
                self.assertTrue(elements[2] == "human")
                self.assertTrue(elements[3] == "identity")
                v_gene, j_gene, vident, jident = vj_tool.assign_vj_genes(light_pc_annotations[i],
                        merged_seqs[i], "human", "identity")
                self.assertTrue(elements[4]==v_gene)
                self.assertTrue(np.allclose(float(elements[5]), vident))
                self.assertTrue(elements[6]==j_gene)
                self.assertTrue(np.allclose(float(elements[7]), jident))
                self.assertTrue("".join(elements[8:-1])==map_sequence(merged_seqs[i],
                    light_pc_annotations[i], light_code_map))
                self.assertTrue(elements[-1]==light_pc_annotations[i][-1])

            self.assertTrue((i+1)==len(light_pc_annotations))

        os.remove("input_fasta_data.fasta")
        os.remove("output_data_heavy.csv")
        os.remove("output_data_light.csv")






    def test_single_chain_csv(self):
        """Take test heavy and light chains,
        then use the CLI to do basic analysis and write
        to fasta, using single chain containing sequences
        as input only."""
        seqs, annotations, sc_aligner = get_unpaired_seqs()
        heavy_code_map = create_manual_mapping([a for a in annotations if
            a[2] == 'H'], sc_aligner)
        light_code_map = create_manual_mapping([a for a in annotations if
            a[2] != 'H'], sc_aligner)

        # Create a temporary fasta input file.
        with open("input_fasta_data.fasta", "w+", encoding="utf-8") as fhandle:
            for i, seq in enumerate(seqs):
                fhandle.write(f">this is sequence # {i}\n{seq}\n")

        _ = subprocess.run(["AntPack-CLI", "input_fasta_data.fasta", "output_data",
            "imgt"])

        self.assertTrue("output_data_heavy.csv" in os.listdir())
        self.assertTrue("output_data_light.csv" in os.listdir())

        ntot = 0
        for chain_type, code_map in zip(["heavy", "light"], [heavy_code_map,
            light_code_map]):
            with open(f"output_data_{chain_type}.csv", "r",
                    encoding="utf-8") as fhandle:
                _ = fhandle.readline()
                for i, line in enumerate(fhandle):
                    ntot += 1
                    elements = line.strip().split(",")
                    code = int(elements[0].split("this is sequence # ")[1])
                    self.assertTrue(annotations[code][-1] == elements[-1])
                    self.assertTrue("".join(elements[2:-1])==map_sequence(seqs[code],
                        annotations[code], code_map))

        self.assertTrue(ntot==len(annotations))

        os.remove("input_fasta_data.fasta")
        os.remove("output_data_heavy.csv")
        os.remove("output_data_light.csv")



    def test_single_chain_vj(self):
        """Take test heavy and light chains,
        then use the CLI to do basic analysis plus vj assignment
        and write to fasta, using single chain containing sequences
        as input only."""
        seqs, annotations, sc_aligner = get_unpaired_seqs()
        vj_tool = VJGeneTool(scheme="imgt", database="imgt")

        heavy_code_map = create_manual_mapping([a for a in annotations if
            a[2] == 'H'], sc_aligner)
        light_code_map = create_manual_mapping([a for a in annotations if
            a[2] != 'H'], sc_aligner)

        # Create a temporary fasta input file.
        with open("input_fasta_data.fasta", "w+", encoding="utf-8") as fhandle:
            for i, seq in enumerate(seqs):
                fhandle.write(f">this is sequence # {i}\n{seq}\n")

        _ = subprocess.run(["AntPack-CLI", "input_fasta_data.fasta", "output_data",
            "imgt", "--vj", "human", "identity"])

        self.assertTrue("output_data_heavy.csv" in os.listdir())
        self.assertTrue("output_data_light.csv" in os.listdir())

        ntot = 0
        for chain_type, code_map in zip(["heavy", "light"], [heavy_code_map,
            light_code_map]):
            with open(f"output_data_{chain_type}.csv", "r",
                    encoding="utf-8") as fhandle:
                _ = fhandle.readline()
                for i, line in enumerate(fhandle):
                    ntot += 1
                    elements = line.strip().split(",")
                    code = int(elements[0].split("this is sequence # ")[1])
                    self.assertTrue(elements[2] == "human")
                    self.assertTrue(elements[3] == "identity")
                    v_gene, j_gene, vident, jident = vj_tool.assign_vj_genes(annotations[code],
                            seqs[code], "human", "identity")
                    self.assertTrue(elements[4]==v_gene)
                    self.assertTrue(np.allclose(float(elements[5]), vident))
                    self.assertTrue(elements[6]==j_gene)
                    self.assertTrue(np.allclose(float(elements[7]), jident))
                    self.assertTrue(annotations[code][-1] == elements[-1])
                    self.assertTrue("".join(elements[8:-1])==map_sequence(seqs[code],
                        annotations[code], code_map))

        self.assertTrue(ntot==len(annotations))

        os.remove("input_fasta_data.fasta")
        os.remove("output_data_heavy.csv")
        os.remove("output_data_light.csv")



def create_manual_mapping(annotation_list, pc_tool,
        scheme = "imgt", add_unobserved_positions = True):
    """Manually creates a dictionary mapping codes to positions."""
    code_set = set()

    for annotation in annotation_list:
        code_set.update(annotation[0])

    if add_unobserved_positions:
        if scheme == "imgt":
            for i in range(1, 129):
                code_set.add(str(i))

    code_set = list(code_set)
    sorted_codes = pc_tool.sort_position_codes(code_set)
    return {k:i for i, k in enumerate(sorted_codes)}



def map_sequence(sequence, annotation, code_dict):
    """Maps an input sequence to a standard length vector to
    create an MSA."""
    output_seq = ["-" for k in code_dict]
    for letter, code in zip(sequence, annotation[0]):
        if code == "-":
            continue
        output_seq[code_dict[code]] = letter
    return "".join(output_seq)


def get_paired_seqs():
    """Get a set of paired chain sequences for testing."""
    project_path = os.path.abspath(os.path.dirname(__file__))
    current_dir = os.getcwd()
    os.chdir(os.path.join(project_path, "test_data"))
    with gzip.open("test_data.csv.gz", "rt") as fhandle:
        _ = fhandle.readline()
        seqs = [line.strip().split(",")[0] for line in fhandle]
    os.chdir(current_dir)

    sc_aligner = SingleChainAnnotator(scheme="imgt")
    pc_aligner = PairedChainAnnotator(scheme="imgt")

    alignments = [sc_aligner.analyze_seq(seq) for seq in seqs]

    heavy_chains = [(seq, a) for a, seq in zip(alignments, seqs)
                if a[2] == "H"]
    light_chains = [(seq, a) for a, seq in zip(alignments, seqs)
                if a[2] != "H"]

    random.seed(0)

    # Join the chains using a very arbitrary linker.
    merged_seqs = ["".join( [light_chain[0], "YYYYGGGGSSSS", heavy_chain[0]] )
                for (heavy_chain, light_chain) in zip(heavy_chains, light_chains)]

    heavy_pc_annotations, light_pc_annotations = [], []

    for merged_seq in merged_seqs:
        h_annot, l_annot = pc_aligner.analyze_seq(merged_seq)
        heavy_pc_annotations.append(h_annot)
        light_pc_annotations.append(l_annot)

    return merged_seqs, heavy_pc_annotations, light_pc_annotations, pc_aligner



def get_unpaired_seqs():
    """Get a set of paired chain sequences for testing."""
    project_path = os.path.abspath(os.path.dirname(__file__))
    current_dir = os.getcwd()
    os.chdir(os.path.join(project_path, "test_data"))
    with gzip.open("test_data.csv.gz", "rt") as fhandle:
        _ = fhandle.readline()
        seqs = [line.strip().split(",")[0] for line in fhandle]
    os.chdir(current_dir)

    sc_aligner = SingleChainAnnotator(scheme="imgt")

    annotations = [sc_aligner.analyze_seq(seq) for seq in seqs]

    return seqs, annotations, sc_aligner


if __name__ == "__main__":
    unittest.main()
