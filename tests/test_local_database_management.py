"""Test database management & search."""
import os
import random
import unittest
from antpack import (build_database_from_fasta,
        SingleChainAnnotator, LocalDBTool, VJGeneTool)
from antpack.antpack_cpp_ext import (SequenceTemplateAligner,
        return_imgt_canonical_numbering_cpp)
from antpack.utilities import read_fasta


AAMAP = {k:i for i,k in enumerate("ACDEFGHIKLMNPQRSTVWY-")}


class TestLocalDBManagement(unittest.TestCase):


    def test_local_db_search(self):
        """Check that searches retrieve the sequences that
        we would expect based on a brute-force exact search."""
        current_dir = os.path.abspath(os.path.dirname(
            __file__))
        data_filepath = os.path.join(current_dir,
                "test_data", "addtnl_test_data.fasta.gz")

        if not os.path.exists('TEMP_DB.db'):
            build_database_from_fasta([data_filepath],
                "TEMP_DB.db", numbering_scheme="imgt",
                cdr_definition_scheme="imgt",
                sequence_type="single", receptor_type="mab",
                pid_threshold=0.7, user_memo="")

        seqs, seqinfos = [], []

        for seqinfo, seq in read_fasta(data_filepath):
            seqs.append(seq)
            seqinfos.append(seqinfo)

        local_db = LocalDBTool("TEMP_DB.db")

        nmbr_scheme = "imgt"
        cdr_scheme = "imgt"

        # First annotate the heavy chains and build an MSA.
        # Use this msa to do brute force searches and compare
        # the results with those we get back from the database.
        sca = SingleChainAnnotator(scheme=nmbr_scheme)
        annotations = sca.analyze_seqs(seqs)

        chains = [(a,s) for (a,s) in zip(annotations,
            seqs) if a[2]=="H"]
        codes, msa = sca.build_msa([t[1] for t in chains],
                [t[0] for t in chains])

        self.random_search_function_test(msa, codes,
                sca, local_db, chains,
                nmbr_scheme, cdr_scheme)

        del local_db
        os.remove("TEMP_DB.db")



    def test_local_db_search_setup(self):
        """Check that the local db management tool is set
        up correctly."""
        current_dir = os.path.abspath(os.path.dirname(
            __file__))
        data_filepath = os.path.join(current_dir,
                "test_data", "addtnl_test_data.fasta.gz")

        if not os.path.exists('TEMP_DB.db'):
            build_database_from_fasta([data_filepath],
                "TEMP_DB.db", numbering_scheme="imgt",
                cdr_definition_scheme="imgt",
                sequence_type="single", receptor_type="mab",
                pid_threshold=0.7, user_memo="")

        seqs, seqinfos = [], []

        for seqinfo, seq in read_fasta(data_filepath):
            seqs.append(seq)
            seqinfos.append(seqinfo)

        local_db = LocalDBTool("TEMP_DB.db")

        # Check that the search tool retrieves the correct metadata.
        metadata = local_db.get_database_metadata()
        self.assertTrue('imgt' == metadata[0])
        self.assertTrue('mab' == metadata[1])
        self.assertTrue('single' == metadata[2])
        self.assertTrue('imgt' == metadata[3])
        self.assertTrue('' == metadata[4])

        # Next check that the position counts table matches
        # expected. First check for heavy. The database
        # returns the tables all together but the Python
        # test code will check them separately.
        mtables, dtables = local_db.get_database_counts()

        sca = SingleChainAnnotator(scheme="imgt")
        annotations = sca.analyze_seqs(seqs)

        numbering = return_imgt_canonical_numbering_cpp()
        sta = SequenceTemplateAligner(numbering, "H", "imgt",
                "imgt")

        heavy_chains = [(a,s) for (a,s) in zip(annotations,
            seqs) if a[2]=="H"]
        msa = [sta.align_sequence(h[1], h[0][0], False) for
                h in heavy_chains]

        labels = sca.assign_cdr_labels(numbering,
                "H", "imgt")
        gt_dtable, gt_ttable = get_kmer_counts(msa, labels)

        for i in range(3):
            self.assertTrue(mtables[0][i]==gt_dtable[i])
            self.assertTrue(dtables[0][i]==gt_ttable[i])

        del numbering, msa, labels, sta



        # ...then light.
        numbering = return_imgt_canonical_numbering_cpp()
        sta = SequenceTemplateAligner(numbering, "L", "imgt",
                "imgt")

        light_chains = [(a,s) for (a,s) in zip(annotations,
            seqs) if a[2]!="H"]
        msa = [sta.align_sequence(h[1], h[0][0], False) for
                h in light_chains]

        labels = sca.assign_cdr_labels(numbering,
                "L", "imgt")
        gt_dtable, gt_ttable = get_kmer_counts(msa, labels)

        for i in range(3):
            self.assertTrue(mtables[1][i]==gt_dtable[i])
            self.assertTrue(dtables[1][i]==gt_ttable[i])

        del numbering, msa, labels, sta

        del local_db
        os.remove("TEMP_DB.db")




    def random_search_function_test(self, msa, codes,
        annotator, db_tool, chains,
        nmbr_scheme = "imgt", cdr_scheme = "imgt"):
        """Runs a test for a supplied set of chains and
        msas using randomly selected search parameters."""
        chain_code = chains[0][0][2]

        labels = annotator.assign_cdr_labels(codes, chain_code, cdr_scheme)
        vj = VJGeneTool(scheme=nmbr_scheme)
        vgenes = [vj.assign_vj_genes(c[0], c[1], "unknown",
            "identity")[0] for c in chains]

        AALIST = list("ACDEFGHIKLMNPQRSTVWY")

        # Randomly select a sequence, randomly mutated it and
        # find its closest matches using a random setting for
        # distance and other parameters.
        random.seed(123)

        for ctr in range(250):
            idx = random.randint(0, len(msa) - 1)
            query_seq = list(msa[idx])
            mut_percentage = random.uniform(0, 0.2)
            max_hits = random.randint(1,100)

            for i, qletter in enumerate(query_seq):
                if qletter == '-':
                    continue
                if labels[i].startswith("cdr") and \
                        random.uniform(0,1) < mut_percentage:
                    query_seq[i] = random.choice(AALIST)

            query_seq = "".join(query_seq)

            if random.randint(0,1) == 1:
                cutoff = 0.33
            else:
                cutoff = 0.25

            if random.randint(0,1) == 1:
                mode = "123"
            else:
                mode = "3"
            if random.randint(0,1) == 1:
                max_cdr_length_shift = 2
            else:
                max_cdr_length_shift = 0
            if random.randint(0,1) == 1:
                retrieve_closest_only = True
            else:
                retrieve_closest_only = False
            if random.randint(0,4) == 4:
                vgene, _, _, _, species = vj.assign_vj_genes(
                        chains[idx][0], chains[idx][1],
                        "unknown", "identity")
            else:
                species, vgene = "", ""

            hit_idx, hit_dists, _ = db_tool.search(
                    query_seq, (codes, 1, chain_code, ""),
                    mode, cutoff, max_cdr_length_shift,
                    max_hits, retrieve_closest_only,
                    vgene, species)

            gt_hit_idx, gt_hit_dists = perform_exact_search(query_seq,
                    msa, labels, mode, cutoff, max_cdr_length_shift,
                    retrieve_closest_only, max_hits,
                    vgenes, vgene)

            # Necessary because the cpp search routine adds 1 to each
            # hit idx.
            hit_idx = [h-1 for h in hit_idx]

            print(f"{ctr}, {mut_percentage} mut percentage, cutoff {cutoff}, "
                    f"max hits {max_hits}, vgene_filter {vgene}, "
                    f"num gt hits was {len(gt_hit_idx)}", flush=True)
            for hdist in list(set(gt_hit_dists)):
                h1 = [h for h,d in zip(gt_hit_idx, gt_hit_dists)
                        if d == hdist]
                h1.sort()
                h2 = [h for h,d in zip(hit_idx, hit_dists)
                        if d == hdist]
                h2.sort()
                if h1 != h2:
                    import pdb
                    pdb.set_trace()
                self.assertTrue(h1==h2)


def perform_exact_search(query, msa, labels, mode, cdr_cutoff,
        max_cdr_length_shift, retrieve_closest_only,
        max_hits, vgenes, vgene_filter=""):
    """Performs an exact search on an input msa, finding
    all sequences that have hamming distance <= specified."""
    retained_dists, hit_idx, cdrlen = [], [], []
    max_hamming = 0

    for j, region in enumerate(["cdr1", "cdr2", "cdr3"]):
        cdrlen.append(len([q for (q,l) in zip(query, labels)
            if q != '-' and l == region]))
        if str(j+1) not in mode:
            continue
        max_hamming += round(cdrlen[-1] * cdr_cutoff)

    max_hamming_cdr3 = cdrlen[-1] * cdr_cutoff

    vgene_cutoff = -1
    for i in range(4, len(vgene_filter)):
        if not vgene_filter[i].isnumeric():
            vgene_cutoff = i
            break

    for i, seq in enumerate(msa):
        net_dist = 0
        if vgene_filter != "":
            if vgenes[i][:vgene_cutoff] != vgene_filter[:vgene_cutoff]:
                continue

        for j, region in enumerate(["cdr1", "cdr2", "cdr3"]):
            if str(j+1) not in mode:
                continue
            region_dist = len([a for (a,b,l) in zip(seq, query,
                labels) if l==region and a != b])
            net_dist += region_dist
            region_len = len([s for (s,l) in zip(seq, labels)
                if l==region and s != '-'])
            # Add an arbitrary large number if the cdr length is
            # unacceptable or if an individual distance is
            # unacceptable so that the sequence is not saved.
            if region_len > cdrlen[j] + max_cdr_length_shift or \
                    region_len < cdrlen[j] - max_cdr_length_shift:
                net_dist += 100
            # Do the same if this is cdr3 and the distance for cdr3
            # only is out of bounds.
            if region == "cdr3" and region_dist > max_hamming_cdr3:
                net_dist += 100

        if net_dist <= max_hamming:
            hit_idx.append(i)
            retained_dists.append(net_dist)

    if len(hit_idx) == 0:
        return [], []

    if retrieve_closest_only:
        best_dist = min(retained_dists)
        output_idx = [h for (h,c) in zip(hit_idx,
            retained_dists) if c == best_dist]
        output_dists = [best_dist for o in output_idx]
        return output_idx, output_dists

    hit_idx = [x for _, x in sorted(zip(retained_dists,
        hit_idx), key=lambda pair: pair[0])]
    retained_dists.sort()

    return hit_idx[:max_hits], retained_dists[:max_hits]


def get_kmer_counts(msa, cdr_labels):
    """Extract kmers from the cdrs of an msa and
    count the number present of each."""
    dimer_count_tables, trimer_count_tables = [], []

    for t, (cdr, next_fmwk) in enumerate([("cdr1", "fmwk2"),
            ("cdr2", "fmwk3"), ("cdr3", "fmwk4")]):
        cdr_start = cdr_labels.index(cdr)
        cdr_end = cdr_labels.index(next_fmwk)
        full_cdrlen = cdr_end - cdr_start

        dimer_count_tables.append(
                [0]*int((full_cdrlen + 1)/2)*21*21)
        trimer_count_tables.append(
                [0]*int((full_cdrlen + 2)/3)*21*21*21)

        # TODO: For now we are indexing on cdr3 only. Update this.
        if cdr != "cdr3":
            continue

        ctr = 0
        for i in range(cdr_start, cdr_end, 2):
            if i < cdr_end - 1:
                for m in msa:
                    kmer_code = AAMAP[m[i]] * 21 + AAMAP[m[i+1]]
                    dimer_count_tables[t][kmer_code+ctr] += 1
            else:
                for m in msa:
                    kmer_code = AAMAP[m[i]] * 21 + 20
                    dimer_count_tables[t][kmer_code+ctr] += 1
            ctr += 21*21

        ctr = 0
        for i in range(cdr_start, cdr_end, 3):
            for m in msa:
                cutpoint = min(cdr_end, i+3)
                trimer = m[i:cutpoint]
                while len(trimer) < 3:
                    trimer += '-'
                kmer_code = AAMAP[trimer[0]] * 21 * 21 + \
                            AAMAP[trimer[1]] * 21 + \
                            AAMAP[trimer[2]]
                trimer_count_tables[t][kmer_code+ctr] += 1
            ctr += 21*21*21

    return dimer_count_tables, trimer_count_tables



if __name__ == "__main__":
    unittest.main()
