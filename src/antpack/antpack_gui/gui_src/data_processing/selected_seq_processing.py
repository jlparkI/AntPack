"""Contains tools for processing a selected seq (i.e. numbering
it, finding adjacent VJ genes and prepping for presentation)."""
from antpack import VJGeneTool



class VJComparisonData:


    def __init__(self, heavy_nmbr, light_nmbr, seq,
            sc_tool, vj_tool, species):
        self.full_sequence = seq

        self.heavy_seqs = []
        self.heavy_pid = []
        self.heavy_err = []
        self.heavy_nmbr = []
        self.heavy_descrip = []

        self.light_seqs = []
        self.light_pid = []
        self.light_err = []
        self.light_nmbr = []
        self.light_descrip = []

        # Note that v and j genes in AntPack are
        # stored pre-numbered using IMGT. TODO:
        # Need to add a converter to switch
        # IMGT-numbered genes to non-IMGT for
        # easier alignment and display.
        if heavy_nmbr is not None:
            self.heavy_seqs, self.heavy_nmbr, self.heavy_pid, self.heavy_err, self.heavy_descrip = \
                    self.process_input_numbering(heavy_nmbr, seq,
                    sc_tool, vj_tool, species)

        if light_nmbr is not None:
            self.light_seqs, self.light_nmbr, self.light_pid, self.light_err, self.light_descrip = \
                    self.process_input_numbering(light_nmbr, seq,
                    sc_tool, vj_tool, species)



    def get_num_heavy(self):
        """Returns number of loaded heavy sequences."""
        return len(self.heavy_seqs)


    def get_num_light(self):
        """Returns number of loaded light sequences."""
        return len(self.light_seqs)


    def get_heavy_data(self):
        """Returns the lists associated with the heavy chain."""
        return self.heavy_seqs, self.heavy_pid, self.heavy_err, self.heavy_descrip

    def get_light_data(self):
        """Returns the lists associated with the light chain."""
        return self.light_seqs, self.light_pid, self.light_err, self.light_descrip

    def get_heavy_numbering(self):
        """Returns the heavy numbering (for the spreadsheet header)."""
        return self.heavy_nmbr

    def get_light_numbering(self):
        """Returns the light numbering (for the spreadsheet header)."""
        return self.light_nmbr




    def process_input_numbering(self, alignment, sequence,
            sc_tool, vj_tool, species):
        """Converts an input sequence and alignment to lists
        of information that can be readily displayed."""
        seqs, numbering, pid, err = [], [], [], []

        seqs.append(sequence)
        numbering.append(alignment)
        pid.append(str(alignment[1]))
        err.append(alignment[3])
        descrip = ["Input sequence"]

        vj = vj_tool.assign_vj_genes(alignment, sequence, species)
        vjgenes = vj[0].split("_") + vj[1].split("_")

        if len(vjgenes) > 0:
            for vgene in vjgenes:
                vseq = vj_tool.get_vj_gene_sequence(vgene, species)
                if vseq != "":
                    seqs.append(vseq)
                    numbering.append( ([str(z) for z in range(1,129)], 0.,
                        alignment[2], "") )
                    pid.append("")
                    err.append("")
                    descrip.append(vgene)

        numbering, seqs = sc_tool.build_msa(seqs, numbering)
        return seqs, numbering, pid, err, descrip



class MultiSequenceData:


    def __init__(self):
        self.full_sequence = seq

        self.heavy_seqs = []
        self.heavy_pid = []
        self.heavy_err = []
        self.heavy_nmbr = []
        self.seq_names = []

        self.light_seqs = []
        self.light_pid = []
        self.light_err = []
        self.light_nmbr = []



    def get_num_heavy(self):
        """Returns number of loaded heavy sequences."""
        return len(self.heavy_seqs)


    def get_num_light(self):
        """Returns number of loaded light sequences."""
        return len(self.light_seqs)


    def get_heavy_data(self):
        """Returns the lists associated with the heavy chain."""
        return self.heavy_seqs, self.heavy_pid, self.heavy_err, self.seq_names

    def get_light_data(self):
        """Returns the lists associated with the light chain."""
        return self.light_seqs, self.light_pid, self.light_err, self.seq_names

    def get_heavy_numbering(self):
        """Returns the heavy numbering (for the spreadsheet header)."""
        return self.heavy_nmbr

    def get_light_numbering(self):
        """Returns the light numbering (for the spreadsheet header)."""
        return self.light_nmbr


    def add_sequence(self, seq, heavy_numbering, light_numbering):
        """Adds a sequence which may contain a heavy chain, a light
        chain or both to the sequence list. Note that unlike the
        VJ comparison, if a sequence contains no heavy or light
        chain, we add a blank sequence (so that heavy and light
        are always the same length)."""





def process_for_vj_comparison(seq, seq_type,
        sc_annotator, pc_annotator,
        vj_tool, liability_tool, scheme,
        pid_threshold = 0.7,
        species = "human"):
    """Get numbering and VJ genes for the input
    sequence. If necessary determine its type.
    Calls to this function should always be
    wrapped in try-except since invalid entries
    (e.g. non amino acid letters) may trigger
    an exception."""
    heavy_nmbr, light_nmbr = None, None
    

    if seq_type == "single":
        nmbr = sc_annotator.analyze_seq(seq)
        if nmbr[1] >= pid_threshold:
            if nmbr[2] == "H":
                heavy_nmbr = nmbr
            else:
                light_nmbr = nmbr
    else:
        heavy_nmbr, light_nmbr = pc_annotator.analyze_seq(seq)
        if heavy_nmbr[1] < pid_threshold:
            heavy_nmbr = None
        if light_nmbr[1] < pid_threshold:
            light_nmbr = None

    if heavy_nmbr is None and light_nmbr is None:
        return None

    ssd = VJComparisonData(heavy_nmbr, light_nmbr,
            seq, sc_annotator, vj_tool, species)
    return ssd
