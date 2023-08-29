"""Defines legal inputs to hmm builder functions."""
genes = {"IG":["HV", "HJ", "KV", "KJ", "LV", "LJ"],
        "TR":["AV", "AJ", "BV", "BJ", "GV", "GJ", "DV", "DJ"]}
genes["ALL"] = genes["IG"] + genes["TR"]


urls = { "HV": "https://www.imgt.org/genedb/GENElect?query=7.3+IGHV&species=%s",
         "HJ": "https://www.imgt.org/genedb/GENElect?query=7.6+IGHJ&species=%s",
         "KV": "https://www.imgt.org/genedb/GENElect?query=7.3+IGKV&species=%s",
         "KJ": "https://www.imgt.org/genedb/GENElect?query=7.6+IGKJ&species=%s",
         "LV": "https://www.imgt.org/genedb/GENElect?query=7.3+IGLV&species=%s",
         "LJ": "https://www.imgt.org/genedb/GENElect?query=7.6+IGLJ&species=%s",
         "AV": "https://www.imgt.org/genedb/GENElect?query=7.3+TRAV&species=%s",
         "AJ": "https://www.imgt.org/genedb/GENElect?query=7.6+TRAJ&species=%s",
         "BV": "https://www.imgt.org/genedb/GENElect?query=7.3+TRBV&species=%s",
         "BJ": "https://www.imgt.org/genedb/GENElect?query=7.6+TRBJ&species=%s",
         "GV": "https://www.imgt.org/genedb/GENElect?query=7.3+TRGV&species=%s",
         "GJ": "https://www.imgt.org/genedb/GENElect?query=7.6+TRGJ&species=%s",
         "DV": "https://www.imgt.org/genedb/GENElect?query=7.3+TRDV&species=%s",
         "DJ": "https://www.imgt.org/genedb/GENElect?query=7.6+TRDJ&species=%s"
         }
allowed_species = {"IG":["Homo+sapiens",
           "Mus",
           "Rattus+norvegicus",
           "Oryctolagus+cuniculus",
           "Sus+scrofa",
           "Macaca+mulatta",
           "Vicugna+pacos",
           "Bos+taurus"],
            "TR":["Homo+sapiens",
           "Mus"]}
allowed_species["ALL"] = allowed_species["IG"]

latin_to_common = {"Homo_sapiens":"human",
           "Mus":"mouse",
           "Rattus_norvegicus":"rat",
           "Oryctolagus_cuniculus":"rabbit",
           "Macaca_mulatta":"rhesus",
           "Sus_scrofa":"pig",
           "Vicugna_pacos":"alpaca",
           "Bos_taurus":"cow"}


############################The following are alignment parameters that
#were selected to try to ensure a good alignment. Making large changes
#to these WILL change the results, possibly substantially.

DEFAULT_QUERY_GAP_PENALTY = -25
DEFAULT_TEMPLATE_GAP_PENALTY = -25

light_blank_positions = {10, 81, 82}

#These are positions that are essentially mandatory. A sequence that deviates
#from one of these is an alignment issue or has a very large deletion. Note that
#IMGT defines position 89 as highly conserved in the sense that a hydrophobic
#amino acid is always observed here, but "hydrophobic amino acid" is not defined
#and a variety of (hydrophobic) aas are observed at that position. Hence we
#handle position 89 separately below.
heavy_conserved_positions = {23:"C", 41:"W", 104:"C", 118:"W", 119:"G", 121:"G"}
light_conserved_positions = {23:"C", 41:"W", 104:"C", 118:"F", 119:"G", 121:"G"}



#These are positions where either A) a blank in the query sequence is very common or
#B) insertions are very common. We want special template and query gap penalties for
#these positions. One or two of these are chain dependent.
heavy_special_positions = {10:[-11.0,-1.0], 73:[-11.0,-1.0],
        33:[-1.0,-1.0], 61:[-1.0,-1.0], 111:[-1.0,-1.0]}

light_special_positions = {10:[-11.0,-1.0], 73:[-11.0,-1.0],
        33:[-1.0,-1.0], 61:[-1.0,-1.0], 81:[-5.0,-1.0],
        82:[-11.0,-1.0], 111:[-1.0,-1.0]}

# These are IMGT defined CDRs. We have a special gap penalty for each, which decreases towards
# the center to ensure that gaps will go there first, and in such a way that the positions
# will be filled in a particular order.
cdrs = {27:-5.1, 28:-4.6, 29:-4.1, 30:-3.6, 31:-3.1, 32:-2.6, 33:-2.5, 34:-3, 35:-3.5, 36:-4, 37:-4.5, 38:-5,
        # CDR 2
        56:-5.1, 57:-4.6, 58:-4.1, 59:-3.6, 60:-3.1, 61:-3, 62:-3.5, 63:-4, 64:-4.5, 65:-5,
        #CDR 3
        105:-5.1, 106:-4.6, 107:-4.1, 108:-3.6, 109:-3.1, 110:-2.6, 111:-1,
        112:-2.5, 113:-3, 114:-3.5, 115:-4, 116:-4.5, 117:-5}
