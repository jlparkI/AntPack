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


insertion_positions = {32, 60, 111}
conserved_positions = {23:["C"], 41:["W"], 104:["C"], 118:["F", "W"]}

#These are positions where either A) a blank in the query sequence is very common or
#B) insertions are very common. We want special template and query gap penalties for
#these positions.
special_positions = {10:[-26.0,-1.0], 73:[-26.0,-1.0],
        32:[-6.1,-1.0], 60:[-7.1,-1.0], 111:[-1.0,-1.0]}

# These are IMGT defined CDRs.
cdrs = {27:-11.1, 28:-10.1, 29:-9.1, 30:-8.1, 31:-7.1, 32:-6.1, 33:-6, 34:-7, 35:-8, 36:-9, 37:-10, 38:-11,
        # CDR 2
        56:-11.1, 57:-10.1, 58:-9.1, 59:-8.1, 60:-7.1, 61:-7.0, 62:-8, 63:-9, 64:-10, 65:-11,
        #CDR 3
        105:-11.1, 106:-10.1, 107:-9.1, 108:-8.1, 109:-7.1, 110:-6.1, 111:-1,
        112:-6, 113:-7, 114:-8, 115:-9, 116:-10, 117:-11}
