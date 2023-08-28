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


insertion_positions = {33, 61, 111}
conserved_positions = {23:["C"], 41:["W"], 104:["C"], 118:["F", "W"]}
# These are IMGT defined CDRs. The key is the position, the value is the gap penalty. If the position is
# also an insertion position, the gap penalty shown here is overriden by the default insertion position gap
# penalty. Penalties decrease towards the middle of the CDR to force insertions there rather than at the
# outer edge wherever possible.
cdrs = {27:-11.0, 28:-10.0, 29:-9.0, 30:-8.0, 31:-7.0, 32:-6.0, 33:-6.0, 34:-7.0, 35:-8.0, 36:-9.0, 37:-10.0, 38:-11.0,
        # CDR 2
        56:-11.0, 57:-10.0, 58:-9.0, 59:-8.0, 60:-7.0, 61:-7.0, 62:-8.0, 63:-9.0, 64:-10.0, 65:-11.0,
        #CDR 3
        105:-11.0, 106:-10.0, 107:-9.0, 108:-8.0, 109:-7.0, 110:-6.0, 111:-6.0,
        112:-6.0, 113:-7.0, 114:-8.0, 115:-9.0, 116:-10.0, 117:-11.0}
