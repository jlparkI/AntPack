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


heavy_conserved_positions = {23:"C", 41:"W", 104:"C", 118:"W", 119:"G"}
light_conserved_positions = {23:"C", 41:"W", 104:"C", 118:"F", 119:"G"}

#These are positions where either A) a blank in the query sequence is very common or
#B) insertions are very common. We want special template and query gap penalties for
#these positions.
special_positions = {10:[-25.0,-1.0], 73:[-25.0,-1.0],
        33:[-1.0,-1.0], 61:[-1.0,-1.0], 111:[-1.0,-1.0]}

# These are IMGT defined CDRs.
cdrs = {27:-5.1, 28:-4.6, 29:-4.1, 30:-3.6, 31:-3.1, 32:-2.6, 33:-2.5, 34:-3, 35:-3.5, 36:-4, 37:-4.5, 38:-5,
        # CDR 2
        56:-5.1, 57:-4.6, 58:-4.1, 59:-3.6, 60:-3.1, 61:-3, 62:-3.5, 63:-4, 64:-4.5, 65:-5,
        #CDR 3
        105:-5.1, 106:-4.6, 107:-4.1, 108:-3.6, 109:-3.1, 110:-2.6, 111:-1,
        112:-2.5, 113:-3, 114:-3.5, 115:-4, 116:-4.5, 117:-5}
