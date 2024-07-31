"""Defines useful defaults for downloading IMGT reference
sequences."""
genes = {"IG":["HV", "HJ", "KV", "KJ", "LV", "LJ"],
        "TR":["AV", "AJ", "BV", "BJ", "GV", "GJ", "DV", "DJ"]}

# For now, we download human and mouse only for IG only.
# We may expand this at a later date.
allowed_species = {"IG":["Homo_sapiens",
           "Mus_musculus"]}

latin_to_common = {"Homo_sapiens":"human",
           "Mus_musculus":"mouse",
           "Rattus_norvegicus":"rat",
           "Oryctolagus_cuniculus":"rabbit",
           "Macaca_mulatta":"rhesus",
           "Sus_scrofa":"pig",
           "Vicugna_pacos":"alpaca",
           "Bos_taurus":"cow"}
