"""Defines legal inputs to parser functions (e.g. what schemes are
allowed."""
allowed_schemes = {"chothia", "kabat", "imgt"}

#Maps the chain label to the class
chain_type_to_class = {"H":"H", "K":"L", "L":"L", "A":"A", "B":"B", "G":"G", "D":"D"}

#The set of allowed amino acids.
allowed_amino_acids = {"A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N",
        "P", "Q", "R", "S", "T", "V", "W", "Y"}

#The list of allowed amino acids. Must be in alphabetical order, very important.
allowed_aa_list = sorted(list(allowed_amino_acids))
