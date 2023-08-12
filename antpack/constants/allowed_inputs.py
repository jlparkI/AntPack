"""Defines legal inputs to parser functions (e.g. what schemes are
allowed."""
#from ..numbering_funcs import number_kabat_heavy, number_kabat_light
#from ..numbering_funcs import number_chothia_heavy, number_chothia_light
#from ..numbering_funcs import number_martin_heavy, number_martin_light
#from ..numbering_funcs import number_wolfguy_heavy, number_wolfguy_light

allowed_schemes = {"martin", "chothia", "kabat", "imgt", "aho", "wolfguy"}

# Maps scheme / chain combinations to functions that implement them.
#scheme_to_fun = {("kabat", "H"):number_kabat_heavy, ("kabat", "KL"):number_kabat_light,
#        ("chothia", "H"):number_chothia_heavy, ("chothia", "KL"):number_chothia_light,
#        ("martin", "H"):number_martin_heavy, ("martin", "KL"):number_martin_light,
#        ("wolfguy", "H"):number_wolfguy_heavy, ("wolfguy", "KL"):number_wolfguy_light}
#
#Maps the chain label to the class
chain_type_to_class = {"H":"H", "K":"L", "L":"L", "A":"A", "B":"B", "G":"G", "D":"D"}

#The list of allowed amino acids.
allowed_amino_acids = {"A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N",
        "P", "Q", "R", "S", "T", "V", "W", "Y"}
