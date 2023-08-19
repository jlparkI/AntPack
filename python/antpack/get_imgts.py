"""
Contains the tools needed to download selected IMGT files and convert them
into an HMM profile that can be used for scoring or numbering.
The parent sequences are found here:

https://www.imgt.org/vquest/refseqh.html     
"""

from html.parser import HTMLParser
from html.entities import name2codepoint
import urllib.request, urllib.parse, urllib.error, os, sys
from .constants import hmmbuild_constants as hmbc




class IMGTBuilder:
    """Class that provides the tools needed to construct HMM profiles
    from IMGT database germline sequences. These are all provided as
    static methods for convenience.
    """


    @staticmethod
    def validate_inputs(species, receptor_type):
        """Checks user specified options for validity, and returns a list
        of urls to retrieve.

        Raises:
            ValueError: A ValueError is raised if user selections are
                not valid.
        """
        if receptor_type in hmbc.genes:
            genes = hmbc.genes[receptor_type]
        else:
            raise ValueError("Unrecognized receptor type supplied. "
                    "Should be one of 'ALL', 'IG', 'TR'.")

        if isinstance(species, str):
            if species == "ALL":
                selected_species = hmbc.allowed_species[receptor_type]
            else:
                raise ValueError("Unexpected species supplied.")
        elif isinstance(species, list):
            for species_name in species:
                if species_name not in hmbc.allowed_species[receptor_type]:
                    raise ValueError(f"For {receptor_type}, only species "
                        f"{hmbc.allowed_species[receptor_type]} are allowed.")
            selected_species = species

        return genes, selected_species


    @staticmethod
    def build_hmm(output_path, species="ALL", receptor_type="ALL"):
        """Constructs a pressed HMM profile file using the IMGT germline data
        for the specified species and chain.

        output_path (str): Filename of a folder where the output HMM
            will be constructed. Should be a directory.
        species (str): Either "ALL" or a list containing one or more of
            the following: "Homo+sapiens", "Mus",
            "Rattus+norvegicus", "Oryctolagus+cuniculus",
            "Macaca+mulatta", "Sus+scrofa", "Vicugna+pacos",
            "Bos+taurus". If "ALL", all species are included in
            the pressed profile file.
        receptor_type (str): One of "ALL", "TR", "IG".
        """
        current_dir = os.getcwd()
        try:
            os.chdir(output_path)
            os.chdir(current_dir)
        except Exception as exc:
            raise ValueError("Invalid output file path supplied.") from exc

        genes, selected_species = IMGTBuilder.validate_inputs(species, receptor_type)
        for gene in genes:
            for target_species in selected_species:
                if gene in hmbc.genes["TR"] and target_species not in hmbc.allowed_species["TR"]:
                    continue
                #A slightly weird quirk of biology -- alpacas lack light chains
                if gene[0] in "KL" and target_species == "Vicugna+pacos":
                    continue
                if retrieve_fasta(target_species, gene):
                    print(f"Error on {species}, {gene}")
                else:
                    print(f"Downloaded {species}, {gene}")





# Html parser class.
class GENEDBParser(HTMLParser):
    """Subclass of HTMLParser.
    """
    currenttag = None
    currentnamedent = None
    _data = []
    def handle_starttag(self, tag, attrs):
        self.currenttag=tag
    def handle_endtag(self, tag):
        self.currenttag=None
    def handle_data(self, data):
        split = data.split("\n")
        start = sum([ 1 if l[0]==">" else 0 for l in split if len(l)])
        if self.currenttag=="pre" and (self.currentnamedent ==">" or start):
            # Two different ways of parsing the html based on how IMGT have formatted the pages.
            # For some reason they format gene db differently sometimes (legacy?) 
            if start > 1: # If you encounter more than one line in the data with a fasta ">" symbol, all sequences will be in the same packet
                name, sequence = None, ""
                for l in split:
                    if not l: continue
                    if l[0]==">":
                        if sequence:
                            self._data.append( (name, sequence) )
                            name, sequence = None, ""
                        name = l
                    else:
                        sequence += l.replace(" ", "")
                if name and sequence:
                    self._data.append( (name, sequence) )
            else: # Otherwise it will be done entry by entry
                print("1")
                try:
                    name = split[0]
                except IndexError:
                    return
                sequence = ("".join( split[1:])).replace(" ", "")
                self._data.append( (name, sequence) )

    def handle_entityref(self, name):
        self.currentnamedent = chr(name2codepoint[name])
    def handle_charref(self, name):
        if name.startswith('x'):
            self.currentnamedent = chr(int(name[1:], 16))
        else:
            self.currentnamedent = chr(int(name))

    def rip_sequences(self,htmlstring):
        """
        Method for this subclass that automates the return of data
        """
        self.reset()
        self._data = []
        self.currenttag = None
        self.currentnamedent = None
        self.feed(htmlstring)
        self._data
        return self._data
        
parser = GENEDBParser()


def get_html(species, gene_type, force = True):
    """
    Get the html from IMGT
    """
    filename = os.path.join(html_outpath,"%s_%s.html"%(species.replace("+", "_"), gene_type) )
    # If html file exists already
    if os.path.isfile(filename):
        return filename
    if urllib.request.urlretrieve( urls[gene_type]%species,  filename ):
        return filename
    else:
        return False

def write_fasta( sequences, species, gene_type ):
    """
    Write a fasta file containing all sequences
    """
    filename = os.path.join(fasta_outpath,"%s_%s.fasta"%(species.replace("+", "_"), gene_type) )
    with open(filename, "w") as outfile:
        for name, sequence in sequences:
            print(">%s"%name, file=outfile)
            print(sequence, file=outfile)

def retrieve_fasta(species, gene_type, force = True):
    """ 
    Retrieve the fasta sequences for a species and gene type from IMGT
    """
    htmlfile = get_html(species, gene_type, force)
    if htmlfile:
        with open(htmlfile) as infile:
            sequences = parser.rip_sequences(infile.read())
        if sequences:
            write_fasta(sequences, species, gene_type )
        else:
            print("Bad parse", end=' ', file=sys.stderr)
            return 1
    else:
        print("Bad Url", end=' ', file=sys.stderr)
        return 1
