"""
Constructs a simple html parser since the IMGT website returns all requested
germline databases in HTML format. The parser extracts the fasta entries
from the HTML returned by requests.
"""

from html.parser import HTMLParser


class IMGT_HTML_Parser(HTMLParser):
    """HTMLParser subclass for parsing the HTML output from requests."""

    def __init__(self, convert_charrefs = True):
        super().__init__(convert_charrefs = convert_charrefs)
        self.current_tag = None
        self._data = []

    def handle_starttag(self, tag, attrs):
        self.current_tag=tag

    def handle_endtag(self, tag):
        self.current_tag=None

    def handle_data(self, data):
        lines = [l for l in data.split("\n") if len(l) > 0]
        num_recs = len([l[0] for l in lines if l[0] == ">"])
        if self.current_tag=="pre" and num_recs > 0 and len(lines) > 0:
            if num_recs > 1:
                name, sequence = None, ""
                for line in lines:
                    if not line:
                        continue
                    if line[0]==">":
                        if sequence:
                            self._data.append( (name, sequence) )
                            sequence = ""
                        name = line
                    else:
                        sequence += line.replace(" ", "")
                if name and sequence:
                    self._data.append( (name, sequence) )
            else: #Precaution for (not often encountered) older format.
                name = lines[0]
                sequence = ("".join( lines[1:])).replace(" ", "")
                self._data.append( (name, sequence) )

    def parse_imgt_html(self,htmlstring):
        """Retrieves the sequences from the html output from the IMGT db."""
        self.reset()
        self._data = []
        self.current_tag = None
        self.feed(htmlstring)
        return self._data
