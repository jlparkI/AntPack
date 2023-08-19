"""Contains the DomainAlignment class, which stores information about
a particular domain alignment."""
from .constants import allowed_inputs
from .germlines import all_germlines
from .numbering_funcs import number_imgt, number_aho

all_reference_states = list(range( 1, 129)) # These are the IMGT reference states (matches)
all_species = list(all_germlines['V']['H'].keys())



class DomainAligment:
    """This class contains information about a particular domain
    alignment.
    """

    def __init__(self, hit, domain, order, ndomains, full_sequence, scheme="imgt",
            assign_germline=False):
        """Class constructor.

        Args:
        """
        self.evalue = hit.evalue
        self.bitscore = hit.score
        self.bias = hit.bias
        self.hit_id = hit.name.decode()
        self.species, self.chain = self.hit_id.split("_")
        self.query_start = domain.alignment.target_from - 1
        self.query_end = domain.alignment.target_to
        self.hmm_start = domain.alignment.hmm_from
        self.hmm_end = domain.alignment.hmm_to

        self.domain = domain
        self.order = order
        self.scheme = scheme
        state_vector = self._alignment_to_states(domain.alignment, ndomains,
                len(full_sequence))

        # Check the numbering for likely very long CDR3s that will have been
        # missed by the first pass. Modify alignments in-place
        self._check_j_region(full_sequence, domain.alignment, scheme)

        # Apply the desired numbering scheme to all sequences
        self.numbering = self._validate_numbering(self.number_sequence_from_alignment(state_vector,
                    full_sequence), full_sequence)


    def _alignment_to_states(self, alignment, ndomains, seq_length):
        """
        Converts an alignment from pyhmmer into a state vector that can be used for numbering the
        sequence.
        """
        reference_string = alignment.hmm_sequence
        state_string = alignment.target_sequence

        assert len(reference_string) == len(state_string), "Aligned reference and state strings had different lengths. Don't know how to handle"

        _hmm_length = self.get_hmm_length(self.species, self.chain)

        # Handle cases where there are n terminal modifications.
        # In most cases the user is going to want these included in the numbered domain even though they are not 'antibody like' and 
        # not matched to the germline. Only allow up to a maximum of 5 unmatched states at the start of the domain
        # Adds a bug here if there is a very short linker between a scfv domains with a modified n-term second domain
        # Thus this is only done for the first identified domain ( hence order attribute on hsp )
        if self.order == 0 and (0 < self.hmm_start < 5):
            n_extend = self.hmm_start
            if self.hmm_start > self.query_start:
                n_extend = min(self.query_start, self.hmm_start - self.query_start)
            state_string = 'Z'*n_extend + state_string
            reference_string = 'Z'*n_extend + reference_string
            self.query_start -= n_extend
            self.hmm_start -= n_extend

        # Handle cases where the alignment should be extended to the end of the j-element
        # This occurs when there a c-terminal modifications of the variable domain that are significantly different to germline
        # Extension is only made when half of framework 4 has been recognised and there is only one domain recognised.
        if ndomains==1 and self.query_end < seq_length and (123 < self.hmm_end < _hmm_length): # Extend forwards
            n_extend = min( _hmm_length - self.hmm_end, seq_length - self.query_end )
            state_string = state_string + 'Z'*n_extend
            reference_string = reference_string + 'Z'*n_extend
            self.query_end += n_extend
            self.hmm_end += n_extend
                    


        # Generate lists for the states and the sequence indices that are included in this alignment
        hmm_states = all_reference_states[ self.hmm_start : self.hmm_end ]
        sequence_indices = list(range(self.query_start,  self.query_end))
        h, s = 0, 0 # initialise the current index in the hmm and the sequence
    
        state_vector = []
        # iterate over the state string (or the reference string)
        for hmm_state, query_state in zip(reference_string, state_string):
            if hmm_state == ".":
                state_type = "i"
            else:
                state_type = "m"
            if query_state == "-":
                state_type = "d"
                sequence_index = None
            else:
                sequence_index = sequence_indices[s]

            state_vector.append(  ((hmm_states[h], state_type),  sequence_index )  )

            # Updates to the indices
            if state_type == "m":
                h+=1
                s+=1
            elif state_type == "i":
                s+=1
            else: # delete state
                h+=1

        return state_vector


    def _validate_numbering(self, numbering, seq):
        """Checks that the numbering does not exhibit obvious errors."""
        last, nseq = -1, ""
        for (index, _), a in numbering:
            assert index >= last, "Numbering was found to decrease along the sequence. Please report."
            last = index
            nseq += a.replace("-","")

        assert nseq in seq.replace("-",""), "The algorithm did not number a contiguous segment for sequence. Please report"

        return numbering



    def get_hmm_length(self, species, ctype ):
        '''
        Get the length of an hmm given a species and chain type. 
        This tells us how many non-insertion positions there could possibly be in a domain (127 or 128 positions under imgt)
        '''
        try:
            return len(list(all_germlines['J'][ctype][species].values())[0].rstrip('-'))
        except KeyError:
            return 128



    def number_sequence_from_alignment(self, state_vector, sequence):
        """Generates numbering for an alignment using the state vector and sequence.

        Args:
            state_vector (list): A list of states (insertion, match etc.) from the initial
                alignment; needs further "tweaking" to generate the IMGT numbering.
            sequence (str): The original sequence.

        Returns:
            numbering (list): A list of tuples of (position id, aa).
        """
        if self.scheme == "imgt":
            return number_imgt(state_vector, sequence)
        elif self.scheme == "aho":
            return number_aho(state_vector, sequence, self.chain)

        scheme_chain = (self.scheme, self.chain)
        if scheme_chain in allowed_inputs.scheme_to_fun:
            return allowed_inputs.scheme_to_fun[scheme_chain](state_vector, sequence)
        raise AssertionError(f"The numbering scheme {self.scheme} is not "
                f"implemented for chain {self.chain}.")






def get_identity( state_sequence, germline_sequence ):
    """
    Get the partially matched sequence identity between two aligned sequences. 
    Partial in the sense that gaps can be in the state_sequence.
    """
    # Ensure that the sequences are the expected length
    assert len( state_sequence) == len(germline_sequence ) == 128
    n, m = 0, 0
    for i in range( 128 ):
        if germline_sequence[i] == "-":continue
        if state_sequence[i].upper() == germline_sequence[i]: m+=1
        n+=1

    if not n:
        return 0    
    return float(m)/n


def run_germline_assignment(state_vector, sequence, chain_type, allowed_species=None ):
    """
    Find the closest sequence identity match.
    """
    genes={'v_gene': [None,None],
           'j_gene': [None,None],
         }


    # Extract the positions that correspond to match (germline) states. 
    state_dict = dict( ((i, 'm'),None) for i in range(1,129))
    state_dict.update(dict(state_vector))
    state_sequence = "".join([ sequence[state_dict[(i, 'm')]] if state_dict[(i,'m')] is not None else "-" for i in range(1,129) ])

    # Iterate over the v-germline sequences of the chain type of interest.
    # The maximum sequence identity is used to assign the germline 
    if chain_type in all_germlines["V"]:
        if allowed_species is not None:
            if not all( [ sp in all_germlines['V'][chain_type] for sp in allowed_species ] ): # Made non-fatal
                return {}
        else:
            allowed_species = all_species
        seq_ids = {}
        for species in allowed_species:
            if species not in all_germlines["V"][ chain_type ]: continue # Previously bug.
            for gene, germline_sequence in all_germlines["V"][ chain_type ][ species ].items():
                seq_ids[ (species, gene) ] = get_identity( state_sequence , germline_sequence )
        genes['v_gene' ][0] = max( seq_ids, key=lambda x: seq_ids[x] )
        genes['v_gene' ][1] = seq_ids[ genes['v_gene' ][0] ]
        
        # Use the assigned species for the v-gene for the j-gene. 
        # This assumption may affect exotically engineered abs but in general is fair.
        species = genes['v_gene' ][0][0]       
        if chain_type in all_germlines["J"]:
            if species in all_germlines["J"][chain_type]:
                seq_ids = {}
                for gene, germline_sequence in all_germlines["J"][ chain_type ][ species ].items():
                    seq_ids[ (species, gene) ] = get_identity( state_sequence , germline_sequence )
                genes['j_gene' ][0] = max( seq_ids, key=lambda x: seq_ids[x] )
                genes['j_gene' ][1] = seq_ids[ genes['j_gene' ][0] ]
     
    return genes



    def check_j_region( full_sequence, alignment, scheme ):
        '''
        As the length of CDR3 gets long (over 30ish) an alignment that does not include the J region becomes more favourable.
        This leads to really long CDR3s not being numberable. 

        To overcome this problem, when no J region is detected we try without the v region.
        '''
        if len(alignments[1]) ==1: # Only do for single domain chains. 

                # Check whether a J region has been identified. If not check whether there is still a considerable amount of sequence   
                # remaining. 
                ali = alignments[i][1][0]

                # Find the last match position. 
                last_state  = ali[-1][0][0]
                last_si     = ali[-1][1]
                if last_state < 120: # No or very little J region
                    if last_si + 30 < len( sequences[i][1] ): # Considerable amount of sequence left...suspicious of a long CDR3
                        # Find the position of the conserved cysteine (imgt 104). 
                        cys_si = dict( ali ).get( (104,'m'), None )
                        if cys_si is not None: # 104 found.

                            # Find the corresponding index in the alignment.
                            cys_ai = ali.index( ((104, 'm'), cys_si) )
                
                            # Try to identify a J region in the remaining sequence after the 104. A low bit score threshold is used.
                            _, re_states, re_details  = run_hmmer( [(sequences[i][0], sequences[i][1][cys_si+1:])], 
                                                               bit_score_threshold=10 )[0] 

                            # Check if a J region was detected in the remaining sequence.
                            if re_states and re_states[0][-1][0][0] >= 126 and re_states[0][0][0][0] <= 117: 

                                # Sandwich the presumed CDR3 region between the V and J regions.

                                vRegion   = ali[:cys_ai+1]
                                jRegion   = [ (state, index+cys_si+1) for state, index in re_states[0] if state[0] >= 117 ]
                                cdrRegion = []
                                next = 105
                                for si in range( cys_si+1, jRegion[0][1] ):
                                    if next >= 116:
                                        cdrRegion.append( ( (116, 'i'), si ) )
                                    else:
                                        cdrRegion.append( ( (next, 'm'), si ) )
                                        next +=1 

                                # Update the alignment entry.
                                alignments[i][1][0] = vRegion + cdrRegion + jRegion
                                alignments[i][2][0]['query_end'] = jRegion[-1][1] + 1
