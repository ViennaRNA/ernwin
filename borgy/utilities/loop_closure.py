class BarnacleEnergy:
    def __init__(self):
        self.prev_r = 1000.
        self.min_r = 1000.
        self.min_contacts = 1000.
    
def build_loop(stem_chain, loop_seq, (a,b,i1,i2), iterations):
    '''
    Build a loop.

    The alignment vectors of the nucleotides a and b in stem_chain
    should align to the alignment vectors of nucleotides i1 and i2
    of the sequence sampled by loop_seq.

    @param stem_chain: A chain containing the assembled stems.
    @param loop_seq: The sequence of the loop region including the
                     adjacent stem segments
    @param (a,b,i1,i2): The numbers of the nucleotides defining
                     where the loop starts in the whole sequence (a and b)
                     and within the loop sequence (i1 and i2)
    @param iterations: The number of MCMC iterations to run

    @return: A Bio.PDB.Chain structure containing the best sampled loop.
    '''
    model = barn.BarnacleCPDB(seq, 2.)



