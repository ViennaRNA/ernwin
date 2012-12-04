
def trim_chain(chain, start_res, end_res):
    '''
    Remove all residues that are not between start_res and end_res.
    '''
    to_detach = []
    for res in chain:
        if res.id[1] < start_res+1 or end_res-1 < res.id[1]:
            to_detach += [res]

    for res in to_detach:
        chain.detach_child(res.id)
