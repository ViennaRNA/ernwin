import Bio.PDB as bpdb
import corgy.utilities.debug as cud
import corgy.utilities.vector as cuv

backbone_atoms = ['P', 'O5*', 'C5*', 'C4*', 'C3*', 'O3*']
ring_atoms = ['C4*', 'C3*', 'C2*', 'C1*', 'O4*']

side_chain_atoms = dict()
side_chain_atoms['U'] = ['N1', 'C2', 'O2', 'N3', 'C4', 'O4', 'C5', 'C6']
side_chain_atoms['C'] = ['N1', 'C2', 'O2', 'N3', 'C4', 'N4', 'C5', 'C6']

side_chain_atoms['A'] = ['N1', 'C2', 'N3', 'C4', 'C5', 'C6', 'N6', 'N7', 'C8', 'N9']
side_chain_atoms['G'] = ['N1', 'C2', 'N2', 'N3', 'C4', 'C5', 'C6', 'O6', 'N7', 'C8', 'N9']

all_rna_atoms = backbone_atoms + ring_atoms
for v in side_chain_atoms.values():
    all_rna_atoms += v
all_rna_atoms = set(all_rna_atoms)

interactions = [('P', 'O5*'),
                ('P', 'OP1'),
                ('P', 'O1P'),
                ('P', 'OP2'),
                ('P', 'O2P'),
                ('C2*', 'O2*'),
                       ('O5*', 'C5*'),
                       ('C5*', 'C4*'),
                       ('C4*', 'O4*'),
                       ('C4*', 'C3*'),
                       ('O4*', 'C1*'),
                       ('C3*', 'C2*'),
                       ('C3*', 'O3*'),
                       ('C2*', 'C1*'),
                       ('C1*', 'N1'),
                       ('N1', 'C2'),
                       ('N1', 'C6'),
                       ('C6', 'C5'),
                       ('C5', 'C4'),
                       ('C4', 'O4'),
                       ('C4', 'N4'),
                       ('C4', 'N3'),
                       ('N3', 'C2'),
                       ('C2', 'O2'),
                       ('C2', 'N2'),
                       ('C1*', 'N9'),
                       ('N9', 'C8'),
                       ('N9', 'C4'),
                       ('C8', 'N7'),
                       ('N7', 'C5'),
                       ('C6', 'O6'),
                       ('C6', 'N6')]

interactions_set = [tuple(sorted(i)) for i in interactions]

def trim_chain(chain, start_res, end_res):
    '''
    Remove all residues that are not between start_res and end_res, inclusive.
    '''
    to_detach = []
    for res in chain:
        if res.id[1] <= start_res or end_res <= res.id[1]:
            to_detach += [res]

    for res in to_detach:
        chain.detach_child(res.id)

def is_covalent(contact):
    '''
    Determine if a particular contact is covalent.

    @param contact: A pair of two Atom objects
    @return True if they are covalently bonded
            False otherwise
    '''
    r1 = contact[0].parent
    r2 = contact[1].parent

    r1a = (r1, contact[0])
    r2a = (r2, contact[1])

    if contact[0].name.find('H') >= 0 or contact[1].name.find('H') >= 0:
        return True

    ((r1, c1), (r2, c2)) = sorted((r1a, r2a), key=lambda x: x[0].id[1])
    
    if r1.id == r2.id:
        if tuple(sorted((c1.name, c2.name))) in interactions_set:
            return True

    if r2.id[1] - r1.id[1] == 1:
        #neighboring residues
        if c1.name == 'O3*' and c2.name == 'P':
            return True

    #cud.pv('((r1.id[1], c1.name), (r2.id[1], c2.name))')
    #cud.pv('r2.id[1] - r1.id[1]')

    return False

def num_noncovalent_clashes(chain):
    '''
    Check if a chain has non-covalent clashes. Non-covalent clashes are found
    when two atoms that aren't covalently linked are within 1.8 A of each other.

    @param chain: The chain to evaluate
    @param return: The number of non-covalent clashes.
    '''
    all_atoms = bpdb.Selection.unfold_entities(chain, 'A')
    ns = bpdb.NeighborSearch(all_atoms)

    contacts = ns.search_all(1.9)

    return len([c for c in contacts if not is_covalent(c)])

def noncovalent_distances(chain, cutoff=0.3):
    '''
    Print out the distances between all non-covalently bonded atoms
    which are closer than cutoff to each other.

    @param chain: The Bio.PDB chain.
    @param cutoff: The maximum distance
    '''
    all_atoms = bpdb.Selection.unfold_entities(chain, 'A')
    ns = bpdb.NeighborSearch(all_atoms)

    contacts = ns.search_all(cutoff)

    return [cuv.magnitude(c[1] - c[0]) for c in contacts if not is_covalent(c)]

