#!/usr/bin/python

from Bio.PDB import *
import os, sys

fragment_dir = 'fragments'

def printIntro():
    print "from pymol.cgo import *"
    print "from pymol import cmd"
    print "from pymol.vfont import plain"
    print ""

    print "obj = ["

def printOutro():
    print " "
    print " ]"

    print "cmd.load_cgo(obj, \'ss\')"

def dot(a, b):
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]

def cross(a, b):
    c = [a[1]*b[2] - a[2]*b[1],
         a[2]*b[0] - a[0]*b[2],
         a[0]*b[1] - a[1]*b[0]]

    return c

def ouputFragment(structure, pdb_name, current_model, prev_chain1, prev_chain2, helixStart1, helixStart2, prev_res1, prev_res2):
    #print >>sys.stderr, pdb_name, prev_chain1, prev_chain2, helixStart1, prev_res1, helixStart2, prev_res2
    start1 = min(helixStart1, prev_res1)
    end1 = max(helixStart1, prev_res1)
    
    start2 = min(helixStart2, prev_res2)
    end2 = max(helixStart2, prev_res2)
    
    assert(end1 - start1 == end2 - start2)
    fragment_length = end1 - start1 + 1

    model = structure[current_model]
    chain = model[prev_chain1]
    res = chain[start1]
    atom = 'C1*'

    if abs(helixStart1 - prev_res1) >= 1:
        res_start_1 = structure[current_model][prev_chain1][start1][atom]
        res_start_2 = structure[current_model][prev_chain2][start2][atom]
        
        res_end_1 = structure[current_model][prev_chain1][end1][atom]
        res_end_2 = structure[current_model][prev_chain2][end2][atom]

        next_res_1 = structure[current_model][prev_chain1][start1+1][atom]
        next_res_2 = structure[current_model][prev_chain2][start2+1][atom]
        
        prev_res_1 = structure[current_model][prev_chain1][end1-1][atom]
        prev_res_2 = structure[current_model][prev_chain2][end2-1][atom]

       
        vs1 = res_start_1.get_vector()
        ve1 = res_end_1.get_vector()

        vs2 = res_start_2.get_vector()
        ve2 = res_end_2.get_vector()
     
        n1 = next_res_1.get_vector()
        p1 = prev_res_1.get_vector()

        n2 = next_res_2.get_vector()
        p2 = prev_res_2.get_vector()
     
        mid1 = (vs1 + ve2) / 2.0
        mid2 = (ve1 + vs2) / 2.0

        start_vec1 = structure[current_model][prev_chain1][helixStart1][atom].get_vector() - structure[current_model][prev_chain2][helixStart2][atom].get_vector()
        end_vec1 = structure[current_model][prev_chain1][prev_res1][atom].get_vector() - structure[current_model][prev_chain2][prev_res2][atom].get_vector()

        start_vec2 = structure[current_model][prev_chain1][helixStart1+1][atom].get_vector() - structure[current_model][prev_chain2][helixStart2-1][atom].get_vector()
        end_vec2 = structure[current_model][prev_chain1][prev_res1-1][atom].get_vector() - structure[current_model][prev_chain2][prev_res2+1][atom].get_vector()


        #print "start_vec:", start_vec
        #print "end_vec:", end_vec

        start_norm_vec = Vector(cross(start_vec1, start_vec2))
        start_norm_vec.normalize()

        end_norm_vec = Vector(cross(end_vec2, end_vec1))
        end_norm_vec.normalize()

        start_vec1 = -start_vec1
        end_vec1 = -end_vec1

        start_axis_vec = start_vec1 + Vector([0.,0.,0.])
        start_axis_vec.normalize()
        
        end_axis_vec = end_vec1 + Vector([0.,0.,0.])
        end_axis_vec.normalize()

        start_origin = structure[current_model][prev_chain1][helixStart1][atom].get_vector()
        end_origin = structure[current_model][prev_chain1][prev_res1][atom].get_vector()

        start_x_norm =  start_origin + start_norm_vec + start_norm_vec + start_norm_vec
        start_y_norm = start_origin + start_axis_vec + start_axis_vec + start_axis_vec

        end_x_norm =  end_origin + end_norm_vec + end_norm_vec + end_norm_vec
        end_y_norm = end_origin + end_axis_vec + end_axis_vec + end_axis_vec

        start_y_vec = Vector(cross(start_norm_vec, start_axis_vec))
        start_y_vec.normalize()
        start_c_vec = (start_axis_vec + start_y_vec) / 2
        start_c_vec.normalize()
        start_c_norm = start_origin + start_c_vec / (1 / 8.04)

        end_y_vec = Vector(cross(end_norm_vec, end_axis_vec))
        end_y_vec.normalize()
        end_c_vec = (end_axis_vec + end_y_vec) / 2
        end_c_vec.normalize()
        end_c_norm = end_origin + end_c_vec / (1 / 8.04)


        print " CYLINDER, %f, %f, %f, %f, %f, %f, 0.3, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0," % (start_origin[0], start_origin[1], start_origin[2], start_x_norm[0], start_x_norm[1], start_x_norm[2])
        print " CYLINDER, %f, %f, %f, %f, %f, %f, 0.3, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0," % (start_origin[0], start_origin[1], start_origin[2], start_y_norm[0], start_y_norm[1], start_y_norm[2])
        print " CYLINDER, %f, %f, %f, %f, %f, %f, 0.3, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0," % (start_origin[0], start_origin[1], start_origin[2], start_c_norm[0], start_c_norm[1], start_c_norm[2])

        print " CYLINDER, %f, %f, %f, %f, %f, %f, 0.3, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0," % (end_origin[0], end_origin[1], end_origin[2], end_x_norm[0], end_x_norm[1], end_x_norm[2])
        print " CYLINDER, %f, %f, %f, %f, %f, %f, 0.3, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0," % (end_origin[0], end_origin[1], end_origin[2], end_y_norm[0], end_y_norm[1], end_y_norm[2])
        print " CYLINDER, %f, %f, %f, %f, %f, %f, 0.3, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0," % (end_origin[0], end_origin[1], end_origin[2], end_c_norm[0], end_c_norm[1], end_c_norm[2])

        mid1 = start_c_norm
        mid2 = end_c_norm


        #print norm_vec

        """
        print " CYLINDER, %f, %f, %f, %f, %f, %f, 0.1, 1.0, 1.0, 1.0, 1.0, 0., 0.," % (vs1[0], vs1[1], vs1[2], ve1[0], ve1[1], ve1[2])
        print " CYLINDER, %f, %f, %f, %f, %f, %f, 0.1, 1.0, 1.0, 1.0, 1.0, 0., 0.," % (vs2[0], vs2[1], vs2[2], ve2[0], ve2[1], ve2[2])
        """

        print " CYLINDER, %f, %f, %f, %f, %f, %f, 1.8, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0," % (mid1[0], mid1[1], mid1[2], mid2[0], mid2[1], mid2[2])

    return

def parseChainBase(chainBase):
    #print >>sys.stderr,"chainBase:", chainBase
    if ord(chainBase[0]) >= ord('A') and ord(chainBase[0]) <= ord('Z'):
        chain = chainBase[0]
        base = int(chainBase[1:])
    else:
        if chainBase[0] == '\'':
            endQuoteIdx = chainBase.find('\'',1)
            chain=chainBase[1:endQuoteIdx]
            base=int(chainBase[endQuoteIdx+1:])
        else:
            chain=''
            base=int(chainBase)

    return (chain, base)

def parseBasePairID(basePairID):
    parts = basePairID.split('-')
    (fromChain, fromBase) = parseChainBase(parts[0].strip())
    (toChain, toBase) = parseChainBase(parts[1].strip())

    return (fromChain, fromBase, toChain, toBase)

def printDistances(pdb_file, base_pairs_file):
    p = PDBParser()
    structure=p.get_structure('temp', pdb_file)
    bp_file = open(base_pairs_file, 'r')
    chains = list(structure.get_chains())
    pdb_name = os.path.basename(pdb_file)
    pdb_name = os.path.splitext(pdb_name)[0]
    residues = Selection.unfold_entities(structure, 'R')
    residueCount = len([r for r in residues if r.get_id()[0] == ' '])

    if len(chains) == 0:
        print >>sys.stderr,"Skipping (no chains):", pdb_name
        return
    
    lines = bp_file.readlines()
    bps = []
    base_pair_line = False
    current_model = -1
    printIntro()

    for line in lines:
        # skip all of the lines that do not describe base pair interactions
        if line.find("Base-pairs ---") == 0:
            base_pair_line = True
            current_model += 1
            bp_found = set()
            continue
        if line.find("Residue conformations") == 0:
            base_pair_line = False
            continue
        if not base_pair_line:
            continue

        #print "model:", current_model, "line:", line.strip()
        line_parts = line.split(':')
        (fromChain, fromBase, toChain, toBase) = parseBasePairID(line_parts[0].strip())

        bond_type = line_parts[1].strip()

        #print "fromChain:", fromChain, "fromBase:", fromBase, "toChain:", toChain, "toBase:", toBase, "bond_type:", bond_type
        bps += [(current_model, fromChain, fromBase, toChain, toBase, bond_type)]

    prev_res1 = -100
    prev_res2 = -100

    prev_chain1 = ''
    prev_chain2 = ''
    
    prev_model = -100
    
    helixStart1 = -100
    helixStart2 = -100

    in_fragment = False

    for i in range(len(bps)):
        (current_model, fromChain, fromBase, toChain, toBase, bond_type) = bps[i]

        #   print "current_model:", current_model, "fromChain:", fromChain, "fromBase:", fromBase, "toChain:", toChain, "toBase:", toBase, "bond_type:", bond_type

        if bond_type.find('Ww/Ww') < 0 and bond_type.find('Wh/Ww') < 0 and bond_type.find('Ws/Ww') < 0 and bond_type.find('Ww/Ws') < 0:
            continue
            
        if in_fragment:
            if current_model == prev_model and fromChain == prev_chain1 and toChain == prev_chain2 and abs(fromBase - prev_res1) == 1 and abs(toBase - prev_res2) == 1:
                prev_res1 = fromBase
                prev_res2 = toBase
                continue
            else:
                ouputFragment(structure, pdb_name, current_model, prev_chain1, prev_chain2, helixStart1, helixStart2, prev_res1, prev_res2)
                
        in_fragment = True
        helixStart1 = fromBase
        helixStart2 = toBase
        prev_model = current_model
        prev_res1 = fromBase
        prev_res2 = toBase
        prev_chain1 = fromChain
        prev_chain2 = toChain
            
    if in_fragment:
        ouputFragment(structure, pdb_name, current_model, prev_chain1, prev_chain2, helixStart1, helixStart2, prev_res1, prev_res2)

    printOutro()
            
BASE_PAIR_DIR = 'base-pairs-mcannotate/'

def getBasePairName(pdb_file):
    pdb_name = os.path.basename(pdb_file)
    pdb_name = os.path.splitext(pdb_name)[0]

    base_pair_name = os.path.join(BASE_PAIR_DIR, pdb_name + ".mcannotate")
    print >>sys.stderr,"base_pair_name:", base_pair_name
    return base_pair_name

if __name__ == "__main__":
    if len(sys.argv) == 3:
        printDistances(sys.argv[1], sys.argv[2])
    else:
        print >>sys.stderr, "PDB name:", sys.argv[1]
        basePairName = getBasePairName(sys.argv[1])
        printDistances(sys.argv[1], basePairName)
else:
    structure = '../structures/2RP1.pdb'
    printDistances(structure, getBasePairName(structure))
