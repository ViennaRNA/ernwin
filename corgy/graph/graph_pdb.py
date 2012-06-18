#!/usr/bin/python

import sys, math
from Bio.PDB import *

from corgy.utilities.vector import dot, cross

def get_centroid(chain, residue_num):
    residue_num = [int(i) for i in residue_num]
    #print >>sys.stderr, "residue_num:", residue_num
    atoms = []
    for i in residue_num:
        try:
            atoms += [chain[i]['C1*']]
        except KeyError:
            # the C1* atom probably doesn't exist
            continue

    vectors = [atom.get_vector() for atom in atoms]

    sum_vec = Vector([0., 0., 0.])
    for i in range(len(vectors)):
        sum_vec += vectors[i]

    sum_vec /= float(len(vectors))
    for i in sum_vec:
        if math.isnan(i):
            raise Exception('nan encountered')
    return sum_vec

def get_bulge_centroid(chain, define):
    i = 0
    res_nums = []
    while i < len(define):
        res_nums += range(int(define[i]), int(define[i+1])+1)
        i += 2

    #print >>sys.stderr, "res_nums:", res_nums
    return get_centroid(chain, res_nums)

def get_mids_core(chain, start1, start2, end1, end2):
    #assert(abs(end1 - start1) == abs(end2 - start2))

    fragment_length = end1 - start1 + 1

    #print "start1:", start1, "end1:", end1, "fragment_length:", fragment_length
    if fragment_length < 2:
        raise Exception("Helix shorter than 1 nucleotide")

    catom_name = 'C1*'

    start_vec1 = chain[start1][catom_name].get_vector() - chain[start2][catom_name].get_vector()
    end_vec1 = chain[end1][catom_name].get_vector() - chain[end2][catom_name].get_vector()

    start_vec2 = chain[start1+1][catom_name].get_vector() - chain[start2-1][catom_name].get_vector()
    end_vec2 = chain[end1-1][catom_name].get_vector() - chain[end2+1][catom_name].get_vector()


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

    start_origin = chain[start1][catom_name].get_vector()
    end_origin = chain[end1][catom_name].get_vector()

    start_x_norm =  start_origin + start_norm_vec + start_norm_vec + start_norm_vec
    start_y_norm = start_origin + start_axis_vec + start_axis_vec + start_axis_vec

    end_x_norm =  end_origin + end_norm_vec + end_norm_vec + end_norm_vec
    end_y_norm = end_origin + end_axis_vec + end_axis_vec + end_axis_vec

    start_y_vec = Vector(cross(start_norm_vec, start_axis_vec))
    start_y_vec.normalize()
    start_c_vec = (start_axis_vec + start_y_vec) / 2
    start_c_vec.normalize()
    start_c_norm = start_origin + start_c_vec / (1 / 8.4)

    end_y_vec = Vector(cross(end_norm_vec, end_axis_vec))
    end_y_vec.normalize()
    end_c_vec = (end_axis_vec + end_y_vec) / 2
    end_c_vec.normalize()
    end_c_norm = end_origin + end_c_vec / (1 / 8.4)


    mid1 = start_c_norm
    mid2 = end_c_norm


    #print " CYLINDER, %f, %f, %f, %f, %f, %f, 1.8, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0," % (mid1[0], mid1[1], mid1[2], mid2[0], mid2[1], mid2[2])
    helix_vector = mid2 - mid1
    #helix_vector.normalize()

    return (mid1, mid2)

def get_mids(chain, define):
    return get_mids_core(chain, int(define[0]), int(define[3]), int(define[1]), int(define[2]))

def get_helix_vector(chain, start1, start2, end1, end2):
    (mid1, mid2) = get_mids(chain, start1, start2, end1, end2)
    return mid2 - mid1

def angle_between_helices(chain, stem1, stem2):
    try:
        vec1 = get_helix_vector(chain, stem1[0][0], stem1[0][1], stem1[1][0], stem1[1][1])
        vec2 = get_helix_vector(chain, stem2[0][0], stem2[0][1], stem2[1][0], stem2[1][1])

        vec3 = get_helix_vector(chain, stem1[1][0], stem1[1][1], stem2[0][0], stem2[0][1])
    except Exception as e:
        raise e
        #print e
        return

    #print "vec1:", vec1, dot(vec1, vec1)
    #print "vec2:", vec2, dot(vec2, vec2)
    #print "vec1 * vec2:", dot(vec1, vec2)

    vec1.normalize()
    vec2.normalize()

    angle = math.acos(dot(vec1, vec2)) * (180.0 / math.pi)
    print "stem1:", stem1, "stem2:", stem2
    print "pymol:\n", "color yellow, resi %d-%d or resi %d-%d" % (stem1[1][0]+1, stem2[0][0]-1, stem2[0][1]+1, stem1[1][1]-1)
    #print "color yellow, resi %d or resi %d" % (stem1[1][0], stem1[1][1])
    #print "pymol:", "color yellow, res %d-%d or resi %d-%d or resi %d-%d" %
    print "angle %d %d %.1f" % (stem2[0][0] - stem1[1][0] - 1, stem1[1][1] - stem2[0][1] - 1, angle)
    print "offset %d %d %.2f" % (stem2[0][0] - stem1[1][0] - 1, stem1[1][1] - stem2[0][1] - 1, math.sqrt(dot(vec3, vec3)))


def main():
    if len(sys.argv) < 10:
        print "Usage: ./helix_angle.py pdb_name sm1st1 sm1st2 sm1ed1 sm1ed2 sm2st1 sm2st2 sm2ed1 sm2ed2"
        print
        print "Where sm means stem, st means start and ed means ed. Essentially, the parameters indicate"
        print "at which residues the stems start and end"
        sys.exit(1)
    
    pdb_name = sys.argv[1]

    stem1_start = (int(sys.argv[2]), int(sys.argv[3]))
    stem1_end = (int(sys.argv[4]), int(sys.argv[5]))

    stem1 = (stem1_start, stem1_end)

    stem2_start = (int(sys.argv[6]), int(sys.argv[7]))
    stem2_end = (int(sys.argv[8]), int(sys.argv[9]))

    stem2 = (stem2_start, stem2_end)
    
    #print "pdb_file", pdb_name
    #print "stem1:", stem1
    #print "stem2:", stem2

    s = PDBParser().get_structure('temp', pdb_name)
    c = list(s.get_chains())[0]
    angle = angle_between_helices(c, stem1, stem2)

if __name__=="__main__":
    main()
