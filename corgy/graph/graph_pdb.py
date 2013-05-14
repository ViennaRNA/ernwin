#!/usr/bin/python

import Bio.PDB as bp
import Bio.PDB as bpdb
import sys, operator

import os.path as op
import corgy.builder.config as cbc
import warnings

import math as m
import collections as co

import numpy as np
import numpy.linalg as nl
import numpy.testing as nt
import numpy.random as nr

import corgy.utilities.debug as cud
import corgy.utilities.my_math as cum
import corgy.utilities.pdb as cup
import corgy.utilities.vector as cuv
import corgy.utilities.average_stem_vres_atom_positions as cua

import scipy.optimize as so
import numpy.linalg as nl

catom_name = 'C1*'

def stem_stem_orientation(bg, s1, s2):
    '''
    Calculate the orientation of stem s2 in relation to stem s1
    as described by 3 parameters:

    1. The distance between the closest points of the two stems.
    2. The angle between s1 and s2 in the plane formed by the axis of
       the first stem and the vector between the two points closest
       to each on both stems.
    3. The angle of s2 out of the plane formed by their axes.

    @param bg: The BulgeGraph containing the stems.
    @param s1: The name of the first stem
    @param s2: The name of the second stem
    @return: (x,y,z) where x,y and z are the parameters described in
        the description above.
    '''
    # shorten the names a little bit
    s1_p0 = bg.coords[s1][0]
    s1_p1 = bg.coords[s1][1]
    
    twist1_vec = bg.twists[s1][0]
    twist2_vec = bg.twists[s2][0]

    s2_p0 = bg.coords[s2][0]
    s2_p1 = bg.coords[s2][1]

    # The vectors of the axes of the cylinders
    s1_vec = bg.coords[s1][1] - bg.coords[s1][0]
    s2_vec = bg.coords[s2][1] - bg.coords[s2][0]

    s1_len = cuv.magnitude(s1_vec)
    s2_len = cuv.magnitude(s2_vec)

    stem1_basis = cuv.create_orthonormal_basis(s1_vec, twist1_vec)
    stem2_basis = cuv.create_orthonormal_basis(s2_vec, twist2_vec)

    stem1_origin = bg.coords[s1][0]
    stem2_origin = bg.coords[s2][0]

    # the minimum distance between the two stems, which are represented
    # as line segments
    (i1, i2) = cuv.line_segment_distance(s1_p0, s1_p1, s2_p0, s2_p1)
    i_vec = i2 - i1
    #cud.pv('s1,s2,cuv.magnitude(i2-i1)')

    s2_on_s1_intersect = cuv.change_basis(i2 - stem1_origin, stem1_basis, cuv.standard_basis)
    s1_on_s2_intersect = cuv.change_basis(i1 - stem2_origin, stem2_basis, cuv.standard_basis)

    if s2_on_s1_intersect[0] > 0.:
        if s2_on_s1_intersect[0] > s1_len:
            s2_on_s1_dist = s2_on_s1_intersect[0] - s1_len
        else:
            s2_on_s1_dist = 0.
    else:
        s2_on_s1_dist = abs(s2_on_s1_intersect[0])

    if s1_on_s2_intersect[0] > 0.:
        if s1_on_s2_intersect[0] > s2_len:
            s1_on_s2_dist = s1_on_s2_intersect[0] - s2_len
        else:
            s1_on_s2_dist = 0.
    else:
        s1_on_s2_dist = abs(s1_on_s2_intersect[0])

    lateral_offset = (s2_on_s1_dist + s1_on_s2_dist) / 2.
    ortho_offset = (abs(s1_on_s2_intersect[1]) + abs(s2_on_s1_intersect[1])) / 2.
        
    offset_vec1 = cuv.vec_angle(i_vec, s1_p1 - s1_p0)
    offset_vec2 = cuv.vec_angle(i_vec, s2_p1 - s2_p0)

    offset_vec1 = min(offset_vec1, m.pi - offset_vec1)
    offset_vec2 = min(offset_vec2, m.pi - offset_vec2)

    offset_vec = (offset_vec1 + offset_vec2) / 2.

    '''
    cud.pv('"***************"')
    cud.pv('(s1,s2)')
    cud.pv('offset_vec')
    cud.pv('(s2_on_s1_dist + s1_on_s2_dist) / 2.')
    cud.pv('(lateral_offset, ortho_offset)')
    '''

    # the normal of the plane defined by the two stem vectors
    plane_norm = np.cross(s1_vec, i_vec)

    # the projects of s2_vec onto the plane
    s2_proj = cuv.vector_rejection(s2_vec, plane_norm) 

    # the angle out of the plane
    ang1 = cuv.vec_angle(s2_proj, s2_vec)

    # the angle in the plane
    #cud.pv('s2_proj')
    #cud.pv('s1_vec')
    #cud.pv('np.dot(cuv.normalize(s2_proj), cuv.normalize(s1_vec))')
    ang2 = cuv.vec_angle(s2_proj, s1_vec)

    #ang1 = min(ang1, m.pi - ang1)
    #ang2 = min(ang2, m.pi - ang2)

    return (cuv.magnitude(i_vec), ang1, ang2, cuv.vec_angle(s1_vec, s2_vec), lateral_offset, ortho_offset)

def get_stem_phys_length(coords):
    '''
    Return the physical length of a stem.

    @param coords: The coordinates of the ends of the helix axis.
    '''

    return cuv.magnitude(coords[1] - coords[0])

def base_normals(pdb_filename):
    '''
    Return a list of the normals for each base in the structure.

    As defined by the average of the cross products between the C2-C5 
    and C2-C6 vectors and the N3-C6 and N3-C5 vectors. The origin of 
    the vector will be the centroid of these four atoms.

    @param pdb_filename: The name of the pdb file containing the structure
    @return: A list of pairs containing the origin the normal as well as the 
        normal itself.
    '''
    chain = list(bp.PDBParser().get_structure('t', pdb_filename).get_chains())[0]
    origin_norms = []

    for res in chain:
        c2 = res['C2'].get_vector().get_array()
        c5 = res['C5'].get_vector().get_array()
        c6 = res['C6'].get_vector().get_array()
        n3 = res['N3'].get_vector().get_array()

        v1 = cuv.normalize(np.cross(c6-c2, c5-c2))
        v2 = cuv.normalize(np.cross(c6-n3, c5-n3))

        # take the average of the two, for accuracy or something
        v_norm = (v1 + v2) / 2

        origin = (c2 + c5 + c6 + n3) / 4
        origin_norms += [(origin, v_norm)]

    return origin_norms


def get_twist_angle(coords, twists):
    ''' 
    Get the angle of the twists with respect to each other.

    @param coords: The coordinates of the ends of the stem.
    @param twists: The two twist vectors.
    @return angle: The angle between the two twist vectors.
    '''

    stem_vec = coords[1] - coords[0]
    basis = cuv.create_orthonormal_basis(stem_vec, twists[0])

    twist2 = cuv.change_basis(twists[1], basis, cuv.standard_basis)
    #assert_allclose(twist2[0], 0., rtol=1e-7, atol=1e-7)

    angle = m.atan2(twist2[2], twist2[1])
    return angle

def twist2_from_twist1(stem_vec, twist1, angle):
    '''
    Get an orientation for the second twist which will place it an
    angle of angle from the first twist.

    @param stem_vec: The vector of the stem.
    @param twist1: The vector of the first twist.
    @param angle: The angular difference between the two twists.
    '''
    basis = cuv.create_orthonormal_basis(stem_vec, twist1)

    twist2_new = np.array([0., m.cos(angle), m.sin(angle)])
    twist2 = np.dot(basis.transpose(), twist2_new)
    #twist2 = cuv.change_basis(twist2_new, cuv.standard_basis, basis)

    return twist2

def get_twist_parameter(twist1, twist2, (u, v)):
    '''
    Calculate how much stem1 must be twisted for its twist vector
    to coincide with that of stem2.

    @param twist1: The twist notator of stem1
    @param twist2: The twist notator of stem2
    @param (u,v): The parameters for rotating stem2 onto stem1
    '''

    rot_mat1 = cuv.rotation_matrix(cuv.standard_basis[2], v)
    rot_mat2 = cuv.rotation_matrix(cuv.standard_basis[1], u - m.pi/2)

    twist2_new = np.dot(rot_mat1, twist2)
    twist2_new = np.dot(rot_mat2, twist2_new)
    
    #print "get_twist_parameter twist2:", twist2_new

    return m.atan2(twist2_new[2], twist2_new[1])

def get_stem_orientation_parameters(stem1_vec, twist1, stem2_vec, twist2):
    '''
    Return a parameterization of the orientation of stem2 with respect to stem1.

    stem1 -> bulge -> stem2

    @param stem1_vec: The vector representing the axis of stem1
    @param twist1: The twist of stem1 closest to the bulge
    @param stem2_vec: The vector representing teh axis of stem2
    '''

    # Since we will denote the orientation of stem2 with respect to stem1
    # We first need to define a new coordinate system based on stem1

    stem1_basis = cuv.create_orthonormal_basis(stem1_vec, twist1)

    # Transform the vector of stem2 to the new coordinate system
    stem2_new_basis = cuv.change_basis(stem2_vec, stem1_basis, cuv.standard_basis) 
    twist2_new_basis = cuv.change_basis(twist2, stem1_basis, cuv.standard_basis)

    # Convert the cartesian coordinates to polar coordinates
    (r, u, v) = cuv.spherical_cartesian_to_polar(stem2_new_basis)
    t = get_twist_parameter(twist1, twist2_new_basis, (u, v))

    return (r, u, v, t)


def get_stem_separation_parameters(stem, twist, bulge):
    '''
    Parameterize the location of the bulge with respect to the stem.

    @param stem: The stem vector.
    @param bulge: the bulge vector.
    '''

    stem_basis = cuv.create_orthonormal_basis(stem, twist)
    bulge_new_basis = cuv.change_basis(bulge, stem_basis, cuv.standard_basis)

    return cuv.spherical_cartesian_to_polar(bulge_new_basis)


def get_stem_twist_and_bulge_vecs(bg, bulge, connections):
    '''
    Return the vectors of the stems and of the twists between which
    we want to calculate angles.

    The two vectors will be defined as follows:

    s1e -> s1b -> b -> s2b -> s2e

    The twists will be the two closest to the bulge.

    @param bulge: The name of the bulge separating the two helices.
    @param connections: The two stems that are connected to this bulge.
    @return: (stem1, twist1, stem2, twist2, bulge)
    '''

    s1 = connections[0]
    s2 = connections[1]

    #s1d = bg.defines[s1]
    #s2d = bg.defines[s2]

    mids1 = bg.coords[s1]
    twists1 = bg.twists[s1]

    mids2 = bg.coords[s2]
    twists2 = bg.twists[s2]

    # find out which sides of the stems are closest to the bulge
    # the results will be indexes into the mids array
    (s1b, s1e) = bg.get_sides(s1, bulge)
    (s2b, s2e) = bg.get_sides(s2, bulge)

    # Create directional vectors for the stems
    stem1_vec = mids1[s1b] - mids1[s1e]
    bulge_vec = mids2[s2b] - mids1[s1b]
    stem2_vec = mids2[s2e] - mids2[s2b]

    #twists1_vec = [twists1[s1b], twists1[s1e]]
    #twists2_vec = [twists2[s2e], twists2[s2b]]

    return (stem1_vec, twists1[s1b], stem2_vec, twists2[s2b], bulge_vec)

def stem2_pos_from_stem1(stem1, twist1, params):
    '''
    Get the starting point of a second stem, given the parameters
    about where it's located with respect to stem1

    @param stem1: The vector representing the axis of stem1's cylinder
    @param twist1: The twist parameter of stem1
    @param params: The parameters describing the orientaiton of stem2 wrt stem1
    '''
    (r, u, v) = params
    stem2 = cuv.spherical_polar_to_cartesian((r, u, v))

    stem1_basis = cuv.create_orthonormal_basis(stem1, twist1)
    stem2_start = np.dot(stem1_basis.transpose(), stem2)

    return stem2_start

def stem2_pos_from_stem1_1(stem1_basis, params):
    '''
    Get the starting point of a second stem, given the parameters
    about where it's located with respect to stem1

    @param stem1: The vector representing the axis of stem1's cylinder
    @param twist1: The twist parameter of stem1
    @param params: The parameters describing the orientaiton of stem2 wrt stem1
    '''
    (r, u, v) = params
    stem2 = cuv.spherical_polar_to_cartesian((r, u, v))
    stem2_start = np.dot(stem1_basis, stem2)

    return stem2_start

def twist2_orient_from_stem1(stem1, twist1, (u, v, t)):
    '''
    Calculate the position of the twist factor of the 2nd stem from its
    parameters and the first stem.

    @param stem1: The vector representing the axis of stem1's cylinder
    @param twist1: The twist factor of stem1.
    @param (u, v, t): The parameters describing how the twist of stem2 is oriented with respect to stem1
    '''

    twist2_new = np.array([0., m.cos(t), m.sin(t)])

    rot_mat1 = cuv.rotation_matrix(cuv.standard_basis[2], v)
    rot_mat2 = cuv.rotation_matrix(cuv.standard_basis[1], u - m.pi/2)

    rot_mat = np.dot(rot_mat2, rot_mat1)
    twist2_new = np.dot(nl.inv(rot_mat), twist2_new)
    
    '''
    twist2_new = dot(inv(rot_mat2), twist2_new)
    twist2_new = dot(inv(rot_mat1), twist2_new)
    '''

    stem1_basis = cuv.create_orthonormal_basis(stem1, twist1)
    twist2_new_basis = cuv.change_basis(twist2_new, cuv.standard_basis, stem1_basis)

    return twist2_new_basis

def twist2_orient_from_stem1_1(stem1_basis, (u, v, t)):
    '''
    Calculate the position of the twist factor of the 2nd stem from its
    parameters and the first stem.

    @param stem1: The vector representing the axis of stem1's cylinder
    @param twist1: The twist factor of stem1.
    @param (u, v, t): The parameters describing how the twist of stem2 is oriented with respect to stem1
    '''

    twist2_new = np.array([0., m.cos(t), m.sin(t)])

    rot_mat1 = cuv.rotation_matrix(cuv.standard_basis[2], v)
    rot_mat2 = cuv.rotation_matrix(cuv.standard_basis[1], u - m.pi/2)

    rot_mat = np.dot(rot_mat2, rot_mat1)
    twist2_new = np.dot(nl.inv(rot_mat), twist2_new)

    twist2_new_basis = np.dot(stem1_basis, twist2_new)

    return twist2_new_basis

def stem2_orient_from_stem1(stem1, twist1, (r, u, v)):
    '''
    Calculate the orientation of the second stem, given its parameterization
    and the parameterization of stem1

    @param stem1: The vector representing the axis of stem1's cylinder
    @param twist1: The twist factor of stem1.
    @param (r,u,v): The orientation of stem2 wrt stem1
    '''
    stem2 = cuv.spherical_polar_to_cartesian((r, u, v))
    stem1_basis = cuv.create_orthonormal_basis(stem1, twist1)
    stem2 = cuv.change_basis(stem2, cuv.standard_basis, stem1_basis)

    return stem2

def stem2_orient_from_stem1_1(stem1_basis, (r, u, v)):
    '''
    Calculate the orientation of the second stem, given its parameterization
    and the parameterization of stem1

    @param stem1: The vector representing the axis of stem1's cylinder
    @param twist1: The twist factor of stem1.
    @param (r,u,v): The orientation of stem2 wrt stem1
    '''
    stem2 = cuv.spherical_polar_to_cartesian((r, u, v))
    #stem1_basis = cuv.create_orthonormal_basis(stem1, twist1)
    #stem2 = cuv.change_basis(stem2, cuv.standard_basis, stem1_basis)
    stem2 = np.dot(stem1_basis, stem2)

    return stem2

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

    vectors = [atom.get_vector().get_array() for atom in atoms]

    return cuv.get_vector_centroid(vectors)

def get_bulge_centroid(chain, define):
    i = 0
    res_nums = []
    while i < len(define):
        res_nums += range(int(define[i]), int(define[i+1])+1)
        i += 2

    #print >>sys.stderr, "res_nums:", res_nums
    return get_centroid(chain, res_nums)

def estimate_mids_core(chain, start1, start2, end1, end2):
    '''
    Get the start and end points of a helix.
    '''
    #assert(abs(end1 - start1) == abs(end2 - start2))

    fragment_length = end1 - start1 + 1

    if fragment_length < 2:
        raise Exception("Helix shorter than 1 nucleotide: start1: %d start2: %d end1: %d end2: %d fragment_length: %d" % (start1, start2, end1, end2, fragment_length))

    # get the vector between the CA atoms of the two starting residues
    # as well as for the two next residues
    start_vec1 = chain[start1][catom_name].get_vector() - chain[start2][catom_name].get_vector()
    end_vec1 = chain[end1][catom_name].get_vector() - chain[end2][catom_name].get_vector()

    # get the vector between the CA atoms of the two ending residues
    # as  well as for the two previous residues
    start_vec2 = chain[start1+1][catom_name].get_vector() - chain[start2-1][catom_name].get_vector()
    end_vec2 = chain[end1-1][catom_name].get_vector() - chain[end2+1][catom_name].get_vector()


    # a vector kind of pointing in the direction of the helix
    start_norm_vec = bp.Vector(np.cross(start_vec1.get_array(), start_vec2.get_array()))
    start_norm_vec.normalize()

    # another vector kind of pointing in the directions of the helix
    end_norm_vec = bp.Vector(np.cross(end_vec2.get_array(), end_vec1.get_array()))
    end_norm_vec.normalize()

    start_vec1 = -start_vec1
    end_vec1 = -end_vec1

    # I guess we're converting them to Vector format in a weird sort of way...
    # otherwise these steps don't really make sense
    start_axis_vec = start_vec1 + bp.Vector([0.,0.,0.])
    start_axis_vec.normalize()
    
    end_axis_vec = end_vec1 + bp.Vector([0.,0.,0.])
    end_axis_vec.normalize()

    start_origin = chain[start1][catom_name].get_vector()
    end_origin = chain[end1][catom_name].get_vector()

    start_y_vec = bp.Vector(np.cross(start_norm_vec.get_array(), start_axis_vec.get_array()))
    start_y_vec.normalize()
    start_c_vec = (start_axis_vec + start_y_vec) / 2
    start_c_vec.normalize()
    start_c_norm = start_origin + start_c_vec / (1 / 8.4)

    end_y_vec = bp.Vector(np.cross(end_norm_vec.get_array(), end_axis_vec.get_array()))
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

def basenormals_mids(chain, start1, start2, end1, end2):
    '''
    Calculate a helix axis based on the base normal vectors.

    See Laederach et al., RNA 2007.
    '''

    # calculate the scatter matrix using the base normal vectors
    scatter = np.zeros((3,3))

    residue_numbers = [i for i in range(start1, end1+1)]
    residue_numbers += [i for i in range(end2, start2+1)]

    for r in residue_numbers:
        c2 = chain[r]['C2'].get_vector().get_array()
        c4 = chain[r]['C4'].get_vector().get_array()
        c6 = chain[r]['C6'].get_vector().get_array()

        xi = cuv.normalize(np.cross(c2 - c4, c4 - c6))
        scatter += np.dot(xi, xi.T)

    scatter /= float(len(residue_numbers))

    # compute the eigenvalues and eigenvectors
    w, v = np.linalg.eig(scatter)
    index, value = max(enumerate(w), key=operator.itemgetter(1))

    axis = v[index]

    cud.pv('w,v, index, value')
    cud.pv('axis')

    # estimate the start and end position, which will be converted scaled
    # to the position of the helix axis
    start_pos = (chain[start1][catom_name].get_vector().get_array() + 
                 chain[start2][catom_name].get_vector().get_array()) / 2.

    end_pos = (chain[start1 + real_stem_length][catom_name].get_vector().get_array() + 
                 chain[start2 - real_stem_length][catom_name].get_vector().get_array()) / 2.


def get_mids_core_a(chain, start1, start2, end1, end2, use_template=True):
    '''
    Estimate the stem cylinder using the old method and then refine it 
    using fitted parameters.
    '''
    if use_template:
        template_stem_length=30
    else:
        template_stem_length = end1 - start1 + 1

    real_stem_length = end1 - start1

    tstart1 = 1
    tstart2 = template_stem_length * 2
    tend1 = template_stem_length
    tend2 = template_stem_length + 1

    template_filename = 'ideal_1_%d_%d_%d.pdb' % (tend1, tend2, tstart2)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        #cud.pv('filename')
        ideal_chain = list(bpdb.PDBParser().get_structure('test', 
                op.join(cbc.Configuration.stem_fragment_dir, template_filename)).get_chains())[0]

        ideal_chain = extract_define_residues([tstart1,tend1,tend2,tstart2], ideal_chain)

    est_mids = estimate_mids_core(ideal_chain, tstart1, tstart2, tend1, tend2)
    est_mids = [est_mids[0].get_array(), est_mids[1].get_array()]

    start_pos = (chain[start1][catom_name].get_vector().get_array() + 
                 chain[start2][catom_name].get_vector().get_array()) / 2.

    end_pos = (chain[start1 + real_stem_length][catom_name].get_vector().get_array() + 
                 chain[start2 - real_stem_length][catom_name].get_vector().get_array()) / 2.

    atom_poss = []
    residue_numbers = [i for i in range(tstart1, tend1+1)]
    residue_numbers += [i for i in range(tend2, tstart2+1)]
    
    #cud.pv('residue_numbers')

    for rn in residue_numbers:
        #atom_poss += [chain[rn]['C1*'].get_vector().get_array()]
        try:
            '''
            for atom in chain[rn].get_list():
                atom_poss += [atom.get_vector().get_array()]
            '''

            atom_poss += [ideal_chain[rn]['P'].get_vector().get_array()]
            atom_poss += [ideal_chain[rn]['O3*'].get_vector().get_array()]
            atom_poss += [ideal_chain[rn]['C3*'].get_vector().get_array()]
            atom_poss += [ideal_chain[rn]['C4*'].get_vector().get_array()]
            atom_poss += [ideal_chain[rn]['C5*'].get_vector().get_array()]
            atom_poss += [ideal_chain[rn]['O5*'].get_vector().get_array()]
            atom_poss += [ideal_chain[rn]['C1*'].get_vector().get_array()]
        except KeyError as ke:
            pass
            #cud.pv('ke')

    mids =  fit_circle(est_mids, np.array(atom_poss), start_pos, end_pos)
    return mids

def get_mids_core(chain, start1, start2, end1, end2, use_template=True):
    ######## Debug function
    '''
    vec = mids[1] - mids[0]
    n1 = []
    n2 = []
    for i in range(0, end1-start1+1):
        notch1 = chain[start1+i][catom_name].get_vector().get_array()
        notch2 = chain[start2-i][catom_name].get_vector().get_array()
        basis = cuv.create_orthonormal_basis(vec)
        notch1_n = cuv.change_basis(notch1, basis, cuv.standard_basis)
        notch2_n = cuv.change_basis(notch2, basis, cuv.standard_basis)

        n1 += [notch1_n[0]]
        n2 += [notch2_n[0]]
    #cud.pv('n')
    dists1 = [j-i for i,j in zip(n1[:-1], n1[1:])]
    dists2 = [j-i for i,j in zip(n2[:-1], n2[1:])]

    for dist in dists1 + dists2: 
        print "ladder", dist
    '''
    ######## End debug

    #filename = 
    stem_length = end1 - start1 + 1
    #cud.pv('stem_length')
    filename = 'ideal_1_%d_%d_%d.pdb' % (stem_length, stem_length+1, stem_length*2)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        #cud.pv('filename')
        ideal_chain = list(bpdb.PDBParser().get_structure('test', 
                op.join(cbc.Configuration.stem_fragment_dir, filename)).get_chains())[0]

        chain = extract_define_residues([start1,end1,end2,start2], chain)

    rotran = cup.pdb_rmsd(chain, ideal_chain, sidechains=False, 
            superimpose=True, apply_sup=False)[2]

    ideal_mids = get_mids_core_a(ideal_chain, 1, stem_length*2, stem_length, stem_length+1, use_template=use_template)

    ideal_new_mids = ideal_mids + rotran[1]
    #av_chain_mids = sum(chain_mids) / 2.
    av_ideal_mids = sum(ideal_mids) / 2.
    chain_new_mids = np.dot(ideal_mids, rotran[0]) + rotran[1]
    #cud.pv('chain_mids')
    #cud.pv('chain_new_mids')

    return (bpdb.Vector(chain_new_mids[0]), bpdb.Vector(chain_new_mids[1]))


def get_twists_core(chain, start1, start2, end1, end2, mids=None, method=cbc.Configuration.mids_method):
    '''
    Get the vectors indicating the twist of the cylinder. In actuality, this will
    be the projection of the (ca_start1 - mid1) onto the plane defined by (mid2 - mid1).
    '''

    #mids = get_mids_core(chain, start1, start2, end1, end2)
    if mids == None:
        mids = get_mids(chain, [start1, end1, end2, start2], method = cbc.Configuration.mids_method)

    '''
    start_vec1 = chain[start1][catom_name].get_vector() - chain[start2][catom_name].get_vector()
    end_vec1 = chain[end1][catom_name].get_vector() - chain[end2][catom_name].get_vector()

    '''
    start_vec1 = chain[start1][catom_name].get_vector() - mids[0]
    end_vec1 = chain[end1][catom_name].get_vector() - mids[1]

    start_vec1a = chain[start2][catom_name].get_vector() - mids[0]
    end_vec1a = chain[end2][catom_name].get_vector() - mids[1]

    notch1 = cuv.vector_rejection(start_vec1.get_array(), (mids[0] - mids[1]).get_array())
    notch2 = cuv.vector_rejection(end_vec1.get_array(), (mids[1] - mids[0]).get_array())

    notch1a = cuv.vector_rejection(start_vec1a.get_array(), (mids[0] - mids[1]).get_array())
    notch2a = cuv.vector_rejection(end_vec1a.get_array(), (mids[1] - mids[0]).get_array())

    #print >>sys.stderr, "twist_angle_1:", cuv.vec_angle(notch1, notch1a)
    #print >>sys.stderr, "twist_angle_2:", cuv.vec_angle(notch2, notch2a)

    return (cuv.normalize(notch1 + notch1a), cuv.normalize(notch2 + notch2a))
    #return (normalize(notch1), normalize(notch2))


def get_mids(chain, define, method = cbc.Configuration.mids_method):
    '''
    Get the mid points of the abstract cylinder which represents a helix.

    @param chain: The Bio.PDB representation of the 3D structure.
    @param define: The define of the helix, as per the BulgeGraph definition standard.
    @return: An array of two vectors representing the two endpoints of the helix.
    '''

    if method == 'template':
        return get_mids_core(chain, int(define[0]), int(define[3]), int(define[1]), int(define[2]))
    elif method == 'fit':
        return get_mids_fit_method(chain, int(define[0]), int(define[3]), int(define[1]), int(define[2]))
    elif method == 'superimpose':
        return get_mids_core(chain, int(define[0]), int(define[3]), int(define[1]), int(define[2]), use_template=False) 
    elif method == 'estimate':
        return estimate_mids_core(chain, int(define[0]), int(define[3]), int(define[1]), int(define[2])) 

def get_twists(chain, define, mids=None, method=cbc.Configuration.mids_method):
    '''
    Get the projection of the (ca - mids) vectors onto the helix axis. This, in a sense
    will define how much the helix twists.

    @param chain: The Bio.PDB representation of the 3D structure.
    @param define: The define of the helix, as per the BulgeGraph definition standard.
    @return: Two vectors which represent the twist of the helix.
    '''

    return get_twists_core(chain, int(define[0]), int(define[3]), int(define[1]), int(define[2]), mids, method)

def get_helix_vector(chain, start1, start2, end1, end2):
    (mid1, mid2) = get_mids(chain, start1, start2, end1, end2)
    return mid2 - mid1

def virtual_res_3d_pos_core(coords, twists, i, stem_len, stem_inv = None):
    '''
    Calculate the virtual position of the i'th nucleotide in the stem.

    The virtual position extrapolates the position of the residues based
    on the twists of the helix.

    @param bg: The BulgeGraph structure
    @param stem: The name of the stem
    @param i: The i'th residue of the stem

    @return: A tuple containing the point located on the axis of the stem
             and a vector away from that point in the direction of the 
             residue.
    '''
    #stem_len = bg.defines[stem][1] - bg.defines[stem][0] + 1
    stem_vec = coords[1] - coords[0]

    # the position of the virtual residue along the axis of
    # the stem
    vres_stem_pos = coords[0] + (i / float(stem_len-1)) * stem_vec

    # the angle of the second twist with respect to the first
    if stem_inv == None:
        stem_basis = cuv.create_orthonormal_basis(stem_vec, twists[0])
        t2 = cuv.change_basis(twists[1], stem_basis, cuv.standard_basis)
    else:
        t2 = np.dot(stem_inv, twists[1])

    ang = cum.atan3(t2[2], t2[1])

    # calculated from an ideal length 30 helix
    average_ang_per_nt = 0.636738030735
    expected_ang = (stem_len-1) * average_ang_per_nt
    expected_dev = expected_ang
    while (expected_dev - (2 * m.pi) > 0):
        expected_dev -= 2 * m.pi


    dev = min(abs(ang - expected_dev), abs(2 * m.pi + ang - expected_dev), 
              abs(2 * m.pi + expected_dev - ang))

    if ang < expected_dev:
        forward = 2 * m.pi + ang - expected_dev
        backward = expected_dev - ang
    else:
        forward = ang - expected_dev
        backward = 2 * m.pi + expected_dev - ang

    if forward < backward:
        ang = expected_ang + forward
    else:
        ang = expected_ang - backward

    ang_per_nt = ang / float(stem_len-1)
    #cud.pv('expected_ang, ang')
    ang = ang_per_nt * i

    # the basis vectors for the helix along which the
    # virtual residues will residue
    u = twists[0]
    v = cuv.normalize(np.cross(stem_vec, twists[0]))
    
    ang_offset = 0.9
    # equation for a circle in 3-space
    return (vres_stem_pos, 
            u * m.cos(ang) + v * m.sin(ang),
            u * m.cos(ang + ang_offset) + v * m.sin(ang + ang_offset),
            u * m.cos(ang - ang_offset) + v * m.sin(ang - ang_offset))

def virtual_res_3d_pos(bg, stem, i, stem_inv = None):
    return virtual_res_3d_pos_core(bg.coords[stem], bg.twists[stem], i, bg.stem_length(stem), stem_inv)

def bg_virtual_residues(bg):
    vress = []

    for s in bg.stems():
        for i in range(bg.stem_length(s)):
            vres = virtual_res_3d_pos(bg, s, i)
            vress += [vres[0] + vres[2], vres[0] + vres[3]]

    return np.array(vress)

def virtual_res_basis_core(coords, twists, i, stem_len, vec = None):
    '''
    Define a basis based on the location of a virtual stem residue.

    The basis will be defined by the direction of the stem, the direction
    of the virtual residue.

    @param bg: The BulgeGraph structure
    @param stem: The name of the stem
    @param i: The i'th residue of the stem

    @return: A 3x3 matrix defining the coordinate system above.
    '''

    if vec == None:
        (pos, vec, vec_l, vec_r) = virtual_res_3d_pos_core(coords, twists, i, stem_len)

    stem_vec = coords[1] - coords[0]

    return cuv.create_orthonormal_basis(stem_vec, vec)

def virtual_res_basis(bg, stem, i, vec = None):
    return virtual_res_basis_core(bg.coords[stem], bg.twists[stem], i, bg.stem_length(stem), vec)


def pos_to_spos(bg, s1, i1, s2, i2):
    '''
    Convert the location of s2, i2 into the coordinate system
    defined by (s1, i1)

    @param bg: The BulgeGraph containing the stems
    @param s1: The basis stem name
    @param i1: The basis res position
    @param s2: The stem containing the nucleotide to be converted
    @param i2: The nucleotide to be converted position
    '''
    sbasis = virtual_res_basis(bg, s1, i1)
    (s1_pos, s1_vec, s1_vec_l, s1_vec_r) = virtual_res_3d_pos(bg, s1, i1)
    (s2_pos, s2_vec, s2_vec_l, s2_vec_r) = virtual_res_3d_pos(bg, s2, i2)

    #rpos = (s2_pos + 7. * s2_vec) - (s1_pos + 7 * s1_vec)
    rpos = (s2_pos + 7. * s2_vec) - (s1_pos)
    #print "sbasis:", sbasis

    spos = cuv.change_basis(rpos, sbasis, cuv.standard_basis)

    '''
    if spos[1] ** 2 + spos[2] ** 2 < 5 and spos[0] > -5 and spos[0] < 5:
        print >>sys.stderr, "spos:", spos, s1, i1, s2, i2
    '''
    return spos

def spos_to_pos(bg, stem, i, spos):
    '''
    Convert the location of spos from the coordinate system
    of (stem, i) into the standard coordinate system.

    @param bg: The BulgeGraph
    @param stem: The name of the stem in the BulgeGraph
    @param i: The i'th residue in 'stem' which will define the coordinate system
    @param spos: The position in the alternate coordinate system

    @return: The coordinates in the cartesian coordinate system of the
        rest of the model.
    '''
    sbasis = virtual_res_basis(bg, stem, i)
    (s1_pos, s1_vec, s1_vec_l, s1_vec_r) = virtual_res_3d_pos(bg, stem, i)
    pos = cuv.change_basis(spos, cuv.standard_basis, sbasis)
    return pos + (s1_pos + s1_vec)

def get_residue_type(i, stem_len):
    '''
    Each nucleotide will be classified according to its position
    within the stem. That way, the distribution of surrounding
    nucleotides will be conditioned on the type of nucleotides.

    This is important due to the fact that nucleotides at the end
    of a stem may have other stem nucleotides in the direction
    of the stem vector. Nucleotides, in the middle shoubulge not due
    to the excluded volume of the stem they occupy.

    @param i: The position of the nucleotide.
    @param stem_len: The length of the stem.

    @return: The type of nucleotide position.
    '''
    assert(i < stem_len)

    return 0

def junction_virtual_res_distance(bg, bulge):
    '''
    Compute the distance between the two virtual residues flanking
    a bulge region.

    @param bg: The BulgeGraph containing the bulge.
    @param bulge: The name of the bulge.
    '''
    connecting_stems = list(bg.edges[bulge])

    (s1b, s1e) = bg.get_sides(connecting_stems[0], bulge)
    (s2b, s2e) = bg.get_sides(connecting_stems[1], bulge)

    if s1b == 1:
        (vr1_p, vr1_v, vr1_v_l, vr1_v_r) = bg.v3dposs[connecting_stems[0]][bg.stem_length(connecting_stems[0]) - 1]
    else:
        (vr1_p, vr1_v, vr1_v_l, vr1_v_r) = bg.v3dposs[connecting_stems[0]][0]

    if s2b == 1:
        (vr2_p, vr2_v, vr2_v_l, vr2_v_r) = bg.v3dposs[connecting_stems[1]][bg.stem_length(connecting_stems[1]) - 1]
    else:
        (vr2_p, vr2_v, vr2_v_l, vr2_v_r) = bg.v3dposs[connecting_stems[1]][0]

    dist2 = cuv.vec_distance((vr1_p + 7 * vr1_v), (vr2_p + 7. * vr2_v))
    return dist2

def get_strand_atom_vrn(bg, s, i):
    '''
    Return the strand and which atom to use for the adjacent
    nucleotide distance calculation.
    '''
    if i == 0:
        return (0, 'P', 0)

    # this might have to just be bg.stem_length(s)
    if i == 1:
        return (0, 'O3*', bg.stem_length(s)-1)
    if i == 2:
        return (1, 'P', bg.stem_length(s) - 1)
    if i == 3:
        return (1, 'O3*', 0)

def junction_virtual_atom_distance(bg, bulge):
    '''
    Compute the distance between the O3' atom and P' atom
    of the two residues that flank the junction segment.

    @param bg: The BulgeGraph containing the bulge.
    @param bulge: The name of the bulge

    @return: A single number corresponding to the distance above.
    '''
    connecting_stems = bg.connections(bulge)

    (i1, k1) = bg.get_sides_plus(connecting_stems[0], bulge)
    (i2, k2) = bg.get_sides_plus(connecting_stems[1], bulge)

    r1 = bg.seq[bg.defines[connecting_stems[0]][i1] - 1]
    r2 = bg.seq[bg.defines[connecting_stems[0]][i2] - 1]

    (strand1, a1, vrn1) = get_strand_atom_vrn(bg, connecting_stems[0], i1)
    (strand2, a2, vrn2) = get_strand_atom_vrn(bg, connecting_stems[1], i2)

    a1_pos = cua.avg_stem_vres_atom_coords[strand1][r1][a1]
    a2_pos = cua.avg_stem_vres_atom_coords[strand2][r2][a2]

    vpos1 = bg.vposs[connecting_stems[0]][vrn1]
    vbasis1 = bg.vbases[connecting_stems[0]][vrn1].transpose()

    vpos2 = bg.vposs[connecting_stems[1]][vrn2]
    vbasis2 = bg.vbases[connecting_stems[1]][vrn2].transpose()

    a1_npos = np.dot(vbasis1, a1_pos) + vpos1
    a2_npos = np.dot(vbasis2, a2_pos) + vpos2

    return cuv.magnitude(a1_npos - a2_npos)

def add_virtual_residues(bg, stem):
    '''
    Create all of the virtual residues and the associated 
    bases and inverses.

    @param bg: The BulgeGraph containing the stem
    @param stem: The name of the stem to be included
    '''
    stem_vec = bg.coords[stem][1] - bg.coords[stem][0]
    stem_basis = cuv.create_orthonormal_basis(stem_vec, bg.get_twists(stem)[0])
    stem_inv = nl.inv(stem_basis.transpose())

    bg.bases[stem] = stem_basis
    bg.stem_invs[stem] = stem_inv

    for i in range(bg.stem_length(stem)):
        vpos = virtual_res_3d_pos(bg, stem, i, stem_inv = stem_inv)
        vbasis = virtual_res_basis(bg, stem, i, vec=vpos[1])
        vinv = nl.inv(vbasis.transpose())

        bg.vposs[stem][i] = vpos[0]
        bg.vvecs[stem][i] = vpos[1]
        bg.v3dposs[stem][i] = vpos
        bg.vbases[stem][i] = vbasis
        bg.vinvs[stem][i] = vinv

def stem_vres_reference_atoms(bg, chain, s, i):
    '''
    Calculate the position of each atom in the reference of the
    stem and virtual residue.

    @param bg: The BulgeGraph
    @param chain: The PDB representation of the chain
    @param s: The stem identifier
    @param i: The i'th base-pair in the stem

    @return (origin, bases, [dict(atoms), dict(atoms)]) 
        The origin of the coordinate system (vpos)
        The basises (one for each nucleotide)
        Two dictionaries containing the positions of each atom in its respective
            coordinate system.
    '''
    coords = [dict(), dict()]
    (vpos, vvec, vvec_l, vvec_r) = virtual_res_3d_pos(bg, s, i)
    vec1 = cuv.normalize(bg.coords[s][1] - bg.coords[s][0])
    vec2 = cuv.normalize(vvec)
    basis = cuv.create_orthonormal_basis(vec1, vec2)

    for k in range(2):


        #r = bg.defines[s][0 + 2*k] + i
        if k == 0:
            r = bg.defines[s][0] + i
        else:
            r = bg.defines[s][3] - i

        for atom in cup.all_rna_atoms:
            try:
                c = chain[r][atom].coord

                new_c =  cuv.change_basis(c - vpos, basis, cuv.standard_basis)
                coords[k][atom] = new_c

            except KeyError:
                continue

    return (vpos, basis, coords)

def bounding_boxes(bg, chain, s, i):
    '''
    Return the bounding boxes of the two nucleotides at the
    i'th position on the stem.

    @param bg: The BulgeGraph
    @param chain: The PDB representation of the chain
    @param s: The stem identifier
    @param i: The i'th base-pair in the stem

    @return (origin, bases, [(c1, c2), (c1, c2)]) The bases (one for each nucleotide) and the corners defining the bounding box
            of the two nucleotides
    '''

    (vpos, bases, atoms) = stem_vres_reference_atoms(bg, chain, s, i)
    corners = []

    for k in range(2):
        basis = bases[k]

        min_c = [10000., 10000., 10000.]
        max_c = [-10000., -10000., -10000.]

        for atom in atoms[k].values():
            for j in range(3):
                min_c[j] = min(min_c[j], atom[j])
                max_c[j] = max(max_c[j], atom[j])
        n = min_c
        x = max_c
        corners += [(n, x)]
    return (vpos, bases, corners)

def virtual_residue_atoms(bg, s, i, strand=0, basis=None, vpos=None, vvec=None):
    '''
    Return two sets of atoms for the virtual residue. One for the nucleotide
    on each strand.

    @param bg: The BulgeGraph
    @param s: The stem
    @param i: The virtual residue number
    @param strand: The strand for which to get the virtual atoms
    '''
    '''
    if vpos == None or vvec == None:
        (vpos, vvec, vvec_l, vvec_r) = virtual_res_3d_pos(bg, s, i)
    if basis == None:
        basis = virtual_res_basis(bg, s, i, vvec).transpose()
    '''

    vpos = bg.vposs[s][i]
    basis = bg.vbases[s][i].transpose()

    rs = (bg.seq[bg.defines[s][0] + i - 1], bg.seq[bg.defines[s][3] - i -1 ])

    new_atoms = dict()

    for a in cua.avg_stem_vres_atom_coords[strand][rs[strand]].items():
        coords = a[1]
        #new_coords = cuv.change_basis(coords, cuv.standard_basis, basis) + vpos

        new_coords = np.dot(basis, coords) + vpos
        #new_coords2 = cuv.change_basis(coords, cuv.standard_basis, basis)

        new_atoms[a[0]] = new_coords

    return new_atoms

def calc_R(xc, yc, p):
    """ calculate the distance of each 2D points from the center (xc, yc) """
    return np.sqrt((p[:,0] - xc) ** 2 + (p[:,1] - yc)**2)

def f_2(c, p):
    """ calculate the algebraic distance between the data points and the mean circle centered at c=(xc, yc) """
    Ri = calc_R(*c, p=p)
    return Ri - Ri.mean()

def circle_fit(p):
    x = p[:,0]
    y = p[:,1]
    x_m = np.mean(x)
    y_m = np.mean(y)
    
    u = x - x_m
    v = y - y_m

        # linear system defining the center (uc, vc) in reduced coordinates:
    #    Suu * uc +  Suv * vc = (Suuu + Suvv)/2
    #    Suv * uc +  Svv * vc = (Suuv + Svvv)/2
    Suv  = sum(u*v)
    Suu  = sum(u**2)
    Svv  = sum(v**2)
    Suuv = sum(u**2 * v)
    Suvv = sum(u * v**2)
    Suuu = sum(u**3)
    Svvv = sum(v**3)

    # Solving the linear system
    A = np.array([ [ Suu, Suv ], [Suv, Svv]])
    B = np.array([ Suuu + Suvv, Svvv + Suuv ])/2.0
    uc, vc = nl.solve(A, B)

    xc_1 = x_m + uc
    yc_1 = y_m + vc

    return (xc_1, yc_1)
    '''
    Ri_1     = sqrt((x-xc_1)**2 + (y-yc_1)**2)
    R_1      = mean(Ri_1)
    residu_1 = sum((Ri_1-R_1)**2)

    return (xc_1, yc_1, R_1)
    '''

def circle_error(c, p):
    errors = f_2(c,p)
    return sum([e ** 2 for e in errors])

def sum_square(nums):
    return sum([n ** 2 for n in nums])

def f_3(vec, points, est):
    """ calculate the optimal circle for the points (p) projected onto
    the plane orthogonal to v """
    basis = cuv.create_orthonormal_basis(vec)
    new_points = cuv.change_basis(points.T, basis, cuv.standard_basis).T
    p = new_points[:,1:]

    center_estimate = est
    #center_2, ier=so.leastsq(f_2, center_estimate,args=p)
    center_2 = circle_fit(p) 
    r = np.mean(calc_R(*center_2, p=p))

    #cud.pv('circle_error(center_2, p)')
    return f_2(center_2, p)

def fit_circle(mids, points, start_pos, end_pos):
    '''
    Calculate the projection of points on the plane normal to
    vec and fit a circle to them.
    '''
    vec = mids[1] - mids[0]
    basis = cuv.create_orthonormal_basis(vec)
    new_points = cuv.change_basis(points.T, basis, cuv.standard_basis).T
    p = new_points[:,1:]

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        v1, ier=so.leastsq(f_3, mids[1] - mids[0], args=(points, mids[0][1:])) #, maxfev=10000)

    basis1 = cuv.create_orthonormal_basis(v1)

    points1 = cuv.change_basis(points.T, basis1, cuv.standard_basis).T
    start_pos1 = cuv.change_basis(start_pos, basis1, cuv.standard_basis)
    end_pos1 = cuv.change_basis(end_pos, basis1, cuv.standard_basis)

    center_5 = circle_fit(points1[:,1:])

    mids_stem_basis = [[start_pos1[0], center_5[0], center_5[1]],
                       [end_pos1[0], center_5[0], center_5[1]]]
    mids_standard_basis = cuv.change_basis(np.array(mids_stem_basis).T,
                                           cuv.standard_basis, basis1).T
    '''
    # works!
    mids_stem_basis = [[nmids[0][0], center_4[0], center_4[1]],
                       [nmids[1][0], center_4[0], center_4[1]]]
    mids_standard_basis = cuv.change_basis(np.array(mids_stem_basis).T,
                                           cuv.standard_basis,
                                           basis).T
    '''
    return mids_standard_basis 

def extract_define_residues(define, chain):
    '''Extract the residues in the define and return them as a new chain.'''
    c = bpdb.Chain.Chain(' ')
    ranges = zip(*[iter(define)]*2)
    for r in ranges:
        for x in range(r[0], r[1]+1):
            c.add(chain[x])
    return c


def fit_circle_old(mids, points, start_pos, end_pos, chain, stem_length, define):
    '''
    Calculate the projection of points on the plane normal to
    vec and fit a circle to them.
    '''
    filename = 'ideal_1_%d_%d_%d.pdb' % (stem_length, stem_length+1, stem_length*2)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        #cud.pv('filename')
        ideal_chain = list(bpdb.PDBParser().get_structure('test', 
                op.join(cbc.Configuration.stem_fragment_dir, filename)).get_chains())[0]

    cud.pv('cup.pdb_rmsd(ideal_chain, chain, sidechains=False, superimpose=False)')
    cud.pv('cup.pdb_rmsd(ideal_chain, chain, sidechains=False, superimpose=True)')
    '''
    rotran = cup.pdb_rmsd(ideal_chain, chain, sidechains=False, 
            superimpose=True, apply_sup=False)[2]
    '''
    rotran = cup.pdb_rmsd(chain, ideal_chain, sidechains=False, 
            superimpose=True, apply_sup=False)[2]

    ideal_mids = get_mids_core(ideal_chain, 1, stem_length*2, stem_length, stem_length+1)
    chain_mids = get_mids_core(chain, define[0], define[3], define[1], define[2])

    ideal_mids = np.array([ideal_mids[0].get_array(), ideal_mids[1].get_array()])
    #chain_mids = np.array([chain_mids[0].get_array(), chain_mids[1].get_array()])

    ideal_new_mids = ideal_mids + rotran[1]
    #av_chain_mids = sum(chain_mids) / 2.
    av_ideal_mids = sum(ideal_mids) / 2.
    chain_new_mids = np.dot(ideal_mids, rotran[0]) + rotran[1]

    cud.pv('ideal_new_mids')
    cud.pv('ideal_mids')
    cud.pv('chain_mids')
    cud.pv('chain_new_mids')

    return chain_new_mids

    vec = mids[1] - mids[0]
    basis = cuv.create_orthonormal_basis(vec)
    new_points = cuv.change_basis(points.T, basis, cuv.standard_basis).T
    start_pos0 = cuv.change_basis(start_pos, basis, cuv.standard_basis)
    end_pos0 = cuv.change_basis(end_pos, basis, cuv.standard_basis)
    p = new_points[:,1:]
    f_3(mids[1] - mids[0], points, mids[0][1:]) 

    v1, ier=so.leastsq(f_3, mids[1] - mids[0], args=(points, mids[0][1:]), maxfev=10000)
    basis1 = cuv.create_orthonormal_basis(v1)

    points1 = cuv.change_basis(points.T, basis1, cuv.standard_basis).T
    start_pos1 = cuv.change_basis(start_pos, basis1, cuv.standard_basis)
    end_pos1 = cuv.change_basis(end_pos, basis1, cuv.standard_basis)

    center_5 = circle_fit(points1[:,1:])
    r5 = np.mean(calc_R(*center_5, p=points1[:,1:]))
    #cud.pv('v1')
    #cud.pv('vec')
    #cud.pv('cuv.vec_angle(v1, vec)')
    #cud.pv('sum_square(f_3(v1,  points))')
    ss1 = sum_square(f_3(v1, points, mids[0][1:]))
    center_estimate = mids[0][1:]
    center_2, ier=so.leastsq(f_2, center_estimate,args=(p))
    center_4 = circle_fit(p)
    r2 = np.mean(calc_R(*center_2, p=p))
    r4 = np.mean(calc_R(*center_4, p=p))
    #cud.pv('(center_2, r, ier)')
    #cud.pv('f_2(center_2, p=p)')
    #cud.pv('f_2(mids[0][1:], p=p)')
    #cud.pv('f_2(center_2, p=p)')
    nmids = cuv.change_basis(np.array(mids).T, basis, cuv.standard_basis).T

    ss2 = sum_square(f_2(center_2, p))
    ss3 = sum_square(f_2(nmids[0][1:], p))
    ss4 = sum_square(f_2(center_4, p))

    #cud.pv('sum_square(f_2(center_2, p))')
    #cud.pv('sum_square(f_2(mids[0][1:], p))')
    #cud.pv('(ss3, ss4, ss1)')
    #cud.pv('(r2, r4, r5)')

    center_x = center_4
    rx = r4
    import matplotlib.pyplot as plt
    import pylab as pl
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1, adjustable='box', aspect=1)
    ax.plot(new_points[:,1], new_points[:,2], "o")
    ax.plot(points1[:,1], points1[:,2], "go")
    ax.plot(center_x[0], center_x[1], 'ro')
    ax.plot(center_5[0], center_5[1], 'co')

    mids = cuv.change_basis(np.array(mids).T, basis, cuv.standard_basis).T
    ax.plot(mids[0][1],mids[0][2], 'yo')

    #ax.plot(mids[0][1], mids[0][2], 'go')
    circle1=plt.Circle(center_x, rx,color='r',alpha=0.5)
    circle2=plt.Circle(mids[0][1:],rx,color='y',alpha=0.5)
    circle3=plt.Circle(center_5, r5, color='g', alpha=0.5)

    fig.gca().add_artist(circle1)
    fig.gca().add_artist(circle2)
    fig.gca().add_artist(circle3)

    mids_stem_basis = [[start_pos1[0], center_5[0], center_5[1]],
                       [end_pos1[0], center_5[0], center_5[1]]]
    #cud.pv('basis1')
    mids_standard_basis = cuv.change_basis(np.array(mids_stem_basis).T,
                                           cuv.standard_basis,
                                           basis1).T
    '''
    # works!
    mids_stem_basis = [[nmids[0][0], center_4[0], center_4[1]],
                       [nmids[1][0], center_4[0], center_4[1]]]
    mids_standard_basis = cuv.change_basis(np.array(mids_stem_basis).T,
                                           cuv.standard_basis,
                                           basis).T
    '''

    #cud.pv('v1')
    #cud.pv('cuv.normalize(mids_standard_basis[1] - mids_standard_basis[0])')
    plt.show()
    #cud.pv('cuv.magnitude(mids_standard_basis[1] - mids_standard_basis[0])')
    return mids_standard_basis 

def get_mids_fit_method(chain, start1, start2, end1, end2):
    '''
    Estimate the endpoints of the cylinder axis by fitting it and using
    the rmsd of the best fit circle as the function to minimize.
    '''
    start_vec = estimate_mids_core(chain, start1, start2, end1, end2)
    
    atom_poss = []
    residue_numbers = [i for i in range(start1, end1+1)]
    residue_numbers += [i for i in range(end2, start2+1)]

    ideal_chain = chain

    for rn in residue_numbers:
        #atom_poss += [chain[rn]['C1*'].get_vector().get_array()]
        try:
            '''
            for atom in chain[rn].get_list():
                atom_poss += [atom.get_vector().get_array()]
            '''

            atom_poss += [ideal_chain[rn]['P'].get_vector().get_array()]
            atom_poss += [ideal_chain[rn]['O3*'].get_vector().get_array()]
            atom_poss += [ideal_chain[rn]['C3*'].get_vector().get_array()]
            atom_poss += [ideal_chain[rn]['C4*'].get_vector().get_array()]
            atom_poss += [ideal_chain[rn]['C5*'].get_vector().get_array()]
            atom_poss += [ideal_chain[rn]['O5*'].get_vector().get_array()]
            atom_poss += [ideal_chain[rn]['C1*'].get_vector().get_array()]
        except KeyError as ke:
            pass
            #cud.pv('ke')

    points = np.array(atom_poss)
    mids = estimate_mids_core(chain, start1, start2, end1, end2)
    mids = np.array([mids[0].get_array(), mids[1].get_array()])
    vec = mids[1] - mids[0]
    v1, ier = so.leastsq(f_3, vec, args=(points, mids[0][1:]), maxfev=10000) 

    start_pos = (chain[start1]['C1*'].get_vector().get_array() +
                 chain[start2]['C1*'].get_vector().get_array()) / 2.
    end_pos = (chain[end1]['C1*'].get_vector().get_array() +
               chain[end2]['C1*'].get_vector().get_array()) / 2.

    basis1 = cuv.create_orthonormal_basis(v1)
    points1 = cuv.change_basis(points.T, basis1, cuv.standard_basis).T
    start_pos1 = cuv.change_basis(start_pos.T, basis1, cuv.standard_basis).T
    end_pos1 = cuv.change_basis(end_pos.T, basis1, cuv.standard_basis).T

    center = circle_fit(points1[:,1:])
    mids_stem_basis = [[start_pos1[0], center[0], center[1]],
                       [end_pos1[0], center[0], center[1]]]
    mids_standard_basis = cuv.change_basis(np.array(mids_stem_basis).T,
                                           cuv.standard_basis,
                                           basis1).T
    return [bpdb.Vector(mids_standard_basis[0]),
            bpdb.Vector(mids_standard_basis[1])]


def stem_vec_from_circle_fit(bg, chain, stem_name='s0'):
    '''
    Attempt to find the stem direcion vector given a set of atom positions.

    This will be done by solving for the stem_vector, then using that
    to project the atoms onto a plane orthogonal to that vector. On that plane,
    a circle will be fit to the positions of the atoms. The stem vector that
    gives a circle with the least residuals will be considered the ideal
    stem vector.

    @return: stem_vector
    '''
    start_vec = nr.random(3)
    atom_poss = []
    #stem_name = 's0'
    for rn in bg.stem_res_numbers(stem_name):
        #atom_poss += [chain[rn]['C1*'].get_vector().get_array()]
        try:
            atom_poss += [chain[rn]['P'].get_vector().get_array()]
            atom_poss += [chain[rn]['O3*'].get_vector().get_array()]
            atom_poss += [chain[rn]['C3*'].get_vector().get_array()]
            atom_poss += [chain[rn]['C4*'].get_vector().get_array()]
            atom_poss += [chain[rn]['C5*'].get_vector().get_array()]
            atom_poss += [chain[rn]['O5*'].get_vector().get_array()]
            atom_poss += [chain[rn]['C1*'].get_vector().get_array()]
        except KeyError as ke:
           # cud.pv('ke')
           pass

    start_pos = (chain[bg.defines[stem_name][0]]['C1*'].get_vector().get_array() +
                 chain[bg.defines[stem_name][3]]['C1*'].get_vector().get_array()) / 2.
    end_pos = (chain[bg.defines[stem_name][1]]['C1*'].get_vector().get_array() +
               chain[bg.defines[stem_name][2]]['C1*'].get_vector().get_array()) / 2.

    mids = get_mids(chain, bg.defines[stem_name])
    # use the original calculation to provide an estimate for the
    # optimized stem position calculation
    stem_chain = extract_define_residues(bg.defines[stem_name], chain)
    mids = (mids[0].get_array(), mids[1].get_array())
    return fit_circle_old(mids, np.array(atom_poss), 
            start_pos, end_pos, stem_chain, bg.stem_length(stem_name), bg.defines[stem_name])

def receptor_angle(bg, l, s):
    (i1, i2) = cuv.line_segment_distance(bg.coords[l][0],
                                         bg.coords[l][1],
                                         bg.coords[s][0],
                                         bg.coords[s][1])

    stem_len = bg.stem_length(s)
    stem_vec = bg.coords[s][1] - bg.coords[s][0]

    res_num = (stem_len - 1.) *  cuv.magnitude(i2 - bg.coords[s][0]) / cuv.magnitude(bg.coords[s][1] - bg.coords[s][0])
    vres = virtual_res_3d_pos(bg, s, res_num)[1]
    if cuv.magnitude(i2 - i1) == 0.:
        return 0.

    incoming_angle = cuv.vector_rejection(i1 - i2, stem_vec)

    return cuv.vec_angle(vres, incoming_angle)
