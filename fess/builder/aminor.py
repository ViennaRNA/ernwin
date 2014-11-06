import forgi.threedee.utilities.graph_pdb as ftug
import forgi.threedee.utilities.vector as ftuv
import forgi.utilities.debug as fud

def get_relative_orientation(cg, l1, l2):
    '''
    Return how l1 is related to l2 in terms of three parameters. l2 should
    be the receptor of a potential A-Minor interaction, whereas l1 should
    be the donor.

        1. Distance between the closest points of the two elements
        2. The angle between l2 and the vector between the two
        3. The angle between the minor groove of l2 and the vector between
           l1 and l2
    '''
    (i1, i2) = ftuv.line_segment_distance(cg.coords[l1][0],
                                          cg.coords[l1][1],
                                          cg.coords[l2][0],
                                          cg.coords[l2][1])

    '''
    angle1 = ftuv.vec_angle(cg.coords[l2][1] - cg.coords[l2][0],
                           i2 - i1)
    '''
    angle1 = ftuv.vec_angle(cg.coords[l2][1] - cg.coords[l2][0],
                            cg.coords[l1][1] - cg.coords[l1][0])
    #fud.pv('angle1')

    tw = cg.get_twists(l2)

    if l2[0] != 's':
        angle2 = ftuv.vec_angle((tw[0] + tw[1]) / 2.,
                               i2 - i1)
    else:
        stem_len = cg.stem_length(l2)

        pos = ftuv.magnitude(i2 - cg.coords[l2][0]) / ftuv.magnitude(cg.coords[l2][1] - cg.coords[l2][0]) * stem_len
        vec = ftug.virtual_res_3d_pos_core(cg.coords[l2], cg.twists[l2], pos, stem_len)[1]
        angle2 = ftuv.vec_angle(vec,
                               i2 - i1)

    dist = ftug.element_distance(cg, l1, l2)

    return (dist, angle1, angle2)
