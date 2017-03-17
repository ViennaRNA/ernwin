import forgi.threedee.utilities.graph_pdb as ftug
import forgi.threedee.utilities.vector as ftuv
import forgi.utilities.debug as fud

import logging
log=logging.getLogger(__name__)

try:
  profile  #The @profile decorator from line_profiler (kernprof)
except:
  def profile(x): 
    return x

@profile
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
    angle1 = ftuv.vec_angle(cg.coords.get_direction(l2),
                            cg.coords.get_direction(l1))
    #fud.pv('angle1')

    tw = cg.get_twists(l2)

    if l2[0] != 's':
        angle2 = ftuv.vec_angle((tw[0] + tw[1]) / 2.,
                               i2 - i1)
    else:
        stem_len = cg.stem_length(l2)

        # Where along the helix our A-residue points to the minor groove.
        # This can be between residues. We express it as floating point nucleotide coordinates.
        # So 0.0 means at the first basepair, while 1.5 means between the second and the third basepair.
        pos = ftuv.magnitude(i2 - cg.coords[l2][0]) / ftuv.magnitude(cg.coords.get_direction(l2)) * (stem_len - 1)

        # The vector pointing to the minor groove, even if we are not at a virtual residue (pos is a float value)
        vec = ftug.virtual_res_3d_pos_core(cg.coords[l2], cg.twists[l2], pos, stem_len)[1]
        angle2 = ftuv.vec_angle(vec,
                               i2 - i1)

    dist = ftuv.vec_distance(i2, i1)

    return (dist, angle1, angle2)
