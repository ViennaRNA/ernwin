from __future__ import absolute_import, division, print_function, unicode_literals
from builtins import (ascii, bytes, chr, dict, filter, hex, input, #pip install future
                      int, map, next, oct, open, pow, range, round,
                      str, super, zip)
from future.builtins.disabled import (apply, cmp, coerce, execfile,
                             file, long, raw_input, reduce, reload,
                             unicode, xrange, StandardError)
__metaclass__=object         



import forgi.threedee.utilities.graph_pdb as ftug
import forgi.threedee.utilities.vector as ftuv
import forgi.utilities.debug as fud
import forgi.graph.bulge_graph as fgb
from collections import namedtuple
import warnings
import logging
import numpy as np
import scipy.stats
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

    return dist, angle1, angle2


    
_AMGeometry = namedtuple("AMGeometry", ["pdb_id", "loop_name", "stem_name", "dist", "angle1", "angle2", "loop_sequence"])
class AMGeometry(_AMGeometry):
    def _asdict(self):
        d = super(AMGeometry, self)._asdict()
        d["loop_type"] = self.loop_name[0]
        return d
    @property
    def loop_type(self):
        return self.loop_name[0]
    
def parse_fred(cutoff_dist, all_cgs, fr3d_out):
    """
    Used by generate_target_distribution.

    :param all_cgs: A dictionary PDB_id (4-letter): [CoarseGrainRNA, ...]
    :param fr3d_out: A file opened for reading.
    """
    geometries = set()
    for line in fr3d_out:
        log.debug("Line '%s'.read", line.strip())
        parts = line.split()
        if len(parts) < 10: 
            continue
        if parts[0] == "Filename": 
            continue
        pdb_id = parts[0]
        if parts[0] in all_cgs:
            cgs = all_cgs[parts[0]]
            if len(cgs)==1:
                cg, = cgs
            else:
                a_res = fgb.resid_from_str("{chain}:{res}".format(chain = parts[8][0], res = parts[3]))
                for cg in cgs:
                    if a_res in cg.seq_ids:
                        break
                    else:
                        pass
                else:
                    log.warning("No CoarseGrainRNA found among those with PDB-ID {} that contains resid {}".format(parts[0], a_res))
                    warnings.warn("No CoarseGrainRNA found among those with PDB-ID {} that contains resid {}".format(parts[0], a_res))
                    continue
        else:
            log.warning("No CoarseGrainRNA found for FR3D annotation with PDB-ID {}. Possible PDB-IDs are {}".format(parts[0], all_cgs.keys()))
            warnings.warn("No CoarseGrainRNA found for FR3D annotation with PDB-ID {}. Possible PDB-IDs are {}".format(parts[0], all_cgs.keys()))
            continue
        try:
            nums = []
            for i in range(3):
                nums.append(fgb.resid_from_str("{chain}:{res}".format(chain = parts[8][i], res = parts[3+2*i])))
            nodes = list(map(lambda x: cg.get_node_from_residue_num(cg.seq_id_to_pos(x)), nums))
        except Exception as e:
            log.exception(e)
            continue
        if nodes[1]!=nodes[2] or nodes[1][0]!="s":
            log.warning("Parse_fred: No canonical stem: {} != {}".format(nodes[1], nodes[2]))
            continue #only look at canonical stems.
        if nodes[0][0]=="s":
            log.warning("Parse_fred: Stem-Stem A-Minor not allowed: {} -> {}".format(nodes[0], nodes[1]))
            continue #No A-Minor between two stems.
        if nodes[0] in cg.edges[nodes[1]]:
            log.warning("Parse_fred:  A-Minor between adjacent elements not allowed: {} -> {}".format(nodes[0], nodes[1]))
            continue #Only look at not connected elements
        dist, angle1, angle2 = get_relative_orientation(cg, nodes[0], nodes[1])
        if dist<=cutoff_dist and "A" in "".join(cg.get_define_seq_str(nodes[0])):
            geometries.add(AMGeometry(pdb_id, nodes[0], nodes[1], dist, angle1, angle2, "&".join(cg.get_define_seq_str(nodes[0]))))
    return geometries



def aminor_probability_function(aminor_geometries, all_geometries, loop_type):
    """
    :param aminor_geometries: A list or iterator of AMGeometry instances which correspond to true interactions.
    :param all_geometries: A list or iterator of AMGeometry instances which correspond to all combinations 
                           of stems and loops, independent if it is an A-Minor-Interaction or not.
    :param loop_type: The type of the loop which donates the Adenine. A single letter string (e.g. "h", "i", ...)
    :returns: A probability function which takes a triple (dist, angle1, angle2) and 
              returns the probability of this being an A-Minor interaction.
    """
    aminor_geometries = np.array([[ x.dist, x.angle1, x.angle2 ] 
                                  for x in aminor_geometries 
                                  if  x.loop_type == loop_type])
    all_geometries = np.array([[ x.dist, x.angle1, x.angle2 ] 
                                  for x in all_geometries 
                                  if  x.loop_type == loop_type])
    log.info("%d interactions and %d geometries given.", len(aminor_geometries), len(all_geometries))
    if len(aminor_geometries) == 0: #E.g. for f0
        return lambda x: np.array([0])
    # Overall Probability/ Frequency of A-Minor interactions
    log.info("len(aminor_geometries)=%d, len(all_geometries)=%d", len(aminor_geometries), len(all_geometries))
    p_interaction = len(aminor_geometries)/ len(all_geometries)
    log.info("p_interaction = %s", p_interaction)
    p_geometry_given_interaction = scipy.stats.gaussian_kde(aminor_geometries.T)
    p_geometry_all = scipy.stats.gaussian_kde(all_geometries.T)
    # According to Peter's Thesis:
    #Gives always 0
    p_interaction_given_geometry = lambda point: (p_geometry_given_interaction(point) / p_geometry_all(point)) * p_interaction 
    
    def p_function(point):
        numerator = p_geometry_given_interaction(point)
        #log.info("Numerator: %s", numerator)
        denom = p_geometry_all(point)
        #log.info("Denominator: %s", denom)
        #log.info("p_interaction %s", p_interaction)
        return numerator/denom*p_interaction

    # The version below was used by Peter to avoid too small denominators in case of pseudoknots.
    # Can give probabilities >1
    p_interaction_given_geometry = lambda point: (p_geometry_given_interaction(point)) / p_geometry_all(point) + p_geometry_given_interaction(point)
    #return p_interaction_given_geometry
    return p_function

def max_prob(loop, cg, prob_fun, cutoff_dist, domain = None):
    """
    Return the maximal probability for the loop form any A-Minor interaction.
    
    .. note::
    
        This is equivalent to what was used in ernwin 0.1. It is now replaced by total_prob.
    
    
    This function tries for interaction of the loop with
    all stems (except adjacent ones) and returns the maximum of 
    the probabilities.
    
    :param loop: The loop name, e.g. "h0"
                 
                 .. warning::
                 
                    The loop type (hairpin/interior) should correspond to the loop-type used 
                    for generating the probability function `prob_fun`
                    
    :param cg:          The CoarseGrainRNA.
    :param prob_fun:    A probability function. A function that takes a triple (distance, angle1, angle2) 
                        as returned by `get_relative_orientation` and returns a probability for 
                        this goemetry to correspond to an A-Minor interaction.
    :param cutoff_dist: Do not consider interactions between elements more 
                        than this many angstroms away.
    :param domain:      A list of element names. Only take these elements into account.
                        None to use the whole graph from cg.
    """
    # Code moved here from fbe.AMinorEnergy.eval_prob
    if domain is not None:
        stems = (s for s in domain if s[0]=="s")
    else:
        stems = cg.stem_iterator()
    probs = []
    for s in stems:
        if s in cg.edges[loop]:
            continue
        if not ftuv.elements_closer_than(cg.coords[loop][0],
                                 cg.coords[loop][1],
                                 cg.coords[s][0],
                                 cg.coords[s][1], 
                                 cutoff_dist):
            continue
            
        point = get_relative_orientation(cg, loop, s)
        p, = prob_fun(point)
        probs.append(p)
    if len(probs) == 0:
        log.debug("max_prob: Returning zero.")
        return 0
    max_prob = max(probs)
    log.debug("max_prob: Returning max(probs): %s", max_prob)
    return max_prob

def total_prob(loop, cg, prob_fun, cutoff_dist, domain = None):
    #return max_prob(loop, cg, prob_fun, cutoff_dist, domain)
    """
    Return the total probability for the loop form any A-Minor interaction.
    
    This function tries for interaction of the loop with
    all stems (except adjacent ones) and returns p1+(1-p1)*p2+... 
    where p1, p2, ... are the individual probabilities.
    
    :param loop: The loop name, e.g. "h0"
                 
                 .. warning::
                 
                    The loop type (hairpin/interior) should correspond to the loop-type used 
                    for generating the probability function `prob_fun`
                    
    :param cg:          The CoarseGrainRNA.
    :param prob_fun:    A probability function. A function that takes a triple (distance, angle1, angle2) 
                        as returned by `get_relative_orientation` and returns a probability for 
                        this goemetry to correspond to an A-Minor interaction.
    :param cutoff_dist: Do not consider interactions between elements more 
                        than this many angstroms away.
    :param domain:      A list of element names. Only take these elements into account.
                        None to use the whole graph from cg.
    """
    # Code moved here from fbe.AMinorEnergy.eval_prob
    if domain is not None:
        stems = (s for s in domain if s[0]=="s")
    else:
        stems = cg.stem_iterator()
    total_prob = 0
    for s in stems:
        if s in cg.edges[loop]:
            continue
        if not ftuv.elements_closer_than(cg.coords[loop][0],
                                 cg.coords[loop][1],
                                 cg.coords[s][0],
                                 cg.coords[s][1], 
                                 cutoff_dist):
            continue
            
        point = get_relative_orientation(cg, loop, s)
        p, = prob_fun(point)
        total_prob += (1-total_prob)*p
    log.debug("total_prob: Returning: %s", total_prob)
    return total_prob
