from __future__ import absolute_import, division, print_function, unicode_literals
from builtins import (ascii, bytes, chr, dict, filter, hex, input, #pip install future
                      map, next, oct, open, pow, range, round,
                      str, super, zip)
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

from logging_exceptions import log_to_exception

log=logging.getLogger(__name__)

try:
  profile  #The @profile decorator from line_profiler (kernprof)
except:
  def profile(x):
    return x

@profile
def get_relative_orientation(cg, loop, stem):
    '''
    Return how loop is related to stemin terms of three parameters.

    The stem is the receptor of a potential A-Minor interaction, whereas the
    loop is the donor.

    The 3 parameters are:

        1. Distance between the closest points of the two elements
        2. The angle between the stem and the vector between the two
        3. The angle between the minor groove of l2 and the vector between
           stem and loop
    '''

    #l1...loop, l2...stem
    (point_on_stem, point_on_loop) = ftuv.line_segment_distance(cg.coords[stem][0],
                                                                cg.coords[stem][1],
                                                                cg.coords[loop][0],
                                                                cg.coords[loop][1])
    conn_vec = point_on_loop-point_on_stem

    angle1 = ftuv.vec_angle(cg.coords.get_direction(stem),
                            conn_vec)
    # The direction of the stem vector is irrelevant, so
    # choose the smaller of the two angles between two lines
    if angle1>np.pi/2:
        angle1 = np.pi-angle1


    #fud.pv('angle1')

    tw = cg.get_twists(stem)

    dist = ftuv.magnitude(conn_vec)

    if dist==0:
        angle2=float("nan")
    else:
        if stem[0] != 's':
            raise ValueError("The receptor needs to be a stem, not {}".format(stem))
        else:
            stem_len = cg.stem_length(stem)

            # Where along the helix our A-residue points to the minor groove.
            # This can be between residues. We express it as floating point nucleotide coordinates.
            # So 0.0 means at the first basepair, while 1.5 means between the second and the third basepair.
            pos = ftuv.magnitude(point_on_stem - cg.coords[stem][0]) / ftuv.magnitude(cg.coords.get_direction(stem)) * (stem_len - 1)

            # The vector pointing to the minor groove, even if we are not at a virtual residue (pos is a float value)
            vec = ftug.virtual_res_3d_pos_core(cg.coords[stem], cg.twists[stem], pos, stem_len)[1]
            try:
                # Note: here the directions of both vectors are well defined,
                # so angles >90 degrees make sense.
                angle2 = ftuv.vec_angle(vec, conn_vec)
            except ValueError:
                if np.all(vec==0):
                    angle2=float("nan")
                else:
                    raise

    return dist, angle1, angle2



_AMGeometry = namedtuple("AMGeometry", ["pdb_id", "loop_name", "stem_name", "dist", "angle1", "angle2", "loop_sequence", "score", "annotation"])
class AMGeometry(_AMGeometry):
    def _asdict(self):
        d = super(AMGeometry, self)._asdict()
        d["loop_type"] = self.loop_name[0]
        return d
    @property
    def loop_type(self):
        return self.loop_name[0]

def _safe_resid_from_chain_res(chain, residue):
    try:
        return fgb.resid_from_str(str("{}:{}".format(chain,residue)))
    except ValueError as e:
        if residue.isdigit():
            with log_to_exception(log, e):
                log.error("Chain is '{}', res is '{}'".format(chain, residue))
            raise
        else:
            warnings.warn("Illegal residue number: '{}'.".format(residue))
            return

def _parse_fred_line(line, all_cgs, current_annotation):
    parts = line.split()
    if len(parts) < 10:
        return
    if parts[0] == "Filename":
        return
    pdb_id = parts[0]
    if parts[0] in all_cgs:
        cgs = all_cgs[parts[0]]
        if len(cgs)==1:
            cg, = cgs
        else:
            a_res = _safe_resid_from_chain_res(chain = parts[8][0], residue = parts[3])
            if a_res is None:
                return
            for cg in cgs:
                #Find the first matching cg.
                if a_res in cg.seq_ids:
                    break
            else:
                warnings.warn("No CoarseGrainRNA found among those with PDB-ID {} that contains resid {}".format(parts[0], a_res))
                return
    else:
        warnings.warn("No CoarseGrainRNA found for FR3D annotation with PDB-ID {}. Possible PDB-IDs are {}".format(parts[0], all_cgs.keys()))
        return

    # We now have the CoarseGrainRNA. Get the interacting cg-elements.
    nums = []
    for i in range(3):
        seq_id = _safe_resid_from_chain_res(chain = parts[8][i], residue = parts[3+2*i])
        if seq_id is None:
            return
        try:
            nums.append(cg.seq_id_to_pos(seq_id))
        except ValueError:
            if seq_id.chain not in cg.chain_ids:
                warnings.warn("Chain {!r} is not part of the cg {}.".format(seq_id.chain, cg.name))
                return
            else:
                raise
    nodes = list(map(lambda x: cg.get_node_from_residue_num(x), nums))
    if nodes[1]!=nodes[2] or nodes[1][0]!="s":
        warnings.warn("Parse_fred: No canonical stem: {} != {}.".format(nodes[1], nodes[2]))
        return #only look at canonical stems.
    if nodes[0][0]=="s":
        warnings.warn("Parse_fred: Stem-Stem A-Minor not allowed: {} -> {}.".format(nodes[0], nodes[1], line))
        return #No A-Minor between two stems.
    if nodes[0] in cg.edges[nodes[1]]:
        warnings.warn("Parse_fred:  A-Minor between adjacent elements not allowed: {} -> {}.".format(nodes[0], nodes[1]))
        return #Only look at not connected elements
    dist, angle1, angle2 = get_relative_orientation(cg, nodes[0], nodes[1])
    if np.isnan(angle1+angle2+dist):
        warnings.warn("Cannot get relative orientation. Zero-length element {}".format(nodes[0]))
        return
    return (AMGeometry(pdb_id, nodes[0], nodes[1], dist, angle1, angle2, "&".join(cg.get_define_seq_str(nodes[0])),float(parts[1]), current_annotation))

def parse_fred(cutoff_dist, all_cgs, fr3d_out):
    """
    Used by generate_target_distribution.

    :param all_cgs: A dictionary PDB_id (4-letter): [CoarseGrainRNA, ...]
    :param fr3d_out: A file opened for reading.

    :returns: A set of AMinor Geometries and a list of PDB-IDs, for which at
              least one line did not lead to an Aminor Geometry
    """
    problematic_pdbids = set()
    geometries = set()
    skipped = 0
    #: What type of AMinor interactions.
    #: Comments in the form "# AMinor 0" have to be added manually
    #: when copying the FR3D-output to a file.They should preced the
    #: actual FR3D-output to which they apply.
    current_annotation = "?"
    for line in fr3d_out:
        line=line.strip()
        if not line: #Empty line
            continue
        if line.startswith("# AMinor"):
            current_annotation = line[9:]
        elif line.startswith("#"):
            current_annotation = "?"
        log.debug("Line '%s'.read", line)
        geometry = _parse_fred_line(line, all_cgs, current_annotation)
        if geometry is None:
            skipped+=1
            if not (line.startswith("Filename") or line.startswith("#")):
                log.warning("Skipping line {!r}".format(line))
                problematic_pdbids.add(line.split()[0])
        elif geometry.dist>cutoff_dist:
            log.info("Skipping because of %f > %f (=cutoff dist): %r",
                     geometry.dist, cutoff_dist, line)
        elif "A" not in geometry.loop_sequence:
            warnings.warn("No adenine in loop %r for line %r", geometry.loop_name, line)
        else:
            geometries.add(geometry)
    return geometries, skipped #problematic_pdbids


def aminor_probability_function(aminor_geometries, non_aminor_geometries, loop_type):
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
    non_aminor_geometries = np.array([[ x.dist, x.angle1, x.angle2 ]
                                  for x in non_aminor_geometries
                                  if  x.loop_type == loop_type])
    log.info("%d interactions and %d non-interactions given for loop type %s.", len(aminor_geometries), len(non_aminor_geometries), loop_type)
    if len(aminor_geometries) < 3: #E.g. for f0
        return lambda x: np.array([0])
    # Overall Probability/ Frequency of A-Minor interactions
    log.info("len(aminor_geometries)=%d, len(non_aminor_geometries)=%d", len(aminor_geometries), len(non_aminor_geometries))
    p_interaction = len(aminor_geometries)/( len(non_aminor_geometries) + len(aminor_geometries))
    log.info("p_interaction = %s", p_interaction)
    p_geometry_given_interaction = scipy.stats.gaussian_kde(aminor_geometries.T)
    log.info("p_geometry_given_interaction done. Calculating p_geometry_all")
    p_geometry_non_interaction = scipy.stats.gaussian_kde(non_aminor_geometries.T)


    def p_function(point):
        numerator = p_geometry_given_interaction(point)*p_interaction
        denom = numerator+p_geometry_non_interaction(point)*(1-p_interaction)
        log.debug("P(I)= %f, P(nonI)= %f, P(geo|I)= %f, P(geo|nonI)=%f", p_interaction, 1-p_interaction, p_geometry_given_interaction(point), p_geometry_non_interaction(point))
        log.debug("prob: %f / %f = %f", numerator, denom, numerator/denom)
        return numerator/denom

    # The version below was used by Peter to avoid too small denominators in case of pseudoknots.
    # Can give probabilities >1
    #p_interaction_given_geometry = lambda point: (p_geometry_given_interaction(point)) / p_geometry_all(point) + p_geometry_given_interaction(point)
    #return p_interaction_given_geometry

    # New version
    return p_function


def total_prob(individual_probs):
    """
    Return the total probability for the loop to form at least one A-Minor interaction.

    This function tries for interaction of the loop with
    all stems (except adjacent ones) and returns p1+(1-p1)*p2+...
    where p1, p2, ... are the individual probabilities.
    """
    log.debug("Entering 'total_prob'")
    total_prob = 0
    for p in individual_probs:
        log.debug("Adding p = %f to total_prob = %f", p, total_prob)
        total_prob += (1-total_prob)*p
    log.debug("total_prob: Returning: %s", total_prob)
    if total_prob>1:
        log.error("Probability >1 for %s %s with domain %s", cg.name, loop, domain)
        assert False
    return total_prob

def iter_probs(loop, cg, prob_fun, cutoff_dist, domain = None):
    """
    Iterate over all stems and yield the probability for an interaction with loop.

    .. warning ::

        Do not rely on len(list(_iter_probs))!
        This function yields a single zero as a last element,
        which does not correspond to any stem.
        This is required to make sure at least one value is yielded in cases
        where no stem is close enough for interactions.

    For params: See the documentation of total_prob
    """
    log.debug("Entering '_iter_probs'")
    if domain is not None:
        stems = (s for s in domain if s[0]=="s")
    else:
        stems = cg.stem_iterator()
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
        if p>1:
            warnings.warn("Probability at %s is %f>1 for %s %s with domain %s" %(point, p, cg.name, loop, domain))
        yield(p)
    yield 0 #Always yield at least one number, so max(_iter_probs(...)) does not raise an error
