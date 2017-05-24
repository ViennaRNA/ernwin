#!/usr/bin/python
from __future__ import absolute_import, division, print_function, unicode_literals
from builtins import (ascii, bytes, chr, dict, filter, hex, input, #pip install future
                      int, map, next, oct, open, pow, range, round,
                      str, super, zip)
from future.builtins.disabled import (apply, cmp, coerce, execfile,
                             file, long, raw_input, reduce, reload,
                             unicode, xrange, StandardError)
__metaclass__=object

from .energy_abcs import EnergyFunction, CoarseGrainEnergy, DEFAULT_ENERGY_PREFACTOR
from collections import defaultdict
import random
import numpy as np
import Bio.KDTree as kd #KD-Trees for distance-calculations in point-cloud.
import warnings
import scipy.stats
import scipy.optimize
import scipy.ndimage
import scipy.misc
import forgi.threedee.utilities.vector as ftuv
import forgi.threedee.model.coarse_grain as ftmc
import forgi.threedee.utilities.graph_pdb as ftug
import fess.builder.aminor as fba
import fess.builder.models as cbm
import forgi.projection.hausdorff as fph
import forgi.projection.projection2d as fpp
import forgi.threedee.model.similarity as ftms
import forgi.threedee.model.descriptors as ftmd
from ..aux.utils import get_all_subclasses, get_version_string
import os.path as op
import os
import pandas as pd
import pkgutil as pu
import StringIO
import math
import re
import sys
import inspect
import itertools
from pprint import pprint
import logging
import time
log = logging.getLogger(__name__)


try:
  profile  #The @profile decorator from line_profiler (kernprof)
except:
  def profile(x):
    return x

INCR = 0.01


def load_local_data(filename):
    '''
    Load a data file that is located within this
    package.

    An example is something like 'stats/longrange.stats'.

    @param: A filename relative to the base directory of the package.
    @return: A generator iterating over the lines in the file.
    '''
    data = pu.get_data('fess', filename)

    return StringIO.StringIO(data)

class RandomEnergy(EnergyFunction):
    _shortname = "RND"
    HELPTEXT = "       {:3}:  Random Energy".format(_shortname)
    def eval_energy(self, cg, background=None, nodes=None, **kwargs):
        return self.prefactor * random.random() + self.adjustment

class ConstantEnergy(EnergyFunction):
    _shortname = "CNST"
    HELPTEXT = "       {:3}:  Constant Energy".format(_shortname)
    @classmethod
    def from_cg(cls, cg, prefactor, adjustment, **kwargs):
        if adjustment is not None:
            warnings.warn("Adjustment {} ignored for constant energy CNST".format(adjustment))
        return cls(prefactor)

    def __init__(self, value = 0):
        super(ConstantEnergy, self).__init__(prefactor = value)
    def eval_energy(self, cg, background=None, nodes=None, **kwargs):
        return self.prefactor + self.adjustment

class CheatingEnergy(EnergyFunction):
    _shortname = "CHE"
    HELPTEXT = "       {:3}:  Cheating Energy. Tries to minimize the RMSD.".format(_shortname)
    @classmethod
    def from_cg(cls, cg, prefactor, adjustment, **kwargs):
        return cls(cg, prefactor, adjustment)
    def __init__(self, ref_cg, prefactor = None, adjustment = None):
        """
        The energy is RMSD**adjustment*prefactor
        """
        super(CheatingEnergy, self).__init__(prefactor = prefactor, adjustment = adjustment)
        self.real_residues = ftug.bg_virtual_residues(ref_cg)
    def eval_energy(self, cg, background=None, nodes=None, **kwargs):
        new_residues = ftug.bg_virtual_residues(cg)
        return  ftms.rmsd(self.real_residues, new_residues)**self.adjustment*self.prefactor

class ProjectionMatchEnergy(EnergyFunction):
    _shortname = "PRO"
    HELPTEXT = ("       {:3}:  Match Projection distances. \n"
               "             Requires the projected distances.".format(_shortname))
    @classmethod
    def from_cg(cls, cg, prefactor, adjustment, pro_distances, **kwargs):
        """
        :param pro_distances: A dictionary tuple(str,str):float that maps pairs
                              of element names to distances in the projected plane (in Angstrom)
        """
        if adjustment is not None:
            warnings.warn("Adjustment {} ignored for PRO energy".format(adjustment))
        return cls(pro_distances, prefactor)

    def __init__(self, distances={}, prefactor=None):
        """
        :param directions: A dict where the keys are tuples of coarse grain element names
                           (e.g.: ("h1","m1")) and the values are the distance
                           IN THE PROJECTED PLANE (i.e. in the micrograph).
        """
        if prefactor is None:
            prefactor = 1.
        super(ProjectionMatchEnergy, self).__init__(prefactor = prefactor)
        self.distances=distances
        #: The optimal projection direction for the last accepted step.
        self.accepted_projDir=np.array([1.,1.,1.])
        self.projDir=np.array([1.,1.,1.])
        self.start_points=self._get_start_points(60)
    def _get_start_points(self, numPoints):
        """
        Return numPoints equally-distributed points on half the unit sphere.

        Implements the 2nd algorithm from https://www.cmu.edu/biolphys/deserno/pdf/sphere_equi.pdf
        """
        numPoints=2*numPoints
        a=4*math.pi/numPoints
        d=math.sqrt(a)
        Mt=int(math.pi/d)
        dt=math.pi/Mt
        df=a/dt
        points=[]
        for m in range(Mt):
            theta=math.pi*(m+0.5)/Mt
            Mf=int(2*math.pi*math.sin(theta)/df)
            for n in range(Mf):
                phi=2*math.pi*n/Mf
                points.append((math.sin(theta)*math.cos(theta),
                               math.sin(theta)*math.sin(phi),
                               math.cos(theta)))
        return [np.array(p) for p in points if p[0]>=0]

    def optimizeProjectionDistance(self, p):
        """
        :param c: A numpy array. The projection vector that HAS TO BE NORMALIZED

        Let theta be the angle between the p, normal vector of the plane of projection,
        and the vector a which is projected onto that plain.
        Then dist_2d=dist_3d*cos(alpha)=dist_3d*cos(90 degrees-theta)
        Then the angle theta is given by: cos(theta)=a*p/len(a) if p is normalized.
        Then dist_2d=dist_3d*sqrt(1-(a*p/len(a))^2)
        This means dist_2d^2=dist3d^2*(1-(a*p/len(a))^2)
        And dist_3d^2-dist2d^2=(a*p/len(a))^2 (Eq.1)
        We now search for the projection angle p which best generates the projected
        lengths given the current 3D structure.
        In the optimal case, Eq 1 holds for some values of p for all vectors a.
        In the suboptimal case, there is a squared deviation d:
        d=(a*p/len(a))^2-(dist_3d^2-dist2d^2)
        We want to minimize the sum over all square deviations for all distances
        considered in this energy function.
        minimize sum(d)=sum((a*p/len(a))^2-(dist_3d^2-dist2d^2)) with respect to p
        for all d, i.e. for all a, dist_3d, dist_2d considered.
        Under the side constraint that p has to be normalized.
        """
        x=0
        #Sum over all given distances
        for (s,e), dist in self.distances.items():
            sc=self.cg.coords[s]
            ec=self.cg.coords[e]
            #The middle point of the cg element
            start=(sc[0]+sc[1])/2
            end=(ec[0]+ec[1])/2
            a=end-start
            lengthDifferenceExperiment=ftuv.magnitude(a)**2-dist**2
            lengthDifferenceGivenP=(p[0]*a[0]+p[1]*a[1]+p[2]*a[2])**2
            # Add the deviation between (3d length-2d length)**2 observed vs calculated
            # for the given projection angle
            x+=abs(lengthDifferenceGivenP-lengthDifferenceExperiment)
        return x

    def eval_energy(self, cg, background=None, nodes=None, **kwargs):
        """
        Returns a measure that is zero for a perfect fit with the target
        projection distances and increases as the fit gets worse.

        This function tries to minimize its value over all projection directions.
        A global optimization is attempted and we are only interested in projection directions
        from the origin to points on half the unit sphere. Thus we first sample the energy for
        multiple, equally distributed directions. From the best one we operform a local optimization.

        For our specific purpose, where our starting points can give a good overview over
        the landscape, this is better and faster than more sophisticated
        global optimization techniques.
        """
        self.cg=cg
        # The projection vector has to be normalized
        c1={'type':'eq', 'fun':lambda x: x[0]**2+x[1]**2+x[2]**2-1}
        best_start=np.array([0.,0.,0.])
        best_score=float('inf')
        for direction in self.start_points:
            score=self.optimizeProjectionDistance(direction)
            if score<best_score:
                best_start=direction
                best_score=score
        opt=scipy.optimize.minimize(self.optimizeProjectionDistance, direction,
                                    constraints=c1, options={"maxiter":200} )
        if opt.success:
            self.projDir=opt.x
            score = math.sqrt(opt.fun)/len(self.distances)
            self._last_measure = score
            return self.prefactor*score
        else:
            self._last_measure = 10**11
            return 10**11

    def accept_last_measure(self):
        super(ProjectionMatchEnergy, self).accept_last_measure()
        self.accepted_projDir=self.projDir

class FPPEnergy(EnergyFunction):
    _shortname = "FPP"
    name = "Four-Point-Projection-Energy"
    HELPTEXT = ("       {:3}:  4 point projection energy.\n"
                "             Select 4 landmarks in the projected image.".format(_shortname))

    @classmethod
    def from_cg(cls, cg, prefactor, adjustment, fpp_landmarks, fpp_ref_image, fpp_scale, **kwargs):
        """
        :param fpp_ref_image: The path to the reference image
        :param fpp_landmarks: A list of 3-tuples (residue_num, x (in pixels), y (in pixels))
        :param fpp_scale: The width of the image in Angstrom
        """
        if adjustment is not None:
            warnings.warn("Adjustment {} ignored for FPP energy".format(adjustment))
        return cls(prefactor, fpp_landmarks, fpp_scale, fpp_ref_image)
    def __init__(self, pre, landmarks, scale, ref_image):
        super(FPPEnergy, self).__init__(prefactor = pre)
        self.landmarks = landmarks
        self.scale = scale
        self.ref_image = fpp.to_grayscale(scipy.ndimage.imread(ref_image))

    @profile
    def eval_energy(self, cg, background=True, nodes=None, **kwargs):
        steplength = self.scale/self.ref_image.shape[0]
        ### Step 1: Calculate projection direction:
        vectors3d, angles, penalties = self._generate_equations(cg)
        a = np.array(vectors3d[:3])
        b = np.array(angles[:3])
        #Solve the equation system ax = b
        try:
            projection_direction = np.linalg.solve(a, b)
        except Exception as e: #Degenerate equations. Estimate the solution from full set of equations
            log.info("{}: USING LSTSQ".format(e), exc_info=True)
            projection_direction = np.linalg.lstsq(np.array(vectors3d), np.array(angles))[0] #lstsq instead of solve, because system may be underdetermined

        ###
        #sm.bg.project_from = ftuv.normalize(-projection_direction)
        #proj =  fpp.Projection2D(sm.bg, project_virtual_atoms = "selected")
        #orig_img, _ = proj.rasterize(self.ref_image.shape[0], warn = False)
        proj, angleA, offset_centroidA, bs = self._find_offset(cg, projection_direction, mirror = True)
        img1, _ = proj.rasterize(self.ref_image.shape[0], bs, rotate = math.degrees(angleA), warn = False, virtual_residues = False)
        scoreA = fph.combined_distance(img1, self.ref_image)
        proj, angleB, offset_centroidB, bs = self._find_offset(cg, projection_direction, mirror = False)
        img2, _ = proj.rasterize(self.ref_image.shape[0], bs, rotate = math.degrees(angleB), warn = False, virtual_residues = False)
        scoreB = fph.combined_distance(img2, self.ref_image)
        if scoreA<scoreB:
            cg.project_from = ftuv.normalize(-projection_direction)
            score, img, params = fph.locally_minimal_distance(self.ref_image, self.scale, cg,
                                                          math.degrees(angleA), offset_centroidA,
                                                          None, distance=fph.combined_distance,
                                                          maxiter=200, virtual_atoms="selected")
        else:
           score, img, params = fph.locally_minimal_distance(self.ref_image, self.scale, cg,
                                                         math.degrees(angleB), offset_centroidB,
                                                         None, distance=fph.combined_distance,
                                                         maxiter=200, virtual_atoms="selected")
        #lmimg = np.zeros_like(img)
        #for l in self.landmarks:
        #    lmimg[l[2],l[1]]=1
        #import matplotlib.pyplot as plt
        #fig, ax  = plt.subplots(2,3)
        #ax[0,0].imshow(self.ref_image, interpolation="none")
        #ax[0,0].set_title("Reference")
        #ax[1,0].imshow(lmimg, interpolation="none")
        #ax[1,0].set_title("Landmarks:")
        #ax[0,1].imshow(img1, interpolation="none")
        #ax[0,1].set_title("IntermediateA {}".format(scoreA))
        #ax[1,1].imshow(img2, interpolation="none")
        #ax[1,1].set_title("IntermediateB {}".format(scoreB))
        #ax[0,2].imshow(img, interpolation="none")
        #ax[0,2].set_title("Final {}".format(score))
        #ax[1,2].imshow(orig_img, interpolation="none")
        #ax[1,2].set_title("Rot/ offset unknown")
        #plt.show()
        self._last_measure = score
        return score*self.prefactor

    @profile
    def _find_offset(self, cg, projection_direction, mirror = False):
        if mirror: projection_direction = -projection_direction
        steplength = self.scale/self.ref_image.shape[0]
        cg.project_from = ftuv.normalize(projection_direction)
        ### Step 2: Find out offset and rotation.
        proj = fpp.Projection2D(cg, project_virtual_residues = [ x[0] for x in self.landmarks], project_virtual_atoms = "selected")
        target = []
        current = []
        for l in self.landmarks:
            target.append(np.array([l[2], l[1]])*steplength)
            current.append(proj.get_vres_by_position(l[0]))

        # The rotation (from optimal superposition)
        target = np.array(target)
        current = np.array(current)
        rotationMatrix = ftms.optimal_superposition(current, target)
        if mirror:
            rotationMatrix[0,1] = -rotationMatrix[0,1]
            rotationMatrix[1,0] = -rotationMatrix[1,0]
        c_rot = np.dot(current, rotationMatrix)
        bs = fph.get_box(proj, self.scale)
        offset_centroid = target - c_rot + np.array((bs[0],bs[2]))
        offset_centroid = ftuv.get_vector_centroid(offset_centroid)
        #The rotation angle in rad"""
        angle = math.atan2(rotationMatrix[1,0], rotationMatrix[0,0])

        #Calculate the bounding square using the offset.
        bs = fph.get_box(proj, self.scale, -offset_centroid)

        return proj, angle, -offset_centroid, bs
    @profile
    def _generate_equations(self, cg):
        penalty = 0
        #Preprocess: Calculate the projection angles for the shorthening of
        # the 3 pairwise distances between landmarks.
        angles = []
        vectors3d = []
        for i, (l0, l1) in enumerate(it.combinations(self.landmarks, 2)):
            vec = cg.get_virtual_residue(l1[0], True) - cg.get_virtual_residue(l0[0], True)
            vectors3d.append(vec)
            distance3d = ftuv.magnitude(vec)
            distance2d = ftuv.vec_distance(np.array([l0[1], l0[2]]), np.array([l1[1], l1[2]]))
            #print ("distance2d = ",distance2d ,"*", self.scale,"/", self.ref_image.shape[0])
            distance2d = distance2d * self.scale/self.ref_image.shape[0] #- (0.5*self.scale/self.ref_image.shape[0])
            try:
                theta = math.acos(distance2d/distance3d)
            except ValueError:
                if distance2d>distance3d:
                    #Projected distance > real distance
                    theta = 0
                    penalty += distance2d-distance3d
                else:
                    raise
            phi = math.pi / 2 - theta # The angle between the proj.-plane normal and the 3D vector
            #print ("distance3d {}, distance 2d {}, phi {}".format(distance3d, distance2d, math.degrees(phi)))
            angles.append(math.cos(phi))
        return vectors3d, angles, penalty

class DistanceExponentialEnergy(EnergyFunction):
    _shortname = "CLA"
    _DEFAULT_CLAMP_DISTANCE = 15.0
    HELPTEXT = ("       {:3}:  Clamp elements (not nucleotides) together \n"
                "             (at ADJ={} Angstrom) with an exponential energy\n"
                "             The prefactor is used as scale (default 1).\n"
                "             Requires the elements which will be clamped together.").format(_shortname, _DEFAULT_CLAMP_DISTANCE)

    @classmethod
    def from_cg(cls, cg, prefactor, adjustment, cla_pairs, **kwargs):
        """
        :param cla_pairs: A list of 2-tuples, containing either residue numbers or cg-element names.
        """
        clamp = []
        for pair in cla_pairs:
            try:
                r1,r2 = pair
            except ValueError as e:
                print("Invalid pair '{}'".format(pair), file=sys.stderr)
                raise
            try: # initially we assume the clamp target are residue numbers
                r1=int(r1)
            except ValueError: #Or they are element names
                e1=r1
            else:
                e1=cg.get_node_from_residue_num(r1)
            try: # initially we assume the clamp target are residue numbers
                r2=int(r2)
            except ValueError: #Or they are element names
                e2=r2
            else:
                e2=cg.get_node_from_residue_num(r2)
            #e1 and e2 are element names
            if e1 not in cg.defines:
                raise ValueError("Invalid clamp pair '{}'['{}' not in cg] -'{}' "
                      ".".format(r1,e1,r2))
            if e2 not in cg.defines:
                raise ValueError("Invalid clamp pair '{}' -'{}'['{}' not in cg] "
                      ".".format(r1,r2, e2))
            if e1==e2:
                warnings.warn("Cannot clamp identical elements "
                              " {} and {} ({}=={})".format(r1,r2,e1,e2))
            else:
                clamp+=[cls(e1, e2, distance = adjustment, scale = prefactor )]
        return CombinedEnergy(clamp)

    def __init__(self, from_elem, to_elem, distance=None, scale=None):
        '''
        Create an exponential distribution for the distance between two elements.

        Here the prefactor gives the steepness of the exponential curve,
        while the adjustment gives the maximal distance that does not incur any penalty.
        '''
        if scale is None:
            scale = 1.
        if distance is None:
            distance = self._DEFAULT_CLAMP_DISTANCE
        super(DistanceExponentialEnergy, self).__init__(prefactor = scale, adjustment = distance)
        self.from_elem = from_elem
        self.to_elem = to_elem
        self.expon = scipy.stats.expon(loc=self.adjustment, scale=self.prefactor)

    def _update_pf(self):
        super(DistanceExponentialEnergy, self)._update_pf()
        self.expon = scipy.stats.expon(loc=self.adjustment, scale=self.prefactor)

    def _update_adj(self):
        super(DistanceExponentialEnergy, self)._update_adj()
        self.expon = scipy.stats.expon(loc=self.adjustment, scale=self.prefactor)

    def get_distance(self, cg):
        from_elem = self.from_elem
        to_elem = self.to_elem

        closest_points = ftuv.line_segment_distance(cg.coords[from_elem][0],
                                                   cg.coords[from_elem][1],
                                                   cg.coords[to_elem][0],
                                                   cg.coords[to_elem][1])

        closest_distance = ftuv.vec_distance(closest_points[1], closest_points[0])
        return closest_distance
    @property
    def shortname(self):
        if self.prefactor==1:
            name = "CLA"
        else:
            name = "{}CLA".format(self.prefactor)
        if self.adjustment != self._DEFAULT_CLAMP_DISTANCE:
            name = "{}{}".format(name, self.adjustment)
        return "{}({},{})".format(name,self.from_elem, self.to_elem)

    def eval_energy(self, cg, background=True, nodes=None, **kwargs):
        "This energy is always >0"
        closest_distance = self.get_distance(cg)

        if closest_distance < self.adjustment:
            return 0

        energy = -np.log(self.expon.pdf(closest_distance)) * 10.

        return energy



class StemVirtualResClashEnergy(EnergyFunction):
    '''
    Determine if the virtual residues clash.
    '''
    _shortname = "CLASH"
    _CLASH_DEFAULT_PREFACTOR = 50000.
    _CLASH_DEFAULT_ATOM_DIAMETER = 1.8
    IS_CONSTRAINT = True
    HELPTEXT = "       {:3}:  Clash constraint energy".format(_shortname)


    def __init__(self, clash_penalty = None, atom_diameter = None):
        """
        :param clash_penalty: The energy attributed to each pair of clashing atoms
        :param atom_radius: The distance between two atoms which counts as a clash
        """
        if clash_penalty is None:
            clash_penalty = self._CLASH_DEFAULT_PREFACTOR
        if atom_diameter is None:
            atom_diameter = self._CLASH_DEFAULT_ATOM_DIAMETER
        # NOTE: The atom-diameter is implemented using the adjustment-variable!
        super(StemVirtualResClashEnergy, self).__init__(prefactor = clash_penalty, adjustment = atom_diameter)
        self.bad_bulges = []
        self.bad_atoms = defaultdict(list)

    def _virtual_residue_atom_clashes_kd(self, cg):
        '''
        Check if any of the virtual residue atoms clash.
        '''
        virtual_atoms = []
        coords = []
        for key1 in self.vras.keys():
            for key2 in self.vras[key1].keys():
                virtual_atoms += [(self.vras[key1][key2], key1, key2)]
                coords += [self.vras[key1][key2]]

        if len(virtual_atoms) == 0:
            return 0

        #coords = np.vstack([p[0] for p in virtual_atoms])
        #coords = np.array([ line for line in np.array(virtual_atoms)[:,0]])
        coords = np.array(coords)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            kdt2 = kd.KDTree(3)
            kdt2.set_coords(coords)
            kdt2.all_search(self.adjustment) #Distance in Angstrom. #1.8

        clashes = 0
        indeces = kdt2.all_get_indices()
        for (ia,ib) in indeces:
            if virtual_atoms[ia][1][0] == virtual_atoms[ib][1][0]:
                continue
            if virtual_atoms[ia][1] == virtual_atoms[ib][1]:
                continue

            key1 = virtual_atoms[ia][1]
            key2 = virtual_atoms[ib][1]

            resn1 = cg.stem_side_vres_to_resn(key1[0], key1[2], key1[1])
            resn2 = cg.stem_side_vres_to_resn(key2[0], key2[2], key2[1])

            #Adjacent residues cannot clash
            if abs(resn1 - resn2) == 1:
                continue
            clash_pair = tuple(sorted([key1[0], key2[0]]))
            if clash_pair not in self.bad_bulges:
                self.bad_bulges.append(clash_pair)
            self.bad_atoms[key1[0]].append(virtual_atoms[ia][0])
            self.bad_atoms[key2[0]].append(virtual_atoms[ib][0])
            clashes += 1

        return clashes

    def _virtual_residue_atom_clashes(self, cg, s1,i1,a1, s2, i2, a2):
        '''
        Check if any of the virtual residue atoms clash.

        This is a reference implementation without KD-Trees.
        Thus it is only used for testing!
        '''
        #(p1, v1, v1_l, v1_r) = ftug.virtual_res_3d_pos(bg, s1, i1)
        #(p2, v2, v2_l, v2_r) = ftug.virtual_res_3d_pos(bg, s2, i2)


        vra1 = self.vras[(s1,i1,a1)]
        vra2 = self.vras[(s2,i2,a2)]

        clashes = 0

        atoms1 = vra1
        atoms2 = vra2

        #for atoms1 in vra1:
            #for atoms2 in vra2:

        for a1 in atoms1.values():
            for a2 in atoms2.values():
                if ftuv.vec_distance(a1,a2) < 1.8:
                #if ftuv.magnitude(a1 - a2) < 1.8:
                    clashes += 1

        return clashes

    def eval_energy(self, cg, background=False, nodes = None, **kwargs):
        '''
        Count how many clashes of virtual residues there are.

        @param sm: The SpatialModel containing the list of stems.
        @param background: Use a background distribution to normalize this one.
                           This should always be false since clashes are independent
                           of any other energies.
        '''
        #: A dict of dicts. The first key is a triple (stem, a, b), e.g.: ('s27', 5, 1)
        #: Where a is the position within the strand and b is the stem (0 or 1)
        #: The key of the inner dict is the atom, e.g. "O3'"
        self.vras = dict()
        self.bases = dict()
        self.bad_bulges = []
        self.bad_atoms = defaultdict(list)
        self.last_clashes=[]
        mult = 8
        points = []
        energy = 0.

        if nodes is None:
            nodes = cg.defines.keys()

        if len([stem for stem in nodes if stem[0]=="s"])<2:
            # Special case, if only one stem is present.
            return 0.

        for d in nodes:
            if d[0] == 's':
                s = d
                s_len = cg.stem_length(s)
                #stem_inv = cg.stem_invs[s]
                log.debug("s %s, slen %s", s, s_len)
                for i in range(s_len):
                    (p, v, v_l, v_r) = cg.v3dposs[d][i]

                    points += [(p+ mult * v_l, d, i, 1)]
                    points += [(p+ mult * v_r, d, i, 0)]
        log.debug("points: %s, nodes %s, defines %s", points, nodes, cg.defines)
        coords = np.vstack([point[0] for point in points])
        clash_pairs = []
        #kk = ss.KDTree(np.array(l))
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            kdt = kd.KDTree(3)
            kdt.set_coords(coords)
            kdt.all_search(10.)
        #print len(kdt.all_get_indices())
        #print len(kk.query_pairs(7.))
        indeces = kdt.all_get_indices()
        for (ia,ib) in indeces:
            (s1,i1,a1) = (points[ia][1], points[ia][2], points[ia][3])
            (s2,i2,a2) = (points[ib][1], points[ib][2], points[ib][3])
            clash_pairs += [((s1,i1,a1), (s2,i2,a2))]
            #print("NO NEW NODES: ", sorted(clash_pairs))

        #potential_clashes = 0
        for (s1, i1, a1), (s2,i2,a2) in clash_pairs:
            #if new_nodes != None:
            #    if s1 not in new_nodes and s2 not in new_nodes:
            #        continue

            if s1 == s2:
                continue

            if len(set.intersection(cg.edges[s1], cg.edges[s2])) > 0:
                # the stems are connected
                continue

            #potential_clashes += 1
            #fud.pv('s1,s2')

            if (s1,i1,a1) not in self.vras.keys():
                self.vras[(s1,i1,a1)] = ftug.virtual_residue_atoms(cg, s1, i1, a1)
            if (s2,i2,a2) not in self.vras.keys():
                self.vras[(s2,i2,a2)] = ftug.virtual_residue_atoms(cg, s2, i2, a2)

        energy += self.prefactor * self._virtual_residue_atom_clashes_kd(cg)
        return energy

class RoughJunctionClosureEnergy(EnergyFunction):
    _shortname = "JUNCT"
    _JUNCTION_DEFAULT_PREFACTOR = 50000.
    IS_CONSTRAINT = True
    HELPTEXT = "       {:3}:  Junction constraint energy".format(_shortname)
    @classmethod
    def from_cg(cls, cg, prefactor, adjustment, **kwargs):
        if adjustment is not None:
            warnings.warn("Adjustment {} ignored for JUNCT JunctionConstraintEnergy".format(adjustment))
        return cls(prefactor)
    def __init__(self, prefactor = None):
        if prefactor is None:
            prefactor = self._JUNCTION_DEFAULT_PREFACTOR
        super(RoughJunctionClosureEnergy, self).__init__(prefactor = prefactor)

    def eval_energy(self, cg, background=True, nodes=None, **kwargs):
        if nodes == None:
            nodes = cg.defines.keys()

        self.bad_bulges = []
        all_bulges = set([d for d in nodes if (d[0] == 'm' and d in nodes)])
        energy = 0.

        for bulge in all_bulges:
            bl = cg.get_bulge_dimensions(bulge)[0]
            dist = ftug.junction_virtual_atom_distance(cg, bulge)
            #
            #cutoff_distance = (bl) * 5.9 + 13.4
            #cutoff_distance = (bl) * 5.908 + 11.309
            #cutoff_distance = (bl) * 6.4 + 6.4
            cutoff_distance = (bl) * 6.22 + 14.0 #Peter's cyclic coordinate descent
            # Note: DOI: 10.1021/jp810014s claims that a typical MeO-P bond is 1.66A long.

            if (dist > cutoff_distance):
                log.info("Junction closure: dist {} > cutoff {} for bulge {} with length {}".format(dist, cutoff_distance, bulge, bl))
                self.bad_bulges += cg.find_bulge_loop(bulge, 200) + [bulge]
                energy += (dist - cutoff_distance) * self.prefactor

        return energy


def _iter_subgraphs(cg, use_subgraphs):
    """
    Generate a number of subgraphs. Used by generate_target_distribution of energies.

    :param cg: The CoarseGrainRNA
    :param use_subgraphs: A bool or int. See documentation of the energy

    :yields: A tuple (domain, nt_length)
    """


    if not use_subgraphs:
        return
    if len(cg.defines)<=5:
        return
    if isinstance(use_subgraphs, bool):
        target_range = it.repeat(range(3, len(cg.defines)-2), 100) # range starts at 3, because this is the
                                                               # smallest length containing 2 stems.
                                                               # For the same reason, we use len(defines)-2
                                                               # to exclude the complete RNA
    else:
        # use_subgraphs is an integer, giving the numbers of subgraphs per structure.
        target_range = []
        for i in range(use_subgraphs):
            target_range.append(random.randint(3, len(cg.defines.keys())-2))

    known_sgs = set([tuple(sorted(cg.defines.keys()))])
    for l in  target_range:
        subgraph =  cg.random_subgraph(l)
        assert len(subgraph) == len(set(subgraph))
        subgraph_t = tuple(sorted(subgraph))
        if subgraph_t not in known_sgs: #No duplicates. We sample without replacement
            known_sgs.add(subgraph_t)
            nt_length = sum(cg.element_length(elem) for elem in subgraph)
            yield subgraph, nt_length

class RadiusOfGyrationEnergy(CoarseGrainEnergy):
    _shortname = "ROG"
    HELPTEXT = "       {:3}:  Radius of gyration energy".format(_shortname)
    real_stats_fn = op.expanduser('stats/rog_target_dist_nr2.110.csv')
    sampled_stats_fn = op.expanduser('stats/rog_reference_dist_nr2.110.csv')

    def __init__(self, rna_length, adjustment=None, prefactor=None):
        """
        :param rna_length: The length in nucleotides of the RNA
        """
        super(RadiusOfGyrationEnergy, self).__init__(rna_length, prefactor=prefactor, adjustment = adjustment)

        if adjustment!=1:
            self._adjust_target_distribution()

    @classmethod
    def generate_target_distribution(cls, cg_filenames, out_filename=None, use_subgraphs = False):
        """
        Generate the target distribution from the given cg-files and write it to a file

        .. warning::

            If out_filename is not given, this overwrites the file with the name given in
            `cls.real_stats_fn`.

        :param cg_filenames: A filename or a list of filenames containing true RNA tertiary structures.
                             Typically these cg files have been generated from the pdb-files
                             using the script `pdb_to_cg.py` provided with forgi.
        :param out_filename: None or a filename relative to fess/
        :param use_subgraphs: Include the radius of subgraphs of the cg-files to get
                            more datapoints.
                            Either a bool or an integer.
                            If it is a bool, generate 100 subgraphs for each number of cg-elements possible.
                            If it is an integer: generate that many subgraphs with random length.
        """
        if out_filename is None:
            out_filename = cls.real_stats_fn
        else:
            out_filename = op.join("fess", out_filename)
        radii = []
        log.info("Generating target distribution for %s", cls.__name__)
        if isinstance(cg_filenames, str):
            cg_filenames = [cg_filenames]
        for fname in cg_filenames:
            log.debug("Processing file %s", fname)
            cg = ftmc.CoarseGrainRNA(fname)
            radii.append((cg.name, cg.seq_length, cg.radius_of_gyration("vres")))
            for subgraph, nt_length in _iter_subgraphs(cg, use_subgraphs):
                rog = ftmd.radius_of_gyration(cg.get_poss_for_domain([stem for stem in subgraph if stem[0]=="s"], mode = "fast"))
                radii.append((cg.name+".subgraph", nt_length, rog))
        log.info("Writing to file %s", out_filename)
        with open(out_filename, "w") as f:
            cls._print_file_header(f, cg_filenames)
            print("# use_subgraphs = {} ({})".format(use_subgraphs, type(use_subgraphs).__name__), file=f)
            for name, nt_len, rog in radii:
                print("{:s} {:d} {:.10f}".format(name, nt_len, rog), file=f)


    def _get_cg_measure(self, cg):
        return cg.radius_of_gyration("fast")

    def _get_distribution_from_file(self, filename, length):
        data = pd.read_csv(load_local_data(filename), delimiter=' ', comment="#")

        rdata = []

        distribution_upper_bound = 1.0
        distribution_lower_bound = 1.0
        target_len=500
        target_len=min(target_len,len(data[:]))
        while (len(rdata) < 500 and len(rdata)<len(data)):
            distribution_lower_bound -= INCR
            distribution_upper_bound += INCR

            mask_upper = data.iloc[:,1] < length * ( distribution_upper_bound )
            mask_lower = data.iloc[:,1] > ( distribution_lower_bound ) * length
            rdata = data.loc[ mask_upper & mask_lower ]

        rogs = rdata.iloc[:,2]
        rogs=np.array(rogs)

        return (self._get_distribution_from_values(rogs), rogs)


class NormalDistributedRogEnergy(RadiusOfGyrationEnergy):
    _shortname = "NDR"
    HELPTEXT = ("       {:3}:  Normal Distributed ROG energy.\n"
               "             Use a normal distribution for the target\n"
               "             radius of gyration with 0.77*ADJ as mean \n"
               "             and 0.23*ADJ stddev\n"
               "             [0.77 is a rough estimate for the relation between\n"
               "             perimeter and radius of gyration]\n".format(_shortname))

    def __init__(self, rna_length, adjustment, prefactor=None):
        """
        A Subclass of the Radius of Gyration energy with a normal distributed target distribution.

        The adjustment is treated as the estimated radius of the RNA. Thus the target distribution
        is a normal distribution with mean=0.77*ADJ and STDDEV=0.23*ADJ.

        Remember that the ROG for polymers is defined as:
        ROG=SQRT(1/N*Sum(r)) for N points with the mean in the origin.

        0.77 is roughly SQRT(3/5), which is the limes of the rog for m->inf many points
        in a 3D unit sphere with the nth point placed at radius (n/m)**1/3
        0.77-0.23 is roughly 1/SQRT(3), which is the limes of the rog for m->inf many points
        equally distributed along one or more diameters (in a star-like fashion)
        0.77+0.23 is 1, corresponding to all points on the surface of the sphere.

        If the radius is estimated from electron microscopy data,
        it is potentially smaller than the radius in 3D space if the shape of the RNA
        is not a perfect sphere.
        This effect is accounted for by allowing values greater than 1*ADJ with a
        non-neglectable probability

        :param rnalength: Used for initial reference distribution
        """
        super(NormalDistributedRogEnergy, self).__init__(rna_length, adjustment, prefactor)
        self.target_distribution = lambda x: np.array([scipy.stats.norm(loc=0.77*self.adjustment, scale=0.23*self.adjustment).pdf(x)])
        self.target_values = None

    def _update_adj(self):
        super(RadiusOfGyrationEnergy, self)._update_adjustment()
        self.target_distribution = lambda x: np.array([scipy.stats.norm(loc=0.77*self.adjustment, scale=0.23*self.adjustment).pdf(x)])


class AMinorEnergy(CoarseGrainEnergy):
    _shortname = "AME"
    HELPTEXT = "       {:3}:  A-Minor energy".format(_shortname)
    real_stats_fn = 'stats/ame_target_dist_nr2.110.csv'
    sampled_stats_fn = 'stats/ame_reference_dist_nr2.110.csv'
    orientation_file = 'stats/ame_orientation_nr2.110.csv'
    cutoff_dist = 30 #Do not consider elements above this distance for interactions.

    @classmethod
    def from_cg(cls, cg, prefactor, adjustment, **kwargs):
        """
        Get the A-minor energies for h- and i-loops.

        :param pre: Energy prefactor
        :param adj: Adjustment

        :returns: A CombinedEnergy
        """
        ame1 = cls(cg.seq_length, loop_type = 'h', prefactor=prefactor, adjustment=adjustment)
        ame2 = cls(cg.seq_length, loop_type = 'i', prefactor=prefactor, adjustment=adjustment)
        energies = []
        if ame1._get_num_loops(cg)>0:
            energies.append(ame1)
        if ame2._get_num_loops(cg)>0:
            energies.append(ame2)
        return CombinedEnergy(energies)

    @classmethod
    def _generate_target_dist_given_orientation(cls, cgs, out_filename, orientation_infile, use_subgraphs = False, cg_filenames = []):
        """
        :param cgs: A dict {pdb_id: [cg, cg,...]}


        """
        if orientation_infile.startswith("fess"):
            orientation_infile = orientation_infile[5:]

        # Create the probability functions
        all_geometries = pd.read_csv(load_local_data(orientation_infile), delimiter=' ', comment="#")
        all_geometries = all_geometries[ all_geometries["dist"] < cls.cutoff_dist ]
        aminor_geometries = all_geometries[all_geometries["is_interaction"]]
        non_ame_geometries = all_geometries[np.invert(all_geometries["is_interaction"])]
        prob_fun= {}
        for loop_type in "himft":
          log.info("Creating probability function for %s", loop_type)
          prob_fun[loop_type] = fba.aminor_probability_function(aminor_geometries.itertuples(),
                                                             non_ame_geometries.itertuples(),
                                                             loop_type)


        log.info("OUT filename is %s", out_filename)
        #Now use prob_fun to evaluate the A-minor probability for all loops in the cgs given.
        with open(out_filename, "w") as f:
            log.info("Writing AMinor target_dist to %s", out_filename)
            cls._print_file_header(f, cg_filenames)
            print("# use_subgraphs = {} ({})".format(use_subgraphs, type(use_subgraphs).__name__), file=f)
            print("# Probabilities from {}".format(orientation_infile), file=f)

            print("# cutoff_dist = {} A".format(cls.cutoff_dist), file=f)
            print("pdb_id rna_length loop_type total_prob max_prob num_interactions", file=f)
            for cgs in cgs.values():
                for cg in cgs:
                    print("# CG:", file=f)
                    rna_length = cg.seq_length
                    cls._print_prob_lines(cg, rna_length, prob_fun, f)
                    if use_subgraphs:
                        print("# Subgraphs:", file=f)
                        log.info("Now generating AMinor for Subgraphs")
                    i=0
                    for subgraph, sg_length in _iter_subgraphs(cg, use_subgraphs):
                        i+=1
                        log.info("Subgraph of length %s (%d th sg)", sg_length, i)
                        cls._print_prob_lines(cg, sg_length, prob_fun, f, subgraph)

        log.info("Successfully generated target distribution for AMinors.")
    @classmethod
    def generate_target_distribution(cls, cg_filenames, fr3d_out, out_filename=None,
                                     orientation_outfile = None, fr3d_query = "",
                                     use_subgraphs = False):
        """
        Generate the target distribution from the given cg-files and write it to a file.
        Additionally generate a file with the relative orientation of all annotated
        A-Minor interactions.

        .. warning::

            If out_filename is not given, this overwrites the file with the name given in
            `cls.real_stats_fn`.

        ..warning::

            The FR3D query needs to contain 3 nucleotides where the first is the Adenine and
            the second and third form a basepair in a stem.

        It needs FR3D output like this::

            Filename Discrepancy       1       2       3 Cha   1-2   1-3   2-3 Con   1-2   1-3   2-3   1-2   1-3   2-1   2-3   3-1   3-2   1-2   1-3   2-1   2-3   3-1   3-2   1-2   1-3   2-3
                1S72      0.0340  A  104  A  957  U 1009 900 ----  ----   cWW  AAA  1909  1885    24                                                                             -     -     -

        :param cg_filenames: A filename or a list of filenames containing true RNA tertiary structures.
                             Typically these cg files have been generated from the pdb-files
                             using the script `pdb_to_cg.py` provided with forgi.
                             .. note::

                                the first 4 characters of the CoarseGrainRNA's name must match the
                                pdb-id as reported by FR3D.

        :param fr3d_out: A file containing the annotations found by FR3D.
        :param out_filename: None or a filename relative to `fess/`
        :param orientation_outfile: Where to store the list of relative orientations of
                                    Adenines and the Stems in A-Minor interactions and
                                    in the background. Relative to "fess/"
        :param fr3d_query: Describe the FR3D-query that you used in a string.
                           It will be added as comment to the output file.
        """
        #Code in part copied from Peter's get_emlements.py

        # If out-filenames are None, use the defaults.
        if out_filename is None:
            out_filename = cls.real_stats_fn
        else:
            out_filename = op.join("fess", out_filename)
        log.info("OUT filename is %s", out_filename)
        if orientation_outfile is None:
            orientation_outfile = cls.orientation_file
        else:
            orientation_outfile = op.join("fess", orientation_outfile)

        #Read the cgs
        all_cgs = defaultdict(list)
        for fn in cg_filenames:
            cg = ftmc.CoarseGrainRNA(fn)
            all_cgs[cg.name[:4]].append(cg)

        #Read the FR3D output
        with open(fr3d_out) as f:
            aminor_geometries = fba.parse_fred(cls.cutoff_dist, all_cgs, f)
        non_ame_geometries = set()
        for pdb_id, cgs in all_cgs.items():
            for cg in cgs:
                for loop in cg.defines:
                    if loop[0]=="s":
                        continue
                    for stem in cg.stem_iterator():
                        if loop in cg.edges[stem]:
                            continue
                        dist, angle1, angle2 = fba.get_relative_orientation(cg, loop, stem)
                        if dist<=cls.cutoff_dist and "A" in "".join(cg.get_define_seq_str(loop)) and not np.isnan(dist+angle1+angle2):
                            geometry = fba.AMGeometry(pdb_id, loop, stem, dist, angle1, angle2, "&".join(cg.get_define_seq_str(loop)))
                            if geometry in aminor_geometries:
                                log.info("Geometry %s is in aminor_geometries", geometry)
                            else:
                                non_ame_geometries.add(geometry)

        #Print orientations to orientation_outfile
        with open(orientation_outfile, "w") as f:
            cls._print_file_header(f, cg_filenames)
            print("# fr3d_out = {}".format(fr3d_out), file=f)
            print("# fr3d_query:", file=f)
            for line in fr3d_query.splitlines():
                line=line.strip()
                print("#    "+line, file = f)
            print("# cutoff_dist = {} A".format(cls.cutoff_dist), file=f)
            # HEADER
            print ("loop_type dist angle1 angle2 is_interaction", file=f)
            for entry in aminor_geometries:
                print("{pdb_id} {loop_type} {dist} {angle1} "
                      "{angle2} {is_interaction}".format(is_interaction = True,
                                                         **entry._asdict()), file = f)
            for entry in non_ame_geometries:
                print("{pdb_id} {loop_type} {dist} {angle1} "
                      "{angle2} {is_interaction}".format(is_interaction = False,
                                                         **entry._asdict()), file = f)
        cls._generate_target_dist_given_orientation(all_cgs, out_filename,
                                                    orientation_outfile, use_subgraphs, cg_filenames)



    def __init__(self, rna_length, loop_type='h', adjustment=None, prefactor=None):
        """
        :param loop_type: A one-letter string giving the loop type. E.g. 'h' or 'i'
        """
        self.loop_type = loop_type
        super(AMinorEnergy, self).__init__(rna_length, prefactor=prefactor, adjustment = adjustment)


        # Load the geometries (aminor and other)
        all_geometries = pd.read_csv(load_local_data(self.orientation_file), delimiter=' ', comment="#")
        all_geometries = all_geometries[ all_geometries["dist"] < self.cutoff_dist ]
        aminor_geometries = all_geometries[all_geometries["is_interaction"]]
        non_ame_geometries = all_geometries[np.invert(all_geometries["is_interaction"])]

        self.prob_function = fba.aminor_probability_function(aminor_geometries.itertuples(),
                                                             non_ame_geometries.itertuples(), #With pandas >=0.19 this yields namedtuples!
                                                             self.loop_type)

        #: The number of coarse grain elements considered in this energy
        self.num_loops=None

    @property
    def shortname(self):
        sn = super(AMinorEnergy, self).shortname
        return sn.replace(self._shortname, "{}({})".format(self._shortname,self.loop_type))

    def _get_distribution_from_file(self, filename, length):
        data = pd.read_csv(load_local_data(filename), delimiter=' ', comment="#")
        data = data[ data["loop_type"]==self.loop_type]
        data = data.as_matrix(["rna_length", "total_prob"])
        rdata = []

        distribution_lower_bound = 1.
        distribution_upper_bound = 1.

        while (len(rdata) < 500 and len(rdata)<len(data)):
            try:
                distribution_lower_bound -= INCR
                distribution_upper_bound += INCR

                rdata = data[np.logical_and( data[:,0] > ( distribution_lower_bound ) * length,
                                             data[:,0] < length * ( distribution_upper_bound ))]
            except KeyboardInterrupt:
                print("len(rdata) is {}, len(data)={}, bound= {}...{}".format(len(rdata),len(data),
                     distribution_lower_bound ) * length, length * ( distribution_upper_bound ))
                raise
        if len(rdata)==0:
            raise ValueError("No data found for distribution in file {}".format(filename))

        rogs = rdata[:,1]
        return (self._get_distribution_from_values(rogs), rogs)

    def eval_prob(self, cg, d):
        return fba.total_prob(d, cg, self.prob_function, self.cutoff_dist)

    def _get_cg_measure(self, cg):
        raise NotImplementedError("Not needed for A-Minor energy")

    @property
    def name(self):
        if self.loop_type == 'i':
            return "A-Minor Energy (interior loops)"
        elif self.loop_type == 'h':
            return "A-Minor Energy (hairpin loops)"
        else:
            return "A-Minor Energy"

    def accept_last_measure(self):
        """
        self._last_measure is more than one
        """
        if self._last_measure is not None:
            self.accepted_measures.extend(self._last_measure)
        self._step_complete()

    def reject_last_measure(self):
        """AMinor.rejectLastMeasure"""
        if len(self.accepted_measures) > 0 and self.num_loops>0:
            self.accepted_measures.extend(self.accepted_measures[-self.num_loops:])

    def _get_num_loops(self, cg):
        possible_loops = [d for d in cg.defines.keys() if d[0] == self.loop_type and 'A' in "".join(cg.get_define_seq_str(d))]
        return len(possible_loops)

    def eval_energy(self, cg, background=True, use_accepted_measure = False, plot_debug=False, **kwargs): #@PROFILE: This takes >50% of the runtime with default energy
        kr = self.target_distribution
        ks = self.reference_distribution

        energy = 0

        if self.num_loops is None:
            self.num_loops=self._get_num_loops(cg)

        self._last_measure = []
        for d in cg.defines.keys():

            # the loop type is encoded as an integer so that the stats file can be
            # loaded using numpy
            if d[0] != self.loop_type or 'A' not in "".join(cg.get_define_seq_str(d)):
                continue


            m = self.eval_prob(cg, d)
            self._last_measure.append(m)

            if background:
                prev_energy = (np.log(kr(m) + 0.00000001 * ks(m)) - np.log(ks(m)))
                self.prev_energy = energy
                self.prev_cg = m
                energy +=  -1 * self.prefactor * prev_energy
            else:
                energy +=  -np.log(kr(m))
        if plot_debug: #For debuging
                import matplotlib.pyplot as plt
                xs=np.linspace(0, 1, 500)
                fig,ax1 = plt.subplots()
                ax2 = ax1.twinx()
                ax1.plot(xs, ks(xs), label="referecne distribution")
                ax1.plot(xs, kr(xs), label="target distribution")
                ax2.plot(xs, -(np.log(kr(xs) + 0.00000001 * ks(xs)) - np.log(ks(xs))), label="energy", color="red")
                ax1.plot(self.accepted_measures, [1]*len(self.accepted_measures), "o", label="Accepted Measures")
                plt.title(self.shortname)
                ax1.legend(loc="lower left")
                ax2.legend()
                plt.show()
        return energy[0]

    @classmethod
    def _print_prob_lines(cls, cg, rna_length, prob_funs, file, domain = None):
        """
        Print the probabilities for all loops in this cg to the file file

        :param cg: A CoarseGrainRNA
        :param prob_funs: A dictionary {loop_type:probability_function} where
                          loop_type is a letter ("i", "m", "h", "f" or "t")
                          and probability function is a function like those
                          returned by fba.aminor_probability_function
        """
        for loop in cg.defines:
            if loop[0] == "s":
                continue
            t_prob = fba.total_prob(loop, cg, prob_funs[loop[0]], cls.cutoff_dist, domain)
            max_prob = fba.total_prob(loop, cg, prob_funs[loop[0]], cls.cutoff_dist, domain)
            num_interactions = fba.total_prob(loop, cg, prob_funs[loop[0]], cls.cutoff_dist, domain)
            if not np.isnan(t_prob*max_prob*num_interactions):
                assert max_prob<=t_prob<=num_interactions, ("{} <=? {} "
                                                           "<=? {}".format(max_prob,
                                                                           t_prob,
                                                                           num_interactions))
                print("{} {} {} {} {} {}".format(cg.name, rna_length,
                                                 loop[0], t_prob,
                                                 max_prob, num_interactions),
                      file = file)

class DoNotContribute(Exception):
    pass


def _minimal_h_h_distance(cg, elem1, elem2_iterator):
    """
    Used by ShortestLoopDistancePerLoop-Energy.

    :param cg: The CoarseGrain RNA
    :param elem1: A STRING. A name of a hairpin loop. e.g. "h1"
    :param elem2_iterator: A ITERATOR/ LIST. Element names to compare h1 with.
    """
    min_dist = float("inf")
    for elem2 in elem2_iterator:
        if elem2==elem1: continue
        (i1,i2) = ftuv.line_segment_distance(cg.coords[elem1][0],
                                            cg.coords[elem1][1],
                                            cg.coords[elem2][0],
                                            cg.coords[elem2][1])
        dist = ftuv.vec_distance(i1, i2)
        if dist < min_dist:
            min_dist = dist
    return min_dist

class ShortestLoopDistancePerLoop(CoarseGrainEnergy):
    _shortname = "SLD"
    HELPTEXT = "       {:3}:  shortest loop distance per loop".format(_shortname)
    real_stats_fn = 'stats/sld_target_dist_nr2.110.csv'
    sampled_stats_fn = 'stats/sld_reference_dist_nr2.110.csv'

    @classmethod
    def from_cg(cls, cg, prefactor, adjustment, **kwargs):
        """
        Get the shortest loopdistance per loop energy for each hloop.

        :param cg: The coarse grained RNA
        :returns: A CombinedEnergy
        """
        energies=[]
        for hloop in cg.hloop_iterator():
            energies+= [cls(rna_length = cg.seq_length, loop_name = hloop, prefactor = prefactor, adjustment = adjustment)]
        return CombinedEnergy(energies)

    @classmethod
    def generate_target_distribution(cls, cg_filenames, out_filename=None, use_subgraphs = False):
        """
        Generate the target distribution from the given cg-files and write it to a file

        .. warning::

            If out_filename is not given, this overwrites the file with the name given in
            `cls.real_stats_fn`.

        :param cg_filenames: A filename or a list of filenames containing true RNA tertiary structures.
                             Typically these cg files have been generated from the pdb-files
                             using the script `pdb_to_cg.py` provided with forgi.
        :param out_filename: None or a filename relative to fess/
        :param use_subgraphs: Include the radius of subgraphs of the cg-files to get
                            more datapoints.
                            Either a bool or an integer.
                            If it is a bool, generate 100 subgraphs for each number of cg-elements possible.
                            If it is an integer: generate that many subgraphs with random length.

        """
        if out_filename is None:
            out_filename = cls.real_stats_fn
        else:
            out_filename = op.join("fess", out_filename)


        distances = []
        log.info("Generating target distribution for %s", cls.__name__)
        if isinstance(cg_filenames, str):
            cg_filenames = [cg_filenames]

        loop_loop_distances = []
        for fname in cg_filenames:
            log.info("Processing file %s", fname)
            cg = ftmc.CoarseGrainRNA(fname)
            pdbid = cg.name
            for h1 in cg.hloop_iterator():
                min_dist = _minimal_h_h_distance(cg, h1, cg.hloop_iterator())
                if min_dist<float("inf"):
                    loop_loop_distances.append((pdbid, cg.seq_length, min_dist))
            pdbid = cg.name+".subgraph"
            for subgraph, nt_length in _iter_subgraphs(cg, use_subgraphs):
                for h1 in subgraph:
                    if h1[0]!="h":
                        continue
                    min_dist = _minimal_h_h_distance(cg, h1, (h for h in subgraph if h[0]=="h"))
                    if min_dist<float("inf"):
                        loop_loop_distances.append((pdbid, nt_length, min_dist))

        log.info("Writing to file %s", out_filename)
        with open(out_filename, "w") as f:
            cls._print_file_header(f, cg_filenames)
            print("# use_subgraphs = {} ({})".format(use_subgraphs, type(use_subgraphs).__name__), file=f)
            for pdbid, nt_len, distance in loop_loop_distances:
                print("{:s} {:d} {:.10f}".format(pdbid, nt_len, distance), file=f)


    def __init__(self, rna_length, loop_name, prefactor=None, adjustment = None):

        #: Add equally distributed points to the target and reference distribution estimation (linspacepoints  lsp)
        #: Weight of the true data compared to the artificial points (integer)
        self._lsp_data_weight = 3
        #: Weight of the artificial data (usually 1 or 0)
        self._lsp_artificial_weight = 1
        #: How many artificial points to add to the reference distribution
        self._lsp_reference_num_points = 30
        #: How many artificial points to add to the target distribution
        self._lsp_target_num_points = 15
        #: Start of the range from which to create artificial ppoints
        self._lsp_min = 3
        #: End of the range for artificial points
        self._lsp_max = 300

        super(ShortestLoopDistancePerLoop, self).__init__(rna_length, prefactor, adjustment)
        self.loop_name = loop_name

    @property
    def shortname(self):
        sn = super(ShortestLoopDistancePerLoop, self).shortname
        return sn.replace(self._shortname, "{}({})".format(self._shortname,self.loop_name))

    def _get_distribution_from_file(self, filename, length):
        data = pd.read_csv(load_local_data(filename), delimiter=' ', comment="#")
        data=data.iloc[:,-1].as_matrix() #Ignore the nt-length and the pdb-id. Only look at the distances.
        if filename == self.sampled_stats_fn:
            lsp=np.linspace(self._lsp_min, self._lsp_max,
                            num=self._lsp_reference_num_points )
        else:
            lsp=np.linspace(self._lsp_min, self._lsp_max,
                            num=self._lsp_target_num_points )

        data=np.concatenate([data]*self._lsp_data_weight+[lsp]*self._lsp_artificial_weight)

        return (self._get_distribution_from_values(data), data)

    def _get_cg_measure(self, cg):
        min_dist = _minimal_h_h_distance(cg, self.loop_name, cg.hloop_iterator())
        return min_dist
    def eval_energy(self, cg, background=True, nodes=None, **kwargs):
        '''
        We want to return an energy of 0. if there's less than two hairpin
        loops.
        '''
        if len(list(cg.hloop_iterator())) < 2:
            return 0.
        else:
            try:
                energy= super(ShortestLoopDistancePerLoop, self).eval_energy(cg,
                                                                         background,
                                                                         nodes, **kwargs)
            except DoNotContribute:
                energy=0
            return energy

class CombinedFunction(object):
    def __init__(self, funcs):
        self._funcs = funcs
    def __call__(self, *args, **kwargs):
        results = []
        for fun in self._funcs:
            results.append(fun(*args, **kwargs))
        log.debug("Combined Function called. Results are {}".format(results))
        if not results or results[0] is None:
            return None
        return results

class CombinedEnergy(object):
    def __init__(self, energies=None, normalize=False):
        """
        :param normalize: Divide the resulting energy by the numbers of contributions
        """
        if energies is not None:
            super(CombinedEnergy, self).__setattr__("energies", energies)
        else:
            super(CombinedEnergy, self).__setattr__("energies", [])
        super(CombinedEnergy, self).__setattr__("constituing_energies", [])
        super(CombinedEnergy, self).__setattr__("normalize", normalize)

    def __setattr__(self, name, val):
        if name not in self.__dict__:
            for e in self.energies:
                setattr(e, name, val)
        else:
            self.__dict__[name]=val

    def __getattr__(self, name):
        #If no energies are present, do nothing (DON'T raise an error)
        if not self.energies:
            return CombinedFunction([])
        log.debug("getattr delegating call to {}".format(name))
        attrs = []
        for e in self.energies:
            try:
                attrs.append(getattr(e, name))
            except AttributeError as e:
                if hasattr(CoarseGrainEnergy, name):
                    # If the function only exists for CoarseGrainedEnergies,
                    # but not for normal Energy functions,
                    # call it only on CoarseGrainedEenrgies.
                    pass
                else:
                    # Else this is an error.
                    raise AttributeError("CombinedEnergy has no attribute '{}', because {}".format(name, e))
        try:
            return sum(attrs)
        except:
            pass
        try:
            return "".join(attrs)
        except:
            pass
        return CombinedFunction(attrs)

    def uses_background(self):
        for e in self.energies:
            if isinstance(e, CoarseGrainEnergy):
                return True
        return False

    def iterate_energies(self):
        """
        Iterate over all member enegies
        in a way that flattens nested CombinedEnergies.
        """
        for e in self.energies:
            if hasattr(e, "iterate_energies"):
                for e2 in e.iterate_energies():
                    yield e2
            else:
                yield e

    @property
    def shortname(self):
        name=[]
        for e in self.energies:
            name.append(e.shortname)
        sn=",".join(name)

        return sn

    @property
    def bad_bulges(self):
        bad_bulges = []
        for e in self.energies:
            if hasattr(e, "bad_bulges"):
                bad_bulges+=e.bad_bulges
        return bad_bulges
    def eval_energy(self, cg, verbose=False, background=True,
                    nodes=None, use_accepted_measure=False, plot_debug=False):
        total_energy = 0.
        self.constituing_energies=[]
        num_contribs=0

        for energy in self.energies:
            contrib = energy.eval_energy(cg, background=background, nodes=nodes,
                                         use_accepted_measure=use_accepted_measure, plot_debug = plot_debug)

            if not np.isscalar(contrib):
                raise TypeError
                contrib, = contrib

            self.constituing_energies.append((energy.shortname, contrib))
            total_energy += contrib
            num_contribs +=1

            if verbose:
                print (energy.name, energy.shortname, contrib)
                if energy.bad_bulges:
                    print("bad_bulges:", energy.bad_bulges)
            log.debug("{} ({}) contributing {}".format(energy.__class__.__name__, energy.shortname, contrib))


        if num_contribs>0:
            if self.normalize:
                total_energy=total_energy/num_contribs
        else:
            assert self.energies == []

        if verbose:
            print ("--------------------------")
            print ("total_energy:", total_energy)
        return total_energy

    def __str__(self):
        out_str = ''
        for en in self.energies:
            out_str += en.__class__.__name__ + " "
        return out_str

    def __len__(self):
        return len(self.energies)

    def hasinstance(self, cls):
        """
        Returns True, if any of the member energies is an instance of the given class
        """
        for e in self.energies:
            if hasattr(e, "hasinstance"):
                if e.hasinstance(cls): return True
            elif isinstance(e, cls):
                return True
        return False
####################################################################################################
## Convenience functions for creating energies
####################################################################################################


def energies_from_string(contribution_string, cg, num_steps = None, **kwargs):
    """
    :param contribution_string: A string with comma-seperated energy contributions.
    :param cg: The CoarseGrainRNA
    :param num_steps: The number of total steps of the simulation. (Used only for simulated annealing. Else None)
    :param kwargs: The keyword args will be passed on to the energie's from_cg methos.
    """
    contributions = contribution_string.split(",")
    energy_classes = { cls._shortname: cls for cls in get_all_subclasses(EnergyFunction) if not cls == CoarseGrainEnergy }
    energies = []
    for contrib in contributions:
        match = re.match(r"([^A-Z]*)([A-Z]+)(.*)", contrib)
        if match is None:
            raise ValueError("Contribution {} not understood".format(contrib))
        cls = energy_classes[match.group(2)]
        pre = _parseEnergyContributionString(match.group(1), num_steps)
        adj = _parseEnergyContributionString(match.group(3), num_steps)
        try:
            energies.append(cls.from_cg(cg, pre, adj, **kwargs))
        except TypeError as e:
            if "arguments" not in e.message:
                raise
            # http://stackoverflow.com/a/2677263
            args = inspect.getargspec(cls.from_cg).args
            missing_arg = set(args)-set(["cls", "cg", "prefactor", "adjustment"]+list(kwargs.keys()))
            raise TypeError("The following required Keyword-Arguments for energy {} are missing: {}".format(match.group(2), list(missing_arg)))
    return CombinedEnergy(energies)

def _parseEnergyContributionString(contrib, num_steps):
    """
    Helper function
    """
    if not contrib:
        return None
    if "_" in contrib:
        a=contrib.split("_")
        if len(a)==2:
            start=float(a[0])
            end=float(a[1])
            if abs(start-end)<1.2:
                step=0.1
            elif abs(start-end)>100:
                step=10
            else:
                step=1
            if end<start:
                step=step*-1
        elif len(a)==3:
            start=float(a[0])
            step=float(a[1])
            end=float(a[2])
        else:
            raise ValueError("Too many underscores in {}".format(contrib))
        frequency=num_steps / (math.ceil((end-start)/step)+1)
        assert frequency>1, numSteps

        if frequency>=num_steps:
            raise ValueError("Could not parse energy program '{}': "
                             "Expected START_STEP_STOP, found {}, {}, {} which would "
                             "lead to a change every {} simulation steps".format(contrib, start, step,
                                                                       end, frequency))
        return (start, step, frequency)
    else:
        return float(contrib)

def get_argparse_help():
    """
    Return a pre-formatted string that can be passed to "help" in argparse.add_argument,
    if the resulting argument is parsed with `energies_from_string`
    """
    help= ["Specify a ','-separated list of energy contributions.\n"
           "Each contribution has the format: [PRE]TYP[ADJ].\n"
           "PRE: optional energy prefactor\n"
           "ADJ: optional adjustment of target distribution\n"
           "     (float), default=1.0\n"
           "For simulated annealing, ADJ and PRE can be changed \n"
           "     during the simulation. To achieve this, give \n"
           "     underscore-seperated ranges START_END or START_STEP_END\n"
           "     as PRE and/or ADJ. E.g. 1.0_0.1_1.4\n"
           "     For each prefactor/ adjustment, equally \n"
           "     many sampling steps are used.\n"
           "     If step is not given, 1.0 or 0.1 is used, \n"
           "     depending on the difference between START and END.\n"
           "TYP: one of the following"]
    for cls in get_all_subclasses(EnergyFunction):
        if inspect.isabstract(cls):
            continue
        if hasattr(cls, "IS_CONSTRAINT") and cls.IS_CONSTRAINT:
            continue
        help.append(cls.HELPTEXT)
    return "\n".join(help)
