
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
from ..aux.utils import get_all_subclasses
import os.path as op
import pandas as pd
import pkgutil as pu
import StringIO
import math
import re
import sys
import inspect
from pprint import pprint
import logging
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


class RadiusOfGyrationEnergy(CoarseGrainEnergy):
    _shortname = "ROG"
    HELPTEXT = "       {:3}:  Radius of gyration energy".format(_shortname)

    def __init__(self, rna_length, adjustment=None, prefactor=None):
        """
        :param rna_length: The length in nucleotides of the RNA
        """
        self.sampled_stats_fn = op.expanduser('stats/subgraph_radius_of_gyration_sampled.csv')        
        self.real_stats_fn = op.expanduser('stats/subgraph_radius_of_gyration.csv')
        super(RadiusOfGyrationEnergy, self).__init__(rna_length, prefactor=prefactor, adjustment = adjustment)

        if adjustment!=1:
            self._adjust_target_distribution()
                    
    def _get_cg_measure(self, cg):
        return cg.radius_of_gyration("vres")

    def _get_distribution_from_file(self, filename, length):
        data = np.genfromtxt(load_local_data(filename), delimiter=' ')

        rdata = []

        distribution_upper_bound = 1.0
        distribution_lower_bound = 1.0
        target_len=500
        target_len=min(target_len,len(data[:]))
        while (len(rdata) < target_len):
            distribution_lower_bound -= INCR
            distribution_upper_bound += INCR

            rdata = data[np.logical_and( data[:,0] > ( distribution_lower_bound ) * length,
                                         data[:,0] < length * ( distribution_upper_bound ))]

        rogs = rdata[:,1]
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

    def __init__(self, rna_length, loop_type='h', adjustment=None, prefactor=None):
        self.real_stats_fn = 'stats/aminors_1s72.csv'
        self.sampled_stats_fn = 'stats/aminors_1jj2_sampled.csv'        
        
        self.types = {'h':0, 'i':1, 'm':2, 's': 3, 'f': 4, 't': 5}
        self.loop_type = self.types[loop_type]

        super(AMinorEnergy, self).__init__(rna_length, prefactor=prefactor, adjustment = adjustment)

        dall = pd.read_csv(load_local_data('stats/tall.csv'), delimiter=' ', 
                           names=['l1','l2','dist','angle', "angle2",'seq']) 
        dall = dall[dall["dist"] < 30]

        # load the loop-anything annotations filter to consider only those
        # that contain an A and only those that are within 30 angstroms
        # of each other
        ael = pd.read_csv(load_local_data('stats/all_elements.csv'), delimiter=' ', 
                          names=['l1','l2','dist','angle',"angle2",'seq'])
        dbg = ael[["A" in x for x in ael["seq"]]]
        dbg_close = dbg[dbg["dist"] < 30.]

        self.dall = dict()
        self.dbg_close = dict()

        self.prob_funcs = dict()
        p_d_a_a2_given_i = dict()
        p_d_a_a2 = dict()
        p_i = dict()

        for lt in ['i', 'h']:
            dl = self.dall[lt] = dall[[lt in x for x in dall['l1']]]
            db = self.dbg_close[lt] = dbg_close[[lt in x for x in dbg_close['l1']]]
            p_i[lt] = len(dl['dist']) / float(len(db['dist']))


            d_for_p_d_a_a2_given_i = np.array(list(zip(dl["dist"], dl["angle"], dl["angle2"])))
            p_d_a_a2_given_i[lt] = scipy.stats.gaussian_kde(d_for_p_d_a_a2_given_i.T)
            p_d_a_a2[lt] = scipy.stats.gaussian_kde(np.array(list(zip(db["dist"], db["angle"], db["angle2"]))).T)
            #self.prob_funcs[lt] = lambda point: (p_d_a_a2_given_i[lt](point) * p_i[lt]) / (p_d_a_a2[lt](point))
            self.prob_funcs[lt] = lambda point: (p_d_a_a2_given_i[lt](point)) / (p_d_a_a2[lt](point) + p_d_a_a2_given_i[lt](point))
        #print ("SELF>MEASURES ({}): {}".format(self.shortname(),len(self.measures)))
        #print (self.measures)

        #: The number of coarse grain elements considered in this energy
        self.num_loops=None

    @property
    def shortname(self):
        sn = super(AMinorEnergy, self).shortname
        return sn.replace(self._shortname, "{}({})".format(self._shortname,self.loop_type))

    def _get_distribution_from_file(self, filename, length):
        data = np.genfromtxt(load_local_data(filename), delimiter=' ')

        rdata = list()

        distribution_lower_bound = 1.
        distribution_upper_bound = 1.

        while (len(rdata) < 500):
            distribution_lower_bound -= INCR
            distribution_upper_bound += INCR

            rdata = data[np.logical_and( data[:,0] > ( distribution_lower_bound ) * length,
                                         data[:,0] < length * ( distribution_upper_bound ))]
        
        srdata = rdata[rdata[:,1] == self.loop_type]
        rogs = srdata[:,2]
        return (self._get_distribution_from_values(rogs), rogs)

    def eval_prob(self, cg, d):
        lt = d[0]
        prob = 0.
        stem_counts = 0
        probs = []

        for s in cg.stem_iterator():
            if s in cg.edges[d]:
                continue

            if not ftuv.elements_closer_than(cg.coords[d][0],
                                             cg.coords[d][1],
                                             cg.coords[s][0],
                                             cg.coords[s][1], 
                                             30):
            #if ftug.element_distance(cg, d, s) > 30:
                continue

            point = fba.get_relative_orientation(cg, d, s)
            stem_counts += 1
            p = self.prob_funcs[lt](point)
            prob += p
            probs += [p]
            #print("dsp:",d,s,p)
        if len(probs) == 0:
            #print("eval_prob: Returning Array with one zero.")
            return np.array([0.])
        #print ("eval_prob: Returning max(probs).", max(probs))
        return max(probs)
        #return (prob, stem_counts)
        #return prob / stem_counts

    def _get_cg_measure(self, cg):
        for d in cg.defines.keys():

            # the loop type is encoded as an integer so that the stats file can be 
            # loaded using numpy
            if self.types[d[0]] != self.loop_type or 'A' not in "".join(cg.get_define_seq_str(d)):
                continue

            m = self.eval_prob(cg, d)[0]

            return m
        
    @property
    def name(self):
        if self.loop_type == self.types['i']:
            return "A-Minor Energy (interior loops)"
        elif self.loop_type == self.types['h']:
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
        possible_loops = [d for d in cg.defines.keys() if self.types[d[0]] == self.loop_type and 'A' in "".join(cg.get_define_seq_str(d))]
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
            if self.types[d[0]] != self.loop_type or 'A' not in "".join(cg.get_define_seq_str(d)):
                continue


            m = self.eval_prob(cg, d)[0]
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

class DoNotContribute(Exception):
    pass

class ShortestLoopDistancePerLoop(CoarseGrainEnergy):
    _shortname = "SLD"
    HELPTEXT = "       {:3}:  shortest loop distance per loop".format(_shortname)

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
    
    def __init__(self, rna_length, loop_name, prefactor=None, adjustment = None):        
        self.real_stats_fn = 'stats/loop_loop3_distances_native.csv'
        self.sampled_stats_fn = 'stats/loop_loop3_distances_sampled.csv'
        
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
        data = np.genfromtxt(load_local_data(filename), delimiter=' ')
        if filename == self.sampled_stats_fn:
            lsp=np.linspace(self._lsp_min, self._lsp_max, 
                            num=self._lsp_reference_num_points )
        else:
            lsp=np.linspace(self._lsp_min, self._lsp_max, 
                            num=self._lsp_target_num_points )
                            
        data=np.concatenate([data]*self._lsp_data_weight+[lsp]*self._lsp_artificial_weight)
        
        return (self._get_distribution_from_values(data), data)

    def _get_cg_measure(self, cg):
        min_dist = float("inf")

        for h in cg.hloop_iterator():
            # don't want the distance to itself
            if h == self.loop_name:
                continue

            (i1,i2) = ftuv.line_segment_distance(cg.coords[self.loop_name][0],
                                              cg.coords[self.loop_name][1],
                                              cg.coords[h][0],
                                              cg.coords[h][1])

            dist = ftuv.vec_distance(i1, i2)
            if dist < min_dist:
                min_dist = dist
        #if min_dist>60:
            #raise DoNotContribute
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
