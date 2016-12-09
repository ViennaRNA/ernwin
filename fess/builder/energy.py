from __future__ import print_function, division


import StringIO
import pickle
import os
import copy
import itertools as it
import warnings
import pkgutil as pu
#import pandas as pa
#import pylab

import Bio.KDTree as kd
import numpy as np
import random as rand

import collections as c
import os.path as op

#import scipy.stats as ss
import forgi.threedee.utilities.vector as ftuv
import forgi.threedee.model.coarse_grain as ftmc
import forgi.threedee.utilities.graph_pdb as ftug
import fess.builder.aminor as fba
import fess.builder.models as cbm
import forgi.utilities.debug as fud
import forgi.projection.hausdorff as fph
import forgi.projection.projection2d as fpp
import forgi.threedee.model.similarity as ftms
import forgi.threedee.model.descriptors as ftmd

from . import config

import scipy.stats as stats
import scipy.stats as ss
import scipy.optimize
import scipy.ndimage
import scipy.misc
import sys, math

import gc
import os.path
import time, random
#from fess.builder.watcher import watcher
import functools
import os

import logging
log = logging.getLogger(__name__)

try:
  profile  #The @profile decorator from line_profiler (kernprof)
except:
  def profile(x): 
    return x


distribution_upper_bound = 1.0
distribution_lower_bound = 1.0

INCR = 0.01

DEFAULT_ENERGY_PREFACTOR = 30

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


class DistanceIterator:
    '''
    A class for iterating over the elements that are a certain
    distance from each other.

    '''
    def __init__(self, min_distance=0., max_distance=sys.float_info.max):
        '''
        @param min_distance: The minimum distance for an interaction.
        @param max_distance: The maximum distance for an interaction.
        '''
        self.min_distance = min_distance
        self.max_distance = max_distance

    def iterate_over_interactions(self, bg):
        '''
        Iterate over the list of elements in a structure that interact
        based on the distance criteria below.

        @param bg: The coarse grain representation of a structure.
        @return: (d1,d2) where d1 and d2 are two interacting elements.
        '''
        keys = bg.defines.keys()

        for i in range(len(keys)):
            for j in range(i + 1, len(keys)):
                d1 = keys[i]
                d2 = keys[j]

                point1 = bg.get_point(d1)
                point2 = bg.get_point(d2)

                dist = ftuv.vec_distance(point1, point2)

                #if dist > 6.0 and dist < 25.0:
                if dist > self.min_distance and dist < self.max_distance:
                    yield tuple(sorted([d1, d2]))

class EnergyFunction(object):
    '''
    The base class for energy functions.
    '''

    def __init__(self):
        #: Used by constraint energies, to store tuples of stems that clash.
        #: Updated every time eval_energy is called.
        self.bad_bulges = []

        self.measures = []
        self.accepted_measures = []
    
    def accept_last_measure(self):
        """EnergyFunction.acceptLastMeasure"""
        if len(self.measures) > 0:        
            self.accepted_measures.append(self.measures[-1])

    def update_adjustment(*args, **kwargs):
        pass

    def reject_last_measure(self):
        """EnergyFunction.rejectLastMeasure"""
        if len(self.accepted_measures) > 0:
            self.accepted_measures += [self.accepted_measures[-1]]

    def eval_energy(self, sm, background=True, nodes=None):
        '''
        The base energy function simply returns a random number.
        '''
        warnings.warn("WARNING: using a random energy.")
        return rand.random()

    def resample_background_kde(self, struct):
        pass

    def calc_fg_parameters(self, bg):
        '''
        Found out the parameters for the target distribution.

        In many cases these will be derived from the list of statistics.
        '''
        pass

    def calc_bg_parameters(self, energy_structs):
        '''
        Adjust the parameters based on the sampling statistics from
        a previous run.

        @param energy_structs: A list of tuples of the form (energy, bg)
        '''
        raise NotImplementedError("TODO")

    def get_energy_name(self):
        '''
        Return the name of the energy.
        '''
        return self.__class__.__name__.lower()

    def shortname(self):
        return self.__class__.__name__.lower()[:9]+"..."
    def dump_measures(self, base_directory, iteration=None):
        '''
        Dump all of the coarse grain measures collected so far
        to a file.
        '''
        output_file = op.join(base_directory, self.get_energy_name()+str(hex(id(self)))+".measures")
        with open(output_file, 'w') as f:
            f.write(" ".join(map("{:.2f}".format,self.accepted_measures)))
            f.write("\n")

        if iteration is not None:
            with open(output_file + ".%d" % (iteration), 'w') as f:
                f.write(" ".join(map("{:.2f}".format,self.accepted_measures)))
                f.write("\n")

class ConstantEnergy(EnergyFunction):
    '''
    An energy function that always just returns a constant value (0.).
    '''
    def __init__(self):
        super(ConstantEnergy, self).__init__()

    
    def eval_energy(self, sm, background=True, nodes=None):
        return 0.

class RandomEnergy(EnergyFunction):
    '''
    An energy function that always just returns a random value.
    '''
    def __init__(self):
        super(RandomEnergy, self).__init__()

    
    def eval_energy(self, sm, background=True, nodes=None):
        return rand.uniform(-5, 3)


class ProjectionMatchEnergy(EnergyFunction):
    def __init__(self, distances={}, prefactor=1):
        """
        :param directions: A dict where the keys are tuples of coarse grain element names 
                           (e.g.: ("h1","m1")) and the values are the distance 
                           IN THE PROJECTED PLANE (i.e. in the micrograph).
        """        
        super(ProjectionMatchEnergy, self).__init__()
        self.distances=distances
        self.prefactor=prefactor
        #: The optimal projection direction for the last accepted step.
        self.accepted_projDir=np.array([1.,1.,1.])
        self.projDir=np.array([1.,1.,1.])
        self.start_points=self.get_start_points(60)
    def get_start_points(self, numPoints):
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

    def shortname(self):
        return "PRO"

    def get_name(self):
        return "Projection Matching Energy"

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

    def eval_energy(self, sm, background=None, nodes=None):
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
        self.cg=sm.bg
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
            self.measures.append(score)
            return self.prefactor*score
        else:            
            self.measures.append(10**11)
            return 10**11

    def accept_last_measure(self):
        super(ProjectionMatchEnergy, self).accept_last_measure()
        self.accepted_projDir=self.projDir


class FPPEnergy(EnergyFunction):
    def __init__(self, pre, landmarks, scale, ref_image):
        super(FPPEnergy, self).__init__()
        self.prefactor = pre
        self.landmarks = landmarks
        self.scale = scale
        self.ref_image = fpp.to_grayscale(scipy.ndimage.imread(ref_image))

    @profile
    def eval_energy(self, sm, background=True, nodes=None):
        steplength = self.scale/self.ref_image.shape[0]
        ### Step 1: Calculate projection direction:
        vectors3d, angles, penalties = self.generate_equations(sm)
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
        proj, angleA, offset_centroidA, bs = self.find_offset(sm, projection_direction, mirror = True)
        img1, _ = proj.rasterize(self.ref_image.shape[0], bs, rotate = math.degrees(angleA), warn = False, virtual_residues = False)
        scoreA = fph.combined_distance(img1, self.ref_image)
        proj, angleB, offset_centroidB, bs = self.find_offset(sm, projection_direction, mirror = False)        
        img2, _ = proj.rasterize(self.ref_image.shape[0], bs, rotate = math.degrees(angleB), warn = False, virtual_residues = False)
        scoreB = fph.combined_distance(img2, self.ref_image)
        if scoreA<scoreB:
            sm.bg.project_from = ftuv.normalize(-projection_direction)
            score, img, params = fph.locally_minimal_distance(self.ref_image, self.scale, sm.bg, 
                                                          math.degrees(angleA), offset_centroidA, 
                                                          None, distance=fph.combined_distance,
                                                          maxiter=200, virtual_atoms="selected")
        else:
           score, img, params = fph.locally_minimal_distance(self.ref_image, self.scale, sm.bg, 
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
        return score
    @profile
    def find_offset(self, sm, projection_direction, mirror = False):
        if mirror: projection_direction = -projection_direction
        steplength = self.scale/self.ref_image.shape[0]
        sm.bg.project_from = ftuv.normalize(projection_direction)
        ### Step 2: Find out offset and rotation.
        proj = fpp.Projection2D(sm.bg, project_virtual_residues = [ x[0] for x in self.landmarks], project_virtual_atoms = "selected")
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
    def generate_equations(self, sm):
        penalty = 0
        #Preprocess: Calculate the projection angles for the shorthening of 
        # the 3 pairwise distances between landmarks.
        angles = []
        vectors3d = []
        for i, (l0, l1) in enumerate(it.combinations(self.landmarks, 2)):
            vec = sm.bg.get_virtual_residue(l1[0], True) - sm.bg.get_virtual_residue(l0[0], True)
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
class CombinedEnergy:
    def __init__(self, energies=None, uncalibrated_energies=None, normalize=False):
        """
        :param normalize: Divide the resulting energy by the numbers of contributions
        """
        if energies is not None:
            self.energies = energies
        else:
            self.energies=[]
        if uncalibrated_energies is not None:
            self.uncalibrated_energies = uncalibrated_energies
        else:
            self.uncalibrated_energies=[]
        self.bad_bulges = []
        self.accepted_measures=[float("nan")]
        self.normalize=normalize
    def uses_background(self):
        for e in it.chain(self.energies, self.uncalibrated_energies):
            if hasattr(e, "sampled_kdes"):
                return True
        return False
    def iterate_energies(self):
        """
        Iterate over all member enegies (calibrated and uncalibrated) 
        in a way that flattens nested CombinedEnergies.
        """
        for e in it.chain(self.energies, self.uncalibrated_energies):
            if isinstance(e, CombinedEnergy):
                for e2 in e.iterate_energies():
                    yield e2
            else:
                yield e
    def shortname(self):
        name=[]
        for e in it.chain(self.energies, self.uncalibrated_energies):
            sn=e.shortname()
            name.append(sn)
        sn=",".join(name)
        if len(sn)>12:
            sn=sn[:9]+"..."
        return sn

    def save_energy(self, energy, directory):
        if not os.path.exists(directory):
            os.makedirs(directory)

        filename = os.path.join(directory,
                                energy.__class__.__name__ + ".energy")
        #print ("saving filename:", filename)
        pickle.dump(energy, open(filename, 'w'))

    def resample_background_kde(self, struct):
        for e in it.chain(self.energies,
                          self.uncalibrated_energies):
            e.resample_background_kde(struct)
        pass
    def update_adjustment(self, step, cg):
        for e in it.chain(self.energies, self.uncalibrated_energies):
            e.update_adjustment(step, cg)
    def accept_last_measure(self):
        """CombinedEnergy.acceptLastMeasure"""
        for e in it.chain(self.energies, self.uncalibrated_energies):
            e.accept_last_measure()
    def reject_last_measure(self):
        """CombinedEnergy.rejectLastMeasure"""
        for e in it.chain(self.energies, self.uncalibrated_energies):
            e.reject_last_measure()

    def dump_measures(self, base_directory, iteration=None):
        for e in it.chain(self.energies,                       
                          self.uncalibrated_energies):
            e.dump_measures(base_directory, iteration)

    def eval_energy(self, sm, verbose=False, background=True,
                    nodes=None, use_accepted_measure=False):
        total_energy = 0.
        self.bad_bulges = []
        self.constituing_energies=[]
        num_contribs=0
        for energy in self.uncalibrated_energies:
            contrib=None
            if use_accepted_measure:
                try:
                    contrib = energy.eval_energy(sm, background=background,
                                                 nodes=nodes, use_accepted_measure=use_accepted_measure)
                except TypeError:
                    pass
            if contrib is None: #In case of TypeError or use_accepted_measure is False
                contrib = energy.eval_energy(sm, background=background,
                                             nodes=nodes)

            if not np.isscalar(contrib):
                contrib, = contrib
            self.constituing_energies.append((energy.shortname(), contrib))
            #try: print("uncalibrated ",energy.shortname(), " contributes ",contrib ,"to the combined energy")
            #except Exception: print("uncalibrated ",energy, " contributes ",contrib ,"to the combined energy")
            total_energy += contrib
            num_contribs +=1
            self.bad_bulges += energy.bad_bulges
            if verbose:
                print (energy.__class__.__name__, energy.shortname(), contrib)            
                if energy.bad_bulges:
                    print("bad_bulges:", energy.bad_bulges)

            log.debug("{} ({}) contributing {}".format(energy.__class__.__name__, energy.shortname(), contrib))
        for energy in self.energies:
            contrib=None
            if use_accepted_measure:
                try:
                    contrib = energy.eval_energy(sm, background=background,
                                                 nodes=nodes, use_accepted_measure=use_accepted_measure)
                except TypeError:
                    pass
            if contrib is None: #In case of TypeError or use_accepted_measure is False
                contrib = energy.eval_energy(sm, background=background,
                                             nodes=nodes)
            sn=energy.shortname()
            self.constituing_energies.append((sn, contrib))
            #try: print(energy.shortname(), " contributes ",contrib ,"to the combined energy")
            #except Exception: print(energy, " contributes ",contrib ,"to the combined energy")
            self.bad_bulges += energy.bad_bulges
            total_energy += contrib
            num_contribs +=1
            if verbose:
                print (energy.__class__.__name__, energy.shortname(), contrib)
                if energy.bad_bulges:
                    if isinstance(energy, StemVirtualResClashEnergy):
                        print(set(tuple(sorted([energy.bad_bulges[i], energy.bad_bulges[i+1]])) for i in range(0, len(energy.bad_bulges),2)))
                    else:
                        print("bad_bulges:", energy.bad_bulges)

            log.debug("{} ({}) contributing {}".format(energy.__class__.__name__, energy.shortname(), contrib))

        if self.normalize:
            total_energy=total_energy/num_contribs

        if verbose:
            print ("--------------------------")
            print ("total_energy:", total_energy)
        return total_energy

    def __str__(self):
        out_str = ''
        for en in it.chain(self.energies, self.uncalibrated_energies):
            out_str += en.__class__.__name__ + " "
        return out_str

class CoarseStemClashEnergy(EnergyFunction):
    '''
    Determine if two stems clash.
    '''

    def __init__(self):
        super(CoarseStemClashEnergy, self).__init__()
    def update_adjustment(self, step, cg): pass
    def eval_energy(self, sm, background=False, nodes=None):
        return 0.
        self.last_clashes=[]
        bg = sm.bg
        min_distance = 8.45
        energy = 0.
        #print list(bg.stems())

        if nodes is None:
            nodes = sm.bg.defines.keys()

        for (s1, s2) in it.combinations(bg.stems(), 2):
            if s1 not in nodes or s2 not in nodes:
                continue

            #print s1, s2
            if bg.are_any_adjacent_stems(s1, s2):
                continue

            closest_points = ftuv.line_segment_distance(bg.coords[s1][0],
                                                       bg.coords[s1][1],
                                                       bg.coords[s2][0],
                                                       bg.coords[s2][1])

            closest_distance = ftuv.vec_distance(closest_points[1], closest_points[0])
            #print "s1, s2", s1, s2, closest_distance

            if closest_distance < min_distance:
                self.last_clashes.append((s1,s2))
                energy += 100000.0
        return energy

class StemVirtualResClashEnergy(EnergyFunction):
    '''
    Determine if the virtual residues clash.
    '''

    def __init__(self):
        super(StemVirtualResClashEnergy, self).__init__()
        self.bad_bulges = []
    def update_adjustment(self, step, cg): pass
    def virtual_residue_atom_clashes_kd(self):
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
            kdt2.all_search(1.8) #Distance in Angstrom. #1.8

        clashes = 0
        indeces = kdt2.all_get_indices()
        for (ia,ib) in indeces:
            if virtual_atoms[ia][1][0] == virtual_atoms[ib][1][0]:
                continue
            if virtual_atoms[ia][1] == virtual_atoms[ib][1]:
                continue

            key1 = virtual_atoms[ia][1]
            key2 = virtual_atoms[ib][1]

            resn1 = self.bg.stem_side_vres_to_resn(key1[0], key1[2], key1[1])
            resn2 = self.bg.stem_side_vres_to_resn(key2[0], key2[2], key2[1])

            if abs(resn1 - resn2) == 1:
                continue
            self.bad_bulges += [key1[0], key2[0]]
            self.bad_atoms[key1[0]].append(virtual_atoms[ia][0])
            self.bad_atoms[key2[0]].append(virtual_atoms[ib][0])
            clashes += 1

        return clashes

    def virtual_residue_atom_clashes(self, bg, s1,i1,a1, s2, i2, a2):
        '''
        Check if any of the virtual residue atoms clash.
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

    def eval_energy(self, sm, background=False, nodes = None):
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
        self.bg = sm.bg
        self.bad_bulges = []
        self.bad_atoms = c.defaultdict(list)
        self.last_clashes=[]
        bg = sm.bg
        mult = 8
        points = []
        energy = 0.

        if nodes is None:
            nodes = sm.bg.defines.keys()

        for d in nodes:
            if d[0] == 's':
                s = d
                s_len = bg.stem_length(s)
                #stem_inv = bg.stem_invs[s]

                for i in range(s_len):
                    (p, v, v_l, v_r) = bg.v3dposs[d][i]

                    points += [(p+ mult * v_l, d, i, 1)]
                    points += [(p+ mult * v_r, d, i, 0)]

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

            if len(set.intersection(bg.edges[s1], bg.edges[s2])) > 0:
                # the stems are connected
                continue

            #potential_clashes += 1
            #fud.pv('s1,s2')

            if (s1,i1,a1) not in self.vras.keys():
                self.vras[(s1,i1,a1)] = ftug.virtual_residue_atoms(bg, s1, i1, a1)
            if (s2,i2,a2) not in self.vras.keys():
                self.vras[(s2,i2,a2)] = ftug.virtual_residue_atoms(bg, s2, i2, a2)

            #energy += 100000. * self.virtual_residue_atom_clashes(sm.bg, s1, i1, a1, s2, i2, a2)
        energy += 100000. * self.virtual_residue_atom_clashes_kd()
        return energy

class DistanceEnergy(EnergyFunction):

    def __init__(self, distance_constraints, multiplier= 10):
        super(DistanceEnergy, self).__init__()
        self.distance_constraints = distance_constraints
        self.multiplier = multiplier

    def eval_energy(self, sm, background=True):
        energy = 0.

        for constraint in self.distance_constraints:
            f = constraint[0]
            t = constraint[1]
            d = float(constraint[2])

            d1 = ftuv.vec_distance(sm.bg.get_point(f), sm.bg.get_point(t))

            energy += abs(d1 - d)

        return energy

class RoughJunctionClosureEnergy(EnergyFunction):
    def __init__(self):
        super(RoughJunctionClosureEnergy, self).__init__()
    def update_adjustment(self, step, cg): pass
    def eval_energy(self, sm, background=True, nodes=None):
        bg = sm.bg
        if nodes == None:
            nodes = bg.defines.keys()

        self.bad_bulges = []
        all_bulges = set([d for d in nodes if (d[0] == 'm' and d in nodes)])
        energy = 0.
        #closed_bulges = all_bulges.difference(sm.sampled_bulges)

        for bulge in all_bulges:
            #bl = bg.defines[bulge][1] - bg.defines[bulge][0] - 1
            bl = bg.get_bulge_dimensions(bulge)[0]
            #dist = ftug.junction_virtual_res_distance(bg, bulge)
            dist = ftug.junction_virtual_atom_distance(bg, bulge)            
            #
            #cutoff_distance = (bl) * 5.9 + 13.4
            #cutoff_distance = (bl) * 5.908 + 11.309
            #cutoff_distance = (bl) * 6.4 + 6.4
            cutoff_distance = (bl) * 6.22 + 14.0 #Peter's cyclic coordinate descent
            #cutoff_distance  = (bl) * 1.60 + 34.02 # 0.99-quantile-regression on nr 2.92 dataset.
            #cutoff_distance  = (bl) * 2.68 + 12.03  # 0.8-quantile-regression on nr 2.92 dataset.
            
            # Note: DOI: 10.1021/jp810014s claims that a typical MeO-P bond is 1.66A long. 
            
            if (dist > cutoff_distance):
                log.info("Junction closure: dist {} > cutoff {} for bulge {} with length {}".format(dist, cutoff_distance, bulge, bl))
                self.bad_bulges += bg.find_bulge_loop(bulge, 200) + [bulge]
                energy += (dist - cutoff_distance) * 10000.

        return energy

def merge(times):
    saved = list(times[0])
    for st, en in sorted([sorted(t) for t in times]):
        if st <= saved[1]:
            saved[1] = max(saved[1], en)
        else:
            yield tuple(saved)
            saved[0] = st
            saved[1] = en

    yield tuple(saved)

class DistanceExponentialEnergy(EnergyFunction):
    def __init__(self, from_elem, to_elem, distance, scale):
        '''
        Create an exponential distribution for the distance between two elements.
        '''
        super(DistanceExponentialEnergy, self).__init__()
        self.from_elem = from_elem
        self.to_elem = to_elem

        self.distance = distance
        self.scale = scale
        self.expon = ss.expon(loc=distance, scale=scale)
        self.bad_bulges = []

    def get_distance(self, sm):
        bg = sm.bg
        from_elem = self.from_elem
        to_elem = self.to_elem

        closest_points = ftuv.line_segment_distance(bg.coords[from_elem][0],
                                                   bg.coords[from_elem][1],
                                                   bg.coords[to_elem][0],
                                                   bg.coords[to_elem][1])

        closest_distance = ftuv.vec_distance(closest_points[1], closest_points[0])
        return closest_distance
    def shortname(self):
        if self.scale==1:
            return "CLA({},{})".format(self.from_elem, self.to_elem)
        else:
            return "{}CLA({},{})".format(self.scale,self.from_elem, self.to_elem) 
    def eval_energy(self, sm, background=True, nodes=None):
        closest_distance = self.get_distance(sm)

        if closest_distance < self.distance:
            return 0

        energy = -np.log(self.expon.pdf(closest_distance)) * 10.

        return energy

class CoarseGrainEnergy(EnergyFunction):
    def __init__(self, energy_prefactor=10):
        """
        A base class for Energy functions that use a background distribution.
        """
        super(CoarseGrainEnergy, self).__init__()

        self.real_kdes = dict()
        self.sampled_kdes = dict()

        self.flocs = []
        self.fscales = []
        self.vals = []
        self.dists = []

        self.dist_type = "kde"

        self.measures = []
        self.prev_energy = 0.

        self.energy_prefactor = energy_prefactor

        self.adjustment=1.0

        self.adjustment_stepwidth=0
        self.adjustment_steps_per_value=float("inf")
        self.current_adjustment_step=0

    def set_dynamic_adjustment(self, step, steps_per_value):
        self.adjustment_stepwidth=step
        self.adjustment_steps_per_value=steps_per_value
    def set_next_adjustment(self, cg):
        del self.real_kdes[self.measure_category(cg)]
        self.adjustment+=self.adjustment_stepwidth
        self.current_adjustment_step+=1
    def update_adjustment(self, step, cg):
        if step>self.adjustment_steps_per_value*(1+self.current_adjustment_step):
            self.set_next_adjustment(cg)
    def resample_background_kde(self, struct):
        values = self.accepted_measures

        if len(values) > 100:
            #print("RESAMPLING OF BACKGROUND")
            #print(self, values)
            new_kde = self.get_distribution_from_values(values)

            if new_kde is not None:
                self.sampled_kdes[self.measure_category(struct)] = new_kde
            else:
                #print ("skipping this time...", file=sys.stderr)
                pass

            #self.vals[1] = values

    def get_distribution_from_values(self, values):
        '''
        Return a probability distribution from the given values.

        @param values: The values to fit a distribution to.
        @return: A probability distribution fit to the values.
        '''
        floc = -0.1
        fscale =  1.5 * max(values)

        log.debug("Getting distribtion from values of shape {}".format(np.shape(values)))
        if self.dist_type == "kde":
            try:
                k = ss.gaussian_kde(values)
            except np.linalg.linalg.LinAlgError:
                return None
        else:
            f = ss.beta.fit(values, floc=floc, fscale=fscale)
            k = lambda x: ss.beta.pdf(x, f[0], f[1], f[2], f[3])

        '''
        self.flocs += [floc]
        self.fscales += [fscale]
        self.vals += [values]
        self.dists += [k]
        '''

        return k

    def measure_category(self, cg):
        '''
        Decide which target function we should use for this structure.

        @param cg: The CoarseGrain graph
        @return: Some sort of identifier to determine which energy distribution
                 to use.
        '''
        return cg.seq_length

    def eval_energy(self, sm, background=True, nodes=None, use_accepted_measure=False):
        '''
        A generic function which simply evaluates the energy based on the
        previously calculated probability distributions.
        '''
        cg = sm.bg

        if self.measure_category(cg) not in self.real_kdes.keys():
            try:
                (self.real_kdes[self.measure_category(cg)], x) = self.get_distribution_from_file(self.real_stats_fn, 
                                                                                             self.measure_category(cg), adjust=self.adjustment)
            except: 
                (self.real_kdes[self.measure_category(cg)], x) = self.get_distribution_from_file(self.real_stats_fn, 
                                                                                             self.measure_category(cg))
        if self.measure_category(cg) not in self.sampled_kdes.keys():
            (self.sampled_kdes[self.measure_category(cg)], self.measures) = self.get_distribution_from_file(self.sampled_stats_fn, 
                                                                                             self.measure_category(cg))
            self.accepted_measures = self.measures[:]

        kr = self.real_kdes[self.measure_category(cg)]
        ks = self.sampled_kdes[self.measure_category(cg)]
        if use_accepted_measure:
            m = self.accepted_measures[-1]
        else:
            m = self.get_cg_measure(sm)

        if False: #For debuging
            import matplotlib.pyplot as plt
            xs=np.linspace(0, 200, 2000)
            fig,ax1 = plt.subplots()
            ax2 = ax1.twinx()
            ax1.plot(xs, ks(xs), label="sampled")
            ax1.plot(xs, kr(xs), label="reference")
            ax2.plot(xs, -(np.log(kr(xs) + 0.00000001 * ks(xs)) - np.log(ks(xs))), label="energy", color="red")
            plt.title(self.shortname())
            ax1.legend(loc="lower left")
            ax2.legend()
            plt.show()

        log.debug("Measure is {:1.4f}".format(m))
        
        self.measures.append(m)
        
        if background:
            energy, = (np.log(kr(m) + 0.00000001 * ks(m)) - np.log(ks(m)))
            self.prev_energy = energy
            self.prev_cg = m
            #if isinstance(self, ShortestLoopDistancePerLoop):
            #    print("Measure {}, energy {}.".format(m, (-1 * self.energy_prefactor * energy)[0]))
            #energy = (np.log(kr.integrate_box_1d(0., m) + 0.0001 * ks.integrate_box_1d(0., m)) - np.log(ks.integrate_box_1d(0., m)))
            return -1 * self.energy_prefactor * energy
        else:
            l = np.log(kr(m))
            log.debug("Energy, = {}".format(l))
            energy, = l
            return -energy


class HausdorffEnergy(CoarseGrainEnergy):
    def __init__(self, img, scale, prefactor):
        """
        :param img: A boolean, square 2D numpy.array
        :param scale: Int or float. How many Angstrom the side length of the image is.
        :param prefactor: Multiply the pseudo-energy with this value.
        """
        super(HausdorffEnergy, self).__init__()        
        self.prefactor=prefactor
        self.ref_img=img

        self.sampled_stats_fn = 'stats/hausdorff_sampled.csv'
        self.sampled_stats_fn = op.expanduser(self.sampled_stats_fn)

        self.real_stats_fn = 'stats/hausdorff.csv'
        self.real_stats_fn = op.expanduser(self.real_stats_fn)
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                self.ref_quarter=(scipy.ndimage.zoom(img, 0.3)>150)
        except RuntimeError:
            print("Image {} could not be read properly ".format(img), file=sys.stderr)
            print("Make sure, it has the right type (Ideally 8 bit WITHOUT alpha)", file=sys.stderr)
            raise
        self.get_refimg_longest_axis()
        self.scale=scale
        self.last_dir=None
        self.accepted_projDir=None
        self.last_img=None
        self.accepted_img=None

    def get_distribution_from_file(self, filename, length, adjust=1.):
        data = np.genfromtxt(load_local_data(filename), delimiter=' ')
        data=np.array([adjust*i for i in data])
        return (self.get_distribution_from_values(data), list(data))
        
    def get_refimg_longest_axis(self):
        max_sdist=0
        best_points=None
        for p1, p2 in it.combinations(np.transpose(np.where(self.ref_img)),2):
            if (p2[0]-p1[0])**2+(p2[1]-p1[1])**2>max_sdist:
                max_sdist=(p2[0]-p1[0])**2+(p2[1]-p1[1])**2
                best_points=(p1,p2)
        deg=math.degrees(math.atan2(best_points[1][0]-best_points[0][0], 
                                    best_points[1][1]-best_points[0][1]))
        deg=90-deg
        self.initial_degrees=[deg%360, (deg+180)%360]
    def shortname(self):
        return "HDE"
    def get_name(self):
        return "Hausdorff-Energy"
    def get_cg_measure(self, sm):
        for residuePos in range(1,sm.bg.total_length()):
            residue=sm.bg.virtual_atoms(residuePos) #To remove this from the profiling of rasterization
        st=time.time()
        if False: #len(np.where(self.ref_quarter)[0])>100: #Only use ref_quarter, if it has enough white.
            print("Using refQuarter")
            s, i, params = fph.globally_minimal_distance(self.ref_quarter, self.scale, sm.bg,  
                                                          start_points=20, 
                                                          starting_rotations=self.initial_degrees,
                                                          virtual_atoms=False)#, distance=fph.modified_hausdorff_distance)         
            score, img, params = fph.locally_minimal_distance(self.ref_img, self.scale, sm.bg, 
                                                          params[1], params[2], params[0], 
                                                          maxiter=200)#, distance=fph.modified_hausdorff_distance)
        else:
            #print("Using Full Resolution")
            s, i, params = fph.globally_minimal_distance(self.ref_img, self.scale, sm.bg,  
                                                          start_points=20,
                                                          starting_rotations=self.initial_degrees,
                                                          virtual_atoms="selected", distance=fph.tp_fp_distance, use_heuristic=False)         
            score, img, params = fph.locally_minimal_distance(self.ref_img, self.scale, sm.bg, 
                                                          params[1], params[2], params[0], 
                                                          maxiter=200, distance=fph.tp_fp_distance)
        self.last_dir=params[0]
        self.measures.append(score)
        #import matplotlib.pyplot as plt
        #fig, ax=plt.subplots(2)
        #ax[0].imshow(self.ref_img, interpolation="none", cmap='gray')
        #ax[1].imshow(img, interpolation="none", cmap='gray')
        #ax[0].set_title("Final Opt Reference")
        #ax[1].set_title("{} distance".format(score))
        #plt.show()
        #print("Time was {:2.2f}".format(time.time()-st))
        return score

    def accept_last_measure(self):
        super(HausdorffEnergy, self).accept_last_measure()
        self.accepted_projDir=fph.from_polar([1,self.last_dir[0], self.last_dir[1]])
        self.accepted_img = self.last_img
        
    def dump_measures(self, base_directory, iteration=None):
        '''
        Dump all of the coarse grain measures collected so far
        to a file.
        '''
        super(HausdorffEnergy, self).dump_measures(base_directory, iteration)
        image_file = op.join(base_directory, "HDE_image."+iteration+".png")
        scipy.misc.imsave(image_file, self.accepted_img)



class CylinderIntersectionEnergy(CoarseGrainEnergy):
    def __init__(self):
        super(CylinderIntersectionEnergy, self).__init__()
        self.max_dist = 30.
        self.min_ratio = 0.
        self.max_ratio = 30.

        self.real_stats_fn = 'stats/cylinder_intersections_native.csv'
        #self.sampled_stats_fn = 'stats/cylinder_intersections_1jj2_rog.csv'
        self.sampled_stats_fn = 'stats/cylinder_intersections_loop_rog.csv'

        self.real_data = None
        self.fake_data = None

        self.real_kdes = dict()
        self.sampled_kdes = dict()

    def get_name(self):
        return "Cylinder Overlap"

    def load_data(self, filename):
        import pandas as pa

        t = pa.read_csv(load_local_data(filename), header=None, sep=' ')

        ratios = t[t.columns[1]]
        ratios = list(ratios[~np.isnan(ratios)])

        #new_ratios = [rand.choice(ratios) for i in range(5000)]
        #new_ratios = ratios

        #new_ratios += list(np.linspace(self.min_ratio, self.max_ratio, 100))
        new_ratios = ratios

        return (new_ratios, stats.gaussian_kde(new_ratios))

    def get_distribution_from_file(self, filename, length):
        '''
        Return the probability distribution of the sum of the cylinder
        intersection lengths.

        @param filename: The filename from which to load the distribution.
        @param length: The length of the RNA molecule
        @return: A probability distribution describing the combined cylinder
                 intersection length.
        '''
        lengths = self.get_measures_from_file(filename, length)

        return (self.get_distribution_from_values(lengths), lengths)

    def get_measures_from_file(self, filename, length):
        all_lengths = self.get_cylinder_intersections_from_file(filename, length)
        lengths = map(sum, all_lengths)

        return lengths

    def get_cylinder_intersections_from_file(self, filename, length):
        '''
        Get the list of cylinder intersection lengths from a file for a
        molecule of a particular length.

        Included in the distribution will be all molecules where

        0.8 * length < len(molecule) < 1.2 * length

        The cylinder intersection lengths for all of these molecules will
        be returned as a list of lists where each sublist contains the
        cylinder intersection lengths for one particular file.

        @param filename: The name of the file containing the cylinder
                         intersection lengths
        @param length: The length of the molecule
        @return: A list of lists containing the cylinder intersection lengths
        '''
        cls = c.defaultdict(list)

        f = load_local_data(filename)
        for line in f:
            parts = line.strip().split()
            mol_size = int(parts[0])
            lengths = map(float, parts[1:])
            cls[mol_size] += [lengths]

        mol_sizes = cls.keys()

        mol_sizes = np.array(mol_sizes)
        mol_sizes = mol_sizes[mol_sizes > length * distribution_lower_bound]
        mol_sizes = mol_sizes[mol_sizes < length * distribution_upper_bound]

        all_lengths = []
        for l in mol_sizes:
            all_lengths += cls[l]

        #all_lengths = [cls[l] for l in mol_sizes]

        return all_lengths

    def calculate_intersection_coverages(self, bg):
        in_cyl_fractions = c.defaultdict(lambda: 0.001)
        covered = c.defaultdict(list)
        cyls = c.defaultdict(lambda: [np.array([0.,0.,0.]), np.array([0.,0.,0.])])
        stem_iloops = list(bg.stem_iterator()) + list(bg.iloop_iterator())

        if len(stem_iloops) == 1:
            return {stem_iloops[0]:0.}

        for (s1, s2) in it.permutations(list(bg.stem_iterator()) + list(bg.iloop_iterator()), 2):
            line = bg.coords[s1]
            cyl = bg.coords[s2]
            extension = 0.

            # the distance between the closest points on the line and cylinder
            (i1, i2) = ftuv.line_segment_distance(bg.coords[s1][0],
                                                bg.coords[s1][1],
                                                bg.coords[s2][0],
                                                bg.coords[s2][1])

            # make sure they're not connected or too far
            dist = ftuv.vec_distance(i1, i2)
            if dist > 30. or dist < 0.01:
                continue

            # extend the cylinder on either side
            cyl_vec = ftuv.normalize(bg.coords[s2][1] - bg.coords[s2][0])
            cyl = [cyl[0] - extension * cyl_vec,
                   cyl[1] + extension * cyl_vec]
            cyls[s2] = cyl

            intersects = ftuv.cylinder_line_intersection(cyl, line,
                                                        self.max_dist)

            if len(intersects) > 0 and np.isnan(intersects[0][0]):
                sys.exit(1)

            if len(intersects) == 0:
                pass
            else:
                #in_cyl_len = ftuv.magnitude(intersects[1] - intersects[0])
                c1 = ftuv.closest_point_on_seg(cyl[0], cyl[1], intersects[0])
                c2 = ftuv.closest_point_on_seg(cyl[0], cyl[1], intersects[1])

                poss = [ftuv.vec_distance(cyl[0], c1), ftuv.vec_distance(cyl[0], c2)]
                poss.sort()

                covered[s2] += [poss]
                '''
                cyl_basis = ftuv.create_orthonormal_basis(cyl_vec)
                intersects_t = ftuv.change_basis((intersects - cyl[0]).T,
                                                cyl_basis,
                                                ftuv.standard_basis).T
                in_cyl_len = abs(intersects_t[1][0] - intersects_t[0][0])
                '''
                #covered[s1] += [(intersects_t[0][0], intersects_t[1][0])]

            #in_cyl_fractions[s1] += in_cyl_len / line_len

        for s in list(bg.stem_iterator()) + list(bg.iloop_iterator()):
            if len(covered[s]) == 0:
                continue

            ms = list(merge(covered[s]))
            cl = sum(map(lambda x: x[1] - x[0], ms))

            #in_cyl_fractions[s] = cl / total_len
            in_cyl_fractions[s] = cl

        return in_cyl_fractions

    def get_cg_measure(self, sm):
        '''
        Calculate the coarse grain measure that is being described by this
        energy function.

        @param sm: The SpatialModel for which to calculate this measure.
        @return: A single floating point number describing this measure.
        '''
        cyl_fractions = self.calculate_intersection_coverages(sm.bg)

        total_cylinder_intersections = sum(cyl_fractions.values())

        return total_cylinder_intersections

class CheatingEnergy(EnergyFunction):
    def __init__(self, real_bg):
        super(CheatingEnergy, self).__init__()
        self.real_bg = copy.deepcopy(real_bg)
        self.real_residues = ftug.bg_virtual_residues(self.real_bg)

    def eval_energy(self, sm, background=True, nodes=None):
        '''
        @param sm: A SpatialModel, which contains a coarse grain model (sm.bg)
        '''
        new_residues = ftug.bg_virtual_residues(sm.bg)

        return  ftms.rmsd(self.real_residues, new_residues)*30
    def shortname(self):
        return "CHE"

def length_and_rog(cg):
    coords = cg.get_ordered_stem_poss()
    rog = ftmd.radius_of_gyration(coords)
    total_length = sum([len(list(cg.define_residue_num_iterator(d))) for d in cg.defines.keys()])
    
    return (total_length, rog)

def length_and_rog_from_file(filename):
    cg = ftmc.CoarseGrainRNA(op.expanduser(filename))
    
    return length_and_rog(cg)

class SimpleRadiusOfGyrationEnergy(EnergyFunction):
    def __init__(self):
        super(SimpleRadiusOfGyrationEnergy, self).__init__()

    def eval_energy(self, sm, background=True, nodes=None):
        cg = sm.bg
        (length, rog) = length_and_rog(cg)

        return -rog
    
class RadiusOfGyrationEnergy(CoarseGrainEnergy):
    def __init__(self, dist_type="kde", adjustment=1., energy_prefactor=DEFAULT_ENERGY_PREFACTOR):
        super(RadiusOfGyrationEnergy, self).__init__(energy_prefactor=energy_prefactor)
        self.sampled_stats_fn = 'stats/subgraph_radius_of_gyration_sampled.csv'
        self.sampled_stats_fn = op.expanduser(self.sampled_stats_fn)

        self.real_stats_fn = 'stats/subgraph_radius_of_gyration.csv'
        self.real_stats_fn = op.expanduser(self.real_stats_fn)

        self.real_rogs = dict()
        self.sampled_rogs = dict()
        self.real_kdes = dict()
        self.sampled_kdes = dict()

        self.background=True
        self.dist_type = 'kde'
        #: the adjustment is used to enlarge or shrink the target distribution
        self.adjustment = adjustment 

    def shortname(self):
        if self.energy_prefactor==DEFAULT_ENERGY_PREFACTOR:
            pre=""
        else:
            pre=self.energy_prefactor
        if self.adjustment==1.0:
            adj=""
        else:
            adj=self.adjustment
        return "{}ROG{}".format(pre,adj)

    def get_name(self):
        return "Radius Of Gyration"

    def get_cg_measure(self, sm):
        (length, rog) = length_and_rog(sm.bg)
        self.measures += [rog]
        log.debug("ROG-energy. get_measure returning a rog of {}".format(rog))
        return rog

    def get_distribution_from_file(self, filename, length, adjust=1.):
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
        rogs=np.array([adjust*i for i in rogs])

        log.info("Radii of gyration [{}, {}]: min-median-max: {}-{}-{}".format(filename, adjust, min(rogs), sorted(rogs)[int(len(rogs)/2)], max(rogs)))
        return (self.get_distribution_from_values(rogs), list(rogs))

class NormalDistributedRogEnergy(RadiusOfGyrationEnergy):
    def __init__(self, dist_type="kde", adjustment=1., energy_prefactor=DEFAULT_ENERGY_PREFACTOR):
        """
        A Subclass of the Radius of Gyration energy with a normal distributed target distribution.

        The adjustment is treated as the estimated perimeter of the RNA. Thus the target distribution
        is a normal distribution with mean=0.77*ADJ and STDDEV=0.23*ADJ.
        
        Remember that the ROG for polymers is defined as:
        ROG=SQRT(1/N*Sum(r)) for N points with the mean in the origin.
  
        0.77 is roughly SQRT(3/5), which is the limes of the rog for m->inf many points
        in a 3D unit sphere with the nth point placed at radius (n/m)**1/3
        0.77-0.23 is roughly 1/SQRT(3), which is the limes of the rog for m->inf many points
        equally distributed along one or more diameters (in a star-like fashion)
        0.77+0.23 is 1, corresponding to all points on the surface of the sphere.

        If the perimeter is estimated from electron microscopy data, 
        it is potentially smaller than the perimeter in 3D space if the shape of the RNA
        is not a perfect sphere.
        This effect is accounted for by allowing values greater than 1*ADJ with a 
        non-neglectable probability
        """
        super(NormalDistributedRogEnergy, self).__init__(dist_type, adjustment, energy_prefactor)
        self.real_stats_fn = None
    def get_distribution_from_file(self, filename, length, adjust=1.):
        if filename is None:
            return lambda x: np.array([stats.norm(loc=0.77*adjust, scale=0.23*adjust).pdf(x)]), None
        else:
            return super(NormalDistributedRogEnergy, self).get_distribution_from_file(filename, length, adjust)
    def shortname(self):
        if self.energy_prefactor==DEFAULT_ENERGY_PREFACTOR:
            pre=""
        else:
            pre=self.energy_prefactor
        if self.adjustment==1.0:
            adj=""
        else:
            adj=self.adjustment
        return "{}NDR{}".format(pre,adj)

class ShortestLoopDistanceEnergy(RadiusOfGyrationEnergy):
    def __init__(self):
        super(ShortestLoopDistanceEnergy, self).__init__()
        self.max_dist = 60
    
        self.real_stats_fn = 'stats/loop_loop2_distances_native.csv'
        self.sampled_stats_fn = 'stats/loop_loop2_distances_sampled.csv'

    def get_name(self):
        return "Loop Distance"

    def get_distribution_from_file(self, filename, length):
        data = np.genfromtxt(load_local_data(filename), delimiter=' ')

        rdata = data[data[:,0] == length]

        rogs = rdata[:,1]
        return (self.get_distribution_from_values(rogs), list(rogs))

    def measure_category(self,cg):
        return len(list(cg.hloop_iterator()))

    def get_shortest_distance_list(self, cg):
        pairs = []

        for (l1, l2) in it.combinations(cg.hloop_iterator(), 2):
            (i1, i2) = ftuv.line_segment_distance(cg.coords[l1][0],
                                                   cg.coords[l1][1],
                                                   cg.coords[l2][0],
                                                   cg.coords[l2][1])

            pairs += [(l1, l2, ftuv.vec_distance(i2, i1))]

        evaluated = set()
        to_eval = []

        pairs.sort(key=lambda x: x[2])
        for p in pairs:
            if p[0] not in evaluated and p[1] not in evaluated:
                to_eval += [p]

                evaluated.add(p[0])
                evaluated.add(p[1])

        all_dists = []
        for (l1, l2, dist) in to_eval:
            
            if dist > self.max_dist:
                continue

            all_dists += [dist]

        return all_dists

    def get_shortest_distances(self, cg):
        all_dists = self.get_shortest_distance_list(cg)
        return sum(all_dists)

    def get_cg_measure(self, sm):
        #import traceback

        #traceback.print_stack()
        return self.get_shortest_distances(sm.bg)

class DoNotContribute(Exception):
  pass

class ShortestLoopDistancePerLoop(ShortestLoopDistanceEnergy):
    def __init__(self, loop_name, prefactor=10):
        super(ShortestLoopDistancePerLoop, self).__init__()
        self.loop_name = loop_name
        self.energy_prefactor=prefactor
        self.real_stats_fn = 'stats/loop_loop3_distances_native.csv'
        self.sampled_stats_fn = 'stats/loop_loop3_distances_sampled.csv'

    def shortname(self):
        if self.energy_prefactor==DEFAULT_ENERGY_PREFACTOR:
            pre=""
        else:
            pre=self.energy_prefactor
        if self.adjustment==1.0:
            adj=""
        else:
            adj=self.adjustment
        return "{}SLD{}".format(pre,adj)

    def measure_category(self, cg):
        return 1

    def get_distribution_from_file(self, filename, length):
        data = np.genfromtxt(load_local_data(filename), delimiter=' ')
        if filename == self.sampled_stats_fn:
            rogs=np.linspace(3, 300, num=30)
            data=np.concatenate([data, data,data,rogs])
        else:
            rogs=np.linspace(3, 300, num=15)
            data=np.concatenate([data, data,data,rogs])
        rogs = data
        return (self.get_distribution_from_values(rogs), list(rogs))

    def get_cg_measure(self, sm):
        min_dist = float("inf")

        cg = sm.bg
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

    def get_energy_name(self):
        '''
        Return the name of the energy.
        '''

        return self.__class__.__name__.lower() + "_" + self.loop_name + ".measures"

    def eval_energy(self, sm, background=True, nodes=None):
        '''
        We want to return an energy of 0. if there's less than two hairpin
        loops.
        '''
        if len(list(sm.bg.hloop_iterator())) < 2:
            return 0.
        else:
            try:
              energy= super(ShortestLoopDistancePerLoop, self).eval_energy(sm,
                                                                         background,
                                                                         nodes)
            except DoNotContribute:
              energy=0
            return energy



        
class AMinorEnergy(CoarseGrainEnergy):
    def __init__(self, dist_type="kde", adjustment=1., energy_prefactor=DEFAULT_ENERGY_PREFACTOR, loop_type='h'):
        import pandas as pd
        super(AMinorEnergy, self).__init__(energy_prefactor=energy_prefactor)
        self.adjustment=adjustment
        self.real_stats_fn = 'stats/aminors_1s72.csv'
        self.sampled_stats_fn = 'stats/aminors_1jj2_sampled.csv'

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
        self.types = {'h':0, 'i':1, 'm':2, 's': 3, 'f': 4, 't': 5}
        self.loop_type = self.types[loop_type]

        self.prob_funcs = dict()
        p_d_a_a2_given_i = dict()
        p_d_a_a2 = dict()
        p_i = dict()

        for lt in ['i', 'h']:
            dl = self.dall[lt] = dall[[lt in x for x in dall['l1']]]
            db = self.dbg_close[lt] = dbg_close[[lt in x for x in dbg_close['l1']]]
            p_i[lt] = len(dl['dist']) / float(len(db['dist']))

            p_d_a_a2_given_i[lt] = ss.gaussian_kde(np.array(zip(dl["dist"], dl["angle"], dl["angle2"])).T)
            p_d_a_a2[lt] = ss.gaussian_kde(np.array(zip(db["dist"], db["angle"], db["angle2"])).T)
            #self.prob_funcs[lt] = lambda point: (p_d_a_a2_given_i[lt](point) * p_i[lt]) / (p_d_a_a2[lt](point))
            self.prob_funcs[lt] = lambda point: (p_d_a_a2_given_i[lt](point)) / (p_d_a_a2[lt](point) + p_d_a_a2_given_i[lt](point))
        #print ("SELF>MEASURES ({}): {}".format(self.shortname(),len(self.measures)))
        #print (self.measures)

        #: The number of coarse grain elements considered in this energy
        self.num_loops=None

    def shortname(self):
        if self.energy_prefactor==DEFAULT_ENERGY_PREFACTOR:
            pre=""
        else:
            pre=self.energy_prefactor
        if self.adjustment==1.0:
            adj=""
        else:
            adj=self.adjustment
        return "{}AME({}){}".format(pre,self.loop_type, adj)

    def get_distribution_from_file(self, filename, length, adjust=1.0):
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
        rogs= np.array([adjust*i for i in rogs])
        return (self.get_distribution_from_values(rogs), list(rogs))

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
                                       cg.coords[s][1], 30):
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

    def get_cg_measure(self, sm):
        for d in sm.bg.defines.keys():

            # the loop type is encoded as an integer so that the stats file can be 
            # loaded using numpy
            if self.types[d[0]] != self.loop_type or 'A' not in "".join(sm.bg.get_define_seq_str(d)):
                continue

            m = self.eval_prob(sm.bg, d)[0]

            return m

    def get_name(self):
        if self.loop_type == self.types['i']:
            return "A-Minor Energy (interior loops)"
        elif self.loop_type == self.types['h']:
            return "A-Minor Energy (hairpin loops)"
        else:
            return "A-Minor Energy"

    def accept_last_measure(self):
        """AMinor.acceptLastMeasure"""
        #The AMinor energy appends more than one measure per step
        #print (hex(id(self)), "Before Accepting measures. length {} Now ".format(self.num_loops), self.accepted_measures)
        if len(self.measures) >0 and self.num_loops>0:        
            self.accepted_measures += self.measures[-self.num_loops:]
        #print (hex(id(self)), "Accepting measures. Now ", self.accepted_measures)

    def reject_last_measure(self):
        """AMinor.rejectLastMeasure"""
        if len(self.accepted_measures) > 0 and self.num_loops>0:
            self.accepted_measures += self.accepted_measures[-self.num_loops:]
    def get_num_loops(self, cg):
        return len([d for d in cg.defines.keys() if self.types[d[0]] == self.loop_type and 'A' in "".join(cg.get_define_seq_str(d))])
    def eval_energy(self, sm, background=True, nodes=None): #@PROFILE: This takes >50% of the runtime with default energy
        cg = sm.bg
        if self.measure_category(cg) not in self.real_kdes.keys():
            (self.real_kdes[self.measure_category(cg)], self.real_measures) = self.get_distribution_from_file(self.real_stats_fn, self.measure_category(cg), adjust=self.adjustment)
        if self.measure_category(cg) not in self.sampled_kdes.keys():
            (self.sampled_kdes[self.measure_category(cg)], self.measures) = self.get_distribution_from_file(self.sampled_stats_fn, self.measure_category(cg))
            self.accepted_measures = self.measures[:]
        kr = self.real_kdes[self.measure_category(cg)]
        ks = self.sampled_kdes[self.measure_category(cg)]

        energy = 0

        if self.num_loops is None:
            self.num_loops=self.get_num_loops(sm.bg)
        for d in sm.bg.defines.keys():

            # the loop type is encoded as an integer so that the stats file can be 
            # loaded using numpy
            if self.types[d[0]] != self.loop_type or 'A' not in "".join(sm.bg.get_define_seq_str(d)):
                continue

            m = self.eval_prob(sm.bg, d)[0]
            self.measures.append(m)

            if background:
                prev_energy = (np.log(kr(m) + 0.00000001 * ks(m)) - np.log(ks(m)))
                self.prev_energy = energy
                self.prev_cg = m
                energy +=  -1 * self.energy_prefactor * prev_energy                
            else:
                energy +=  -np.log(kr(m))
        if False: #For debuging
                import matplotlib.pyplot as plt
                xs=np.linspace(0, 1, 500)
                fig,ax1 = plt.subplots()
                ax2 = ax1.twinx()
                ax1.plot(xs, ks(xs), label="sampled")
                ax1.plot(xs, kr(xs), label="reference")
                ax2.plot(xs, -(np.log(kr(xs) + 0.00000001 * ks(xs)) - np.log(ks(xs))), label="energy", color="red")
                ax1.plot(self.measures, [1]*len(self.measures), "o", label="Measures")
                plt.title(self.shortname())
                ax1.legend(loc="lower left")
                ax2.legend()
                plt.show()
        return energy[0]

"""
class SpecificAMinorEnergy(AMinorEnergy):
    def __init__(self, dist_type="kde", adjustment=1., energy_prefactor=DEFAULT_ENERGY_PREFACTOR, loop_name='h0'):
        super(SpecificAMinorEnergy, self).__init__(energy_prefactor=energy_prefactor)

        self.real_stats_fn = 'real'
        self.sampled_stats_fn = 'sampled'
        self.loop_name = loop_name
        self.adjustment=adjustment
    def get_distribution_from_file(self, filename, category, adjust=1.0):
        if filename == 'real':
            rogs = list(np.linspace(0, 1, 100)) + [1.1] * 400
        else:
            rogs = np.linspace(0, 1, 100)
 
        return (self.get_distribution_from_values(rogs), list(rogs))

    def get_cg_measure(self, cg):
        d = self.loop_name

        return self.eval_prob(cg, d)

    def eval_energy(self, sm, background=True, nodes=None, new_nodes=None):
        '''
        A generic function which simply evaluates the energy based on the
        previously calculated probability distributions.
        '''
        cg = sm.bg

        if self.measure_category(cg) not in self.real_kdes.keys():
            (self.real_kdes[self.measure_category(cg)], self.real_measures) = self.get_distribution_from_file(self.real_stats_fn, self.measure_category(cg))
            (self.sampled_kdes[self.measure_category(cg)], self.measures) = self.get_distribution_from_file(self.sampled_stats_fn, self.measure_category(cg))

            self.accepted_measures = self.measures[:]


        kr = self.real_kdes[self.measure_category(cg)]
        ks = self.sampled_kdes[self.measure_category(cg)]

        d = self.loop_name

        m = self.eval_prob(sm.bg, d)[0]
        self.measures.append(m)

        energy = 0
        if background:
            prev_energy = (np.log(kr(m) + 0.00000001 * ks(m)) - np.log(ks(m)))
            self.prev_energy = energy
            self.prev_cg = m
            #energy = (np.log(kr.integrate_box_1d(0., m) + 0.0001 * ks.integrate_box_1d(0., m)) - np.log(ks.integrate_box_1d(0., m)))
            energy +=  -1 * self.energy_prefactor * prev_energy
        else:
            energy +=  -np.log(kr(m))
        
        return energy
"""
