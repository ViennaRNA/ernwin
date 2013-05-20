#!/usr/bin/python

import time
import pdb
import pickle, os
import Bio.PDB as bpdb
import copy
import itertools as it
import math
import warnings

import scipy.ndimage as sn
import scipy.spatial as ss
import Bio.KDTree as kd
import numpy as np
import numpy.linalg as nl
import random as rand

import collections as c

import corgy.utilities.vector as cuv
import corgy.graph.bulge_graph as cgb
import corgy.graph.graph_pdb as cgg
import corgy.exp.kde as cek
import corgy.builder.models as cbm
import corgy.builder.reconstructor as rtor
import corgy.utilities.statistics as cus
import corgy.builder.sampling as cbs
import corgy.builder.config as cbc
import corgy.utilities.debug as cud
import corgy.builder.rmsd as cbr

import scipy.stats as stats
#import scipy.stats as ss
import sys

def my_log(x):
    return np.log(x + 1e-200)

class MissingTargetException(Exception):
    def __init__(self, message):
        Exception.__init__(self, message)

class EnergyFunction(object):
    '''
    The base class for energy functions.
    '''

    def __init__(self):
        self.interaction_energies = None

    def eval_energy(self, sm, background=True):
        '''
        The base energy function simply returns a random number.
        '''
        return rand.random()

    def iterate_over_interaction_energies(self, bg, background):
        sm = cbm.SpatialModel(bg)

        for stem in bg.stems():
            cgg.add_virtual_residues(bg, stem)


        self.eval_energy(sm, background)
        for key in self.interaction_energies.keys():
            if abs(self.interaction_energies[key]) > 0.0000001:
                yield (key, self.interaction_energies[key])

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
        pass

    def calibrate(self, sm, iterations = 10, bg_energy = None):
        '''
        Calibrate this energy function.

        This is done by sampling structures given a background energy function.
        The sampled structures are used to normalize the energy of this
        function.
        '''
        self.calc_fg_parameters(sm.bg)

        stats = cbs.SamplingStatistics(sm)
        stats.silent = True

        # if not background energy function is provided, then the background
        # distribution is simply the proposal distribution implicit in 
        # cbs.GibbsBGSampler

        if bg_energy == None:
            bg_energy = EnergyFunction()
        
        gs = cbs.GibbsBGSampler(copy.deepcopy(sm), bg_energy, stats)
        for _ in range(iterations):
            gs.step()

        # Get the set of sampled structures
        ser_structs = sorted(stats.energy_rmsd_structs, key=lambda x: x[0])

        # I only want the top 2/3 of the sampled structures
        selected = ser_structs[:2 * (len(ser_structs) / 3)]
        selected_structs = [s[2] for s in selected]

        self.calc_bg_parameters(selected_structs)
             
        
class RandomEnergy(EnergyFunction):
    '''
    An energy function that always just returns a random value.
    '''
    def __init__(self):
        super(RandomEnergy, self).__init__()

    def eval_energy(self, sm, background=True):
        return rand.uniform(-5, 3)

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
            for j in range(i+1, len(keys)):
                d1 = keys[i]
                d2 = keys[j]

                point1 = bg.get_point(d1)
                point2 = bg.get_point(d2)

                dist = cuv.vec_distance(point1, point2)

                #if dist > 6.0 and dist < 25.0:
                if dist > self.min_distance and dist < self.max_distance:
                    yield tuple(sorted([d1, d2]))

lri_iter = DistanceIterator(6., 25.)

class CombinedEnergy:
    def __init__(self, energies=[], uncalibrated_energies=[]):
        self.energies = energies
        self.uncalibrated_energies = uncalibrated_energies

    def save_energy(self, energy, directory):
        if not os.path.exists(directory):
            os.makedirs(directory)

        filename = os.path.join(directory, energy.__class__.__name__ + ".energy")
        print "saving filename:", filename
        pickle.dump(energy, open(filename, 'w'))

    def calibrate(self, sm, iterations=40, bg_energy=None, output_dir='/home/mescalin/pkerp/projects/ernwin/energies'):
        '''
        Calibrate each of the energies by taking into account the
        background distribution induced by non-energy directed 
        sampling.
        '''
        self.energies[0].calibrate(sm, iterations)
        filename = os.path.join(output_dir, str(sm.bg.name))
        filename = os.path.join(filename, str(iterations))

        self.save_energy(self.energies[0], filename)
        filename = os.path.join(filename, self.energies[0].__class__.__name__)

        for i in range(1, len(self.energies)):
            ce = CombinedEnergy(self.energies[:i])
            self.energies[i].calibrate(sm, iterations, ce)

            self.save_energy(self.energies[i], filename)
            filename = os.path.join(filename, self.energies[i].__class__.__name__)

        self.save_energy(self, filename)

    def eval_energy(self, sm, verbose=False, background=True):
        total_energy = 0.

        for energy in self.uncalibrated_energies:
            total_energy += energy.eval_energy(sm)
    
        for energy in self.energies:
            contrib = energy.eval_energy(sm, background)

            total_energy += contrib

            if verbose:
                print energy.__class__.__name__, contrib

        return total_energy

class SkewNormalInteractionEnergy(EnergyFunction):
    '''
    This energy will assume that all elements need to be a certain
    distance from each other to interact.

    This distance is centered at around 15 angstroms.

    The distribution of distances is modeled by a skew-normal-distribution and
    the total energy will be the sum of log probabilities for the interactions.
    '''
    def __init__(self):
        super(SkewNormalInteractionEnergy, self).__init__()
        self.fg = None
        self.bgs = dict()

    def get_target_distribution(self, long_range_stats_fn='../fess/stats/temp.longrange.contact'):
        '''
        Get the target distribution of long range interaction
        lengths.
        '''
        f = open(long_range_stats_fn, 'r')
        lengths = []

        #length = list(np.linspace(0, 200, 200))
        for line in f:
            parts = line.strip().split(' ')
            lengths += [float(parts[2])]

        print "len(lengths):", len(lengths)
        lengths = lengths[::len(lengths)/100]
        self.fg = cus.interpolated_kde(lengths)

    def calc_fg_parameters(self, bg):
        self.get_target_distribution(cbc.Configuration.longrange_contact_stats_fn)

    def calc_bg_parameters(self, structs):
        '''
        Calculate the energy parameters of a given distribution.

        In this case, the structs parameter contains a list of structures. These structures
        will have a particular distribution of this energy function. The distribution of 
        energies of these structures will be the background distribution.
        
        @param structs: The structures to used to define the background energy distribution.
        '''
        interaction_distances = c.defaultdict(list)

        for bg in structs:

            defines = list(bg.defines.keys())
        
            for j in range(len(defines)):
                for k in range(j+1, len(defines)):
                    if defines[j] not in bg.edges[defines[k]]:
                        interaction = tuple(sorted([defines[j], defines[k]]))

                        distance = cuv.vec_distance(bg.get_point(interaction[0]), bg.get_point(interaction[1]))
                        interaction_distances[interaction] += [distance]

        for interaction in interaction_distances.keys():
            interaction_distances[interaction] += list(np.linspace(0, 200, 100))
            interactions = interaction_distances[interaction][::len(interaction_distances[interaction])/100]
            self.bgs[interaction] = cus.interpolated_kde(interactions)

        #self.plot_energies()

    def plot_energies(self):
        '''
        Make plots of the foreground and background energies.
        '''
        import pylab as pl

        for interaction in self.bgs.keys():
            bg = self.bgs[interaction]
            fg = self.fg

            pl.clf()
            distance_range = np.linspace(0, 100., 300)

            pl.ylim((-20, 10))
            pl.plot(distance_range, fg(distance_range), 'ro')
            pl.plot(distance_range, fg(distance_range) - bg(distance_range), 'go')
            pl.plot(distance_range, bg(distance_range), 'bo')
            figname = 'figs/%s.png' % ("-".join(interaction))
            print >> sys.stderr, "saving: %s..." % (figname)

            pl.savefig(figname, bbox_inches=0)

    def get_energy_contribution(self, bg, interaction, background=True):
        '''
        Get the energy term for an interaction.
        '''

        fg = self.fg
        distance = cuv.vec_distance(bg.get_point(interaction[0]), bg.get_point(interaction[1]))

        bgf = self.bgs[interaction]
        bgp = 1.

        if background:
            #print >>sys.stderr, "distance;", distance, "fg:", fg, "interaction:", interaction


            try:
                fgp = fg(distance)
                #print "distance: ", distance, "fgp:", fgp, "interaction:", interaction
            except FloatingPointError:
                fgp = 1e-200
            bgp = bgf(distance)
        else:
            fgp = fg(distance)
        
        #energy = my_log(fgp) - my_log(bgp)
        energy = fgp - bgp

        return energy

    def iterate_over_interactions(self, bg, background=True):
        defines = list(bg.defines.keys())

        for j in range(len(defines)):
            for k in range(j+1, len(defines)):
                if defines[j] not in bg.edges[defines[k]]:
                    # Ignore interactions with elements that are only length 1
                    if ((bg.defines[defines[j]][1] - bg.defines[defines[j]][0] == 1) or
                        (bg.defines[defines[k]][1] - bg.defines[defines[k]][0] == 1)):
                        continue

                    # Ignore interactions between extremely close elements
                    if bg.bp_distances == None:
                        bg.calc_bp_distances()
                    if  bg.bp_distances[defines[j]][defines[k]] < 10:
                        continue
                
                    interaction = tuple(sorted([defines[j], defines[k]]))

                    energy = self.get_energy_contribution(bg, interaction, background)

                    yield (interaction, energy)

    def prune_energies(self, energies):
        '''
        Take only the three most favorable energy contributions for any
        element.

        if s1 has four interactions ('s3':5, 'b3':4, 'x3':7, 'x4':8) then
        its total energy would be 8 + 7 + 5 = 20.

        @param energies: A dictionary with the interactions (e1, e2) as the key and an
            energy as the value.
        '''
        sorted_interactions_falling = sorted(energies, key=lambda key: -energies[key])
        sorted_interactions_rising = sorted(energies, key=lambda key: energies[key])

        energy_counts = c.defaultdict(int)
        new_energies = dict()

        num_best = 1
        num_worst = 0

        for interaction in sorted_interactions_falling:
            if energy_counts[interaction[0]] < num_best and energy_counts[interaction[1]] < num_best:
                energy_counts[interaction[0]] += 1 
                energy_counts[interaction[1]] += 1
                new_energies[interaction] = energies[interaction]

        energy_counts = c.defaultdict(int)
        for interaction in sorted_interactions_rising:
            if energy_counts[interaction[0]] < num_worst and energy_counts[interaction[1]] < num_worst:
                energy_counts[interaction[0]] += 1 
                energy_counts[interaction[1]] += 1
                new_energies[interaction] = energies[interaction]

        return new_energies
        
    def eval_energy(self, sm, background=True):
        energy_total = 0.
        bg = sm.bg

        energies = dict()

        for (interaction, energy) in self.iterate_over_interactions(bg, background):
            energies[interaction] = energy

        new_energies = self.prune_energies(energies)
        self.interaction_energies = new_energies

        for energy in new_energies.values():
            energy_total += energy

        return -energy_total

class JunctionClosureEnergy(EnergyFunction):
    def __init__(self):
        super(JunctionClosureEnergy, self).__init__()
        self.name = 'jce'
        self.fgs = dict()
        self.bgs = dict()

    def get_target_distributions(self, angle_stats, length):
        '''
        Fit a skew-normal distribution to the distribution of bulge
        lengths.

        @param angle_stats: The statistics file.
        '''

        # we only care about single stranded bulges
        angle_stats = angle_stats[0]

        k = length
        stats = [s.r1 for s in angle_stats[k]]

        if len(stats) < 4:
            fit = [np.mean(stats), 1.0, 0.0]
        else:
            fit = cus.interpolated_kde(stats)

        return fit

    def calc_fg_parameters(self, bg):
        all_bulges = set([d for d in bg.defines.keys() if d[0] != 's' and len(bg.edges[d]) == 2])
        sm = cbm.SpatialModel(bg)

        # build the structure to see if there are any closed bulges
        sm.traverse_and_build()
        closed_bulges = all_bulges.difference(sm.sampled_bulges)

        for bulge in closed_bulges:
            fg_fit = self.get_target_distributions(sm.angle_stats, abs(bg.defines[bulge][1] - bg.defines[bulge][0]))
            self.fgs[abs(bg.defines[bulge][1] - bg.defines[bulge][0])] = fg_fit

    def calc_bg_parameters(self, structs):
        '''
        Calculate the energy parameters of a given distribution.

        In this case, the structs parameter contains a list of structures. These structures
        will have a particular distribution of this energy function. The distribution of 
        energies of these structures will be the background distribution.
        
        @param structs: The structures used to define the background energy distribution.
        '''
        bg = structs[0]
        sm = cbm.SpatialModel(copy.deepcopy(bg))

        sm.traverse_and_build()

        distances = c.defaultdict(list)

        all_bulges = set([d for d in bg.defines.keys() if d[0] != 's' and len(bg.edges[d]) == 2])
        closed_bulges = all_bulges.difference(sm.sampled_bulges)

        for bg in structs:
            for bulge in closed_bulges:
                bl = abs(bg.defines[bulge][1] - bg.defines[bulge][0])
                distance = cuv.vec_distance(bg.coords[bulge][1], bg.coords[bulge][0])
                distances[bulge] += [distance]

        for bulge in closed_bulges:
            bg_fit = cus.interpolated_kde(distances[bulge])

            self.bgs[abs(bg.defines[bulge][1] - bg.defines[bulge][0])] = bg_fit

            '''
            ds = np.array(distances[bulge])

            fg = fg_fit(ds)
            bg = bg_fit(ds)


            pl.plot(ds, fg, 'bo')
            pl.plot(ds, bg, 'ro')
            pl.plot(ds, fg - bg, 'go')
            #show()
            '''

    def eval_energy(self, sm, background=True):
        #bulge = 'b5' 
        bg = sm.bg
        all_bulges = set([d for d in bg.defines.keys() if d[0] != 's' and len(bg.edges[d]) == 2])
        closed_bulges = all_bulges.difference(sm.sampled_bulges)

        energy = np.array([0.])

        for bulge in closed_bulges:
            bl = abs(bg.defines[bulge][1] - bg.defines[bulge][0])

            fgd = self.fgs[bl]
            bgd = self.bgs[bl]

            dist = cuv.vec_distance(bg.coords[bulge][1], bg.coords[bulge][0])
            #print "bl:", bl, "dist:", dist

            if background:
                energy += -(fgd(dist) - bgd(dist))
            else:
                energy += -fgd(dist)
        
        #print "energy:", energy
        #print "energy[0]:", energy[0]

        return energy[0]

class LongRangeInteractionCount(EnergyFunction):
    '''
    An energy function to keep track of how many elements are within
    a certain distance of each other.
    '''

    def __init__(self, di = lri_iter):
        super(LongRangeInteractionCount, self).__init__()
        self.distance_iterator = di
        self.target_interactions = None

    def get_target_interactions(self, bg, filename):
        '''
        Calculate the linear regression of interaction counts.

        @param bg: The BulgeGraph of the target structure
        @param filename: The filename of the statistics file
        '''
        import scipy.stats as ss

        f = open(filename, 'r')
        long_range = []
        all_range = []
        for line in f:
            parts = line.strip().split(' ')

            if float(parts[1]) < 400:
                long_range += [float(parts[0])]
                all_range += [np.sqrt(float(parts[1]))]

        gradient, intercept, r_value, p_value, std_err = ss.linregress(all_range, long_range)

        di = self.distance_iterator
        self.distance_iterator = DistanceIterator()
        total_interactions = self.count_interactions(bg)
        target_interactions = gradient * np.sqrt(total_interactions) + intercept
        self.distance_iterator = di
        
        return target_interactions

    def calc_fg_parameters(self, bg):
        self.target_interactions = self.get_target_interactions(bg, cbc.Configuration.lric_stats_fn)

    def calc_bg_parameters(self, structs):
        '''
        Calculate the energy parameters of a given distribution.

        In this case, the structs parameter contains a list of structures. These structures
        will have a particular distribution of this energy function. The distribution of 
        energies of these structures will be the background distribution.
        
        @param structs: The structures to used to define the background energy distribution.
        '''
        interactions = [self.count_interactions(struct) for struct in structs]
        rand.shuffle(interactions)
        self.bgf = cus.interpolated_kde([float(interaction) for interaction in interactions])
            
    def count_interactions(self, bg):
        '''
        Count the number of long range interactions that occur in the structure.
        '''

        count = 0

        for inter in self.distance_iterator.iterate_over_interactions(bg):
            count += 1

        return count

    def eval_energy(self, sm, background=True):
        bg = sm.bg
        self.distance_iterator = lri_iter
        count = self.count_interactions(bg)
        import scipy.stats as ss

        if self.target_interactions == None:
            raise MissingTargetException("LongRangeInteractionEnergy target_interaction is not defined. This energy probably hasn't been calibrated")

        #return float(count)
        contrib = -(my_log(ss.norm.pdf(float(count), self.target_interactions, 8.)) - self.bgf(float(count)))

        return contrib
        #return -(log(ss.norm.pdf(float(count), 89., 8.)) - log(skew(count, self.skew_fit[0], self.skew_fit[1], self.skew_fit[2])))

class CoarseStemClashEnergy(EnergyFunction):
    '''
    Determine if two stems clash.
    '''
    
    def __init__(self):
        super(CoarseStemClashEnergy, self).__init__()

    def eval_energy(self, sm, background=False):
        bg = sm.bg
        min_distance = 8.45
        energy = 0.
        #print list(bg.stems())

        for (s1, s2) in it.combinations(bg.stems(), 2):
            #print s1, s2
            if bg.are_any_adjacent_stems(s1, s2):
                continue

            closest_points = cuv.line_segment_distance(bg.coords[s1][0],
                                                       bg.coords[s1][1],
                                                       bg.coords[s2][0],
                                                       bg.coords[s2][1])

            closest_distance = cuv.magnitude(closest_points[1] - closest_points[0])
            #print "s1, s2", s1, s2, closest_distance

            if closest_distance < min_distance:
                energy += 100000.0

        return energy

class StemVirtualResClashEnergy(EnergyFunction):
    '''
    Determine if the virtual residues clash.
    '''

    def __init__(self):
        super(StemVirtualResClashEnergy, self).__init__()

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
            kdt2.all_search(1.8)

        clashes = 0
        indeces = kdt2.all_get_indices()
        for (ia,ib) in indeces:
            '''
            cud.pv('(virtual_atoms[ia][1], virtual_atoms[ib][1])')
            print >> sys.stderr, "----------------"
            '''
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

            clashes += 1

        return clashes

    def virtual_residue_atom_clashes(self, bg, s1,i1,a1, s2, i2, a2):
        '''
        Check if any of the virtual residue atoms clash.
        '''
        #(p1, v1, v1_l, v1_r) = cgg.virtual_res_3d_pos(bg, s1, i1)
        #(p2, v2, v2_l, v2_r) = cgg.virtual_res_3d_pos(bg, s2, i2)


        vra1 = self.vras[(s1,i1,a1)]
        vra2 = self.vras[(s2,i2,a2)]

        clashes = 0

        atoms1 = vra1
        atoms2 = vra2

        #for atoms1 in vra1:
            #for atoms2 in vra2:

        for a1 in atoms1.values():
            for a2 in atoms2.values():
                if np.dot(a1-a2, a1-a2) < 1.8 ** 2:
                #if cuv.magnitude(a1 - a2) < 1.8:
                    clashes += 1

        print >>sys.stderr, "clashes1", clashes
        return clashes

    def eval_energy(self, sm, background=False):
        '''
        Cound how many clashes of virtual residues there are.

        @param sm: The SpatialModel containing the list of stems.
        @param background: Use a background distribution to normalize this one.
                           This should always be false since clashes are independent
                           of any other energies.
        '''
        self.vras = dict()
        self.bases = dict()
        self.bg = sm.bg

        l = []
        bg = sm.bg
        mult = 8
        points = []
        energy = 0.

        for d in sm.bg.defines.keys():
            if d[0] == 's':
                s = d
                s_len = bg.stem_length(s)
                #stem_inv = bg.stem_invs[s]

                for i in range(s_len):
                    (p, v, v_l, v_r) = bg.v3dposs[d][i]

                    points += [(p+ mult * v_l, d, i, 1)]
                    points += [(p+ mult * v_r, d, i, 0)]

        coords = np.vstack([p[0] for p in points])

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

            if s1 == s2:
                continue

            if (s1,i1,a1) not in self.vras.keys():
                self.vras[(s1,i1,a1)] = cgg.virtual_residue_atoms(bg, s1, i1, a1)
            if (s2,i2,a2) not in self.vras.keys():
                self.vras[(s2,i2,a2)] = cgg.virtual_residue_atoms(bg, s2, i2, a2)

            #energy += 100000. * self.virtual_residue_atom_clashes(sm.bg, s1, i1, a1, s2, i2, a2)
        energy += 100000. * self.virtual_residue_atom_clashes_kd()

        return energy

class StemClashEnergy(EnergyFunction):
    '''
    Determine if there's any atom clashes in the structures.
    '''

    def __init__(self):
        super(StemClashEnergy, self).__init__()
        self.stem_library = dict()

    def eval_energy(self, sm, background=True):
        '''
        The base energy function simply returns a random number.
        '''
        #chain = rtor.reconstruct_stems(sm, self.stem_library)
        chain = sm.chain
        atoms = bpdb.Selection.unfold_entities(chain, 'A')

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            ns = bpdb.NeighborSearch(atoms)
            contacts1 = len(ns.search_all(0.8))

        return contacts1 * 1000.

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

            d1 = cuv.magnitude(sm.bg.get_point(f) - sm.bg.get_point(t))

            energy += abs(d1 - d)
        
        return energy

class GaussianHelixOrientationEnergy(EnergyFunction):
    def __init__(self):
        super(GaussianHelixOrientationEnergy, self).__init__()
        self.real_kde = self.load_stem_orientation_data('fess/stats/stem_nt.stats')
        self.fake_kde = self.load_stem_orientation_data('fess/stats/stem_nt_sampled.stats')

    def load_stem_orientation_data(self, filename):
        import pandas as pa
        stats = pa.read_csv(filename,header=None, sep=' ')
        t = stats
        points = stats[[t.columns[2], t.columns[3], t.columns[4]]].as_matrix()
        
        return cek.gaussian_kde(points.T)


    def eval_energy(self, sm, background=True):
        bg = sm.bg
        stems = [d for d in bg.defines.keys() if d[0] == 's']
        score = 0

        for s1 in stems:
            s1_len = bg.defines[s1][1] - bg.defines[s1][0]
            for s2 in stems:
                if s1 != s2:
                    s2_len = bg.defines[s2][1] - bg.defines[s2][0]
                    for l in range(s1_len):
                        for k in range(s2_len):
                            r2_spos = cgg.pos_to_spos(bg, s1, k, s2, l)

                            score_incr = my_log(self.real_kde(r2_spos)) - my_log(self.fake_kde(r2_spos))
                            #print

                            score += score_incr
        return -score

class ImgHelixOrientationEnergy(EnergyFunction):
    def __init__(self):
        super(ImgHelixOrientationEnergy).__init__()
        self.res = 2.
        self.real_img, self.real_min_dims = self.load_stem_orientation_data('fess/stats/stem_bulge_nt.stats')
        self.fake_img, self.fake_min_dims = self.load_stem_orientation_data('fess/stats/stem_bulge_nt_sampled.stats')
        '''
        self.real_img, self.real_min_dims = self.load_stem_orientation_data('fess/stats/stem_bulge_nt_truncated.stats')
        self.fake_img, self.fake_min_dims = self.load_stem_orientation_data('fess/stats/stem_bulge_nt_sampled_truncated.stats')
        '''

    def load_stem_orientation_data(self, filename):
        import pandas as pa
        stats = pa.read_csv(filename,header=None, sep=' ')
        t = stats
        points = stats[[t.columns[2], t.columns[3], t.columns[4]]].as_matrix()

        min_dims = np.array([min(points[:,j]) for j in xrange(points.shape[1])])
        max_dims = np.array([max(points[:,j]) for j in xrange(points.shape[1])])

        n_points = [int((max_dims[j] - min_dims[j]) / float(self.res))+1 for j in range(points.shape[1])]

        img = np.zeros(n_points)
        for p in points:
            ixs = [int((p[j] - min_dims[j]) / self.res) for j in xrange(points.shape[1])]
            img[ixs[0],ixs[1],ixs[2]] += 1

        img = sn.gaussian_filter(img, (3,3,3))

        return (img, min_dims)

    def get_img_score(self, points):
        score = 0.
        points = np.array(points)

        for p in points:
            ixs_real = [int((p[j] - self.real_min_dims[j]) / self.res) for j in xrange(points.shape[1])]
            ixs_fake = [int((p[j] - self.fake_min_dims[j]) / self.res) for j in xrange(points.shape[1])]

            try:
                #val_real = my_log(self.real_img[ixs_real[0], ixs_real[1], ixs_real[2]])
                #val_fake = my_log(self.fake_img[ixs_fake[0], ixs_fake[1], ixs_fake[2]])
                val_real = np.log(self.real_img[ixs_real[0], ixs_real[1], ixs_real[2]])
                val_fake = np.log(self.fake_img[ixs_fake[0], ixs_fake[1], ixs_fake[2]])
                #val_fake = 0.

                score += val_real - val_fake
            except IndexError:
                score += -350.

        return score
                    

    def eval_energy(self, sm, background=True):
        bg = sm.bg
        stems = [d for d in bg.defines.keys() if d[0] == 's']
        #stems = [d for d in bg.defines.keys() if (bg.weights[d] == 2 or bg.weights[d] == 0)]
        score = 0.
        points = []
        s1_start = np.zeros(3)
        s1_end = np.zeros(3)
        r2_spos = np.zeros(3)

        vposs = sm.bg.vposs
        vbasis = sm.bg.vbases
        invs = sm.bg.vinvs

        max_distance = 500.

        starts = c.defaultdict( dict )
        ends = c.defaultdict( dict )
        points = []

        self.interaction_energies = c.defaultdict(int)

        if len(bg.vposs.keys()) == 0:
            for stem in bg.stems():
                cgg.add_virtual_residues(bg, stem)

        # pre-calculate all virtual positions and their change of basis matrices
        for s in stems:

            for i in range(bg.stem_length(s)):
                points += [[vposs[s][i], s, i]]

        # pre-calculate the start and end positions of each virtual res
        for s in stems:
            s_len = bg.stem_length(s) # bg.defines[s][1] - bg.defines[s][0] + 1
            s1_0_pos = vposs[s][0]
            s1_len_pos = vposs[s][s_len - 1]

            for i in range(s_len):
                s1_pos = vposs[s][i]

                starts[s][i] = np.dot(invs[s][i], s1_0_pos - s1_pos)
                ends[s][i] = np.dot(invs[s][i], s1_len_pos - s1_pos)

        coords = np.vstack([p[0] for p in points])

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            kdt = kd.KDTree(3)
            kdt.set_coords(coords)

            kdt.all_search(max_distance)
        indices = kdt.all_get_indices()

        energy1 = 0.
        energy2 = 0.
        stem_interactions = dict()
        count = 0

        for s1 in stems:
            for s2 in stems:
                stem_interactions[(s1, s2)] = 0.

        for (ia,ib) in indices:
            point_energies = []

            for (i1, i2) in [(ia, ib), (ib, ia)]:
                point_energy = 0
                s1,l = points[i1][1:]
                s2,k = points[i2][1:]

                s1_pos = vposs[s1][l]

                if s1 != s2 and not bg.are_adjacent_stems(s1, s2):
                    s2_pos = vposs[s2][k]

                    s1_end = ends[s1][l]
                    s1_start = starts[s1][l]

                    np.dot(invs[s1][l], s2_pos - s1_pos, out=r2_spos)

                    #if cuv.magnitude(r2_spos) < max_distance and r2_spos[0] > s1_start[0] and r2_spos[0] < s1_end[0]:
                    #if cuv.magnitude(r2_spos) < max_distance and r2_spos[0] > -3 and r2_spos[0] < 3:
                    if True:
                        point_score = self.get_img_score([r2_spos])
                        #print "point_score:", point_score
                        energy2 += point_score
                        #point_energy += point_score
                        stem_interactions[(s1,s2)] += point_score
                        count += 1
                        #self.interaction_energies[tuple(sorted([s1, s2]))] += -point_score

        energy1 = 0.
        ses = []
        for (s1, s2) in stem_interactions:
            se = min(stem_interactions[(s1, s2)], stem_interactions[(s2, s1)])
            #se = (stem_interactions[(s1,s2)] + stem_interactions[(s2,s1)]) / 2
            energy1 += se
            self.interaction_energies[tuple(sorted([s1,s2]))] = se
            if abs(se) > 0.00001:
                ses += [se]

        score = energy1
        #score = energy2
        if np.isnan(score):
            return 1000000
        return -score

class RoughJunctionClosureEnergy(EnergyFunction):
    def __init__(self):
        super(RoughJunctionClosureEnergy, self).__init__()

    def eval_energy(self, sm, background=True):
        bg = sm.bg
        all_bulges = set([d for d in bg.defines.keys() if d[0] != 's' and (len(bg.edges[d]) == 2 and bg.weights[d] == 1)])
        energy = 0.
        #closed_bulges = all_bulges.difference(sm.sampled_bulges)

        for bulge in all_bulges:
            bl = bg.defines[bulge][1] - bg.defines[bulge][0] - 1
            #dist = cgg.junction_virtual_res_distance(bg, bulge)
            dist = cgg.junction_virtual_atom_distance(bg, bulge)
            
            # 
            #cutoff_distance = (bl) * 5.9 + 13.4
            #cutoff_distance = (bl) * 5.908 + 11.309
            cutoff_distance = (bl) * 6.4 + 6.4


            if (dist > cutoff_distance):
                #print "bulge:", bulge, "bl:", bl, "cutoff_distance:", cutoff_distance, "dist:", dist
                energy += (dist - cutoff_distance) * 10000.

        return energy

class StemStemOrientationEnergy(EnergyFunction):
    def __init__(self, cols=[2]):
        super(StemStemOrientationEnergy, self).__init__()
        self.max_dist = 30
        self.max_lateral_dist = 13.
        self.sample_num = 1000
        #self.col = col
        self.cols = cols

        self.real_data = None
        self.fake_data = None

        '''
        import matplotlib.pyplot as plt
        xs = np.linspace(0,1.57,1000)
        plt.plot(xs, self.real_data(xs), 'g')
        plt.plot(xs, self.fake_data(xs), 'r')
        plt.show() 
        '''

    def load_stem_stem_data(self, filename):
        import pandas as pa
        t = pa.read_csv(filename, header=None, sep=' ')
        sampled_angles = []
        orig_angles = []
        for col in self.cols:
            angles = t[np.all([t[t.columns[0]] < self.max_dist, t[t.columns[4]] < self.max_lateral_dist], axis=0)][t.columns[col]].values
            #sampled_angles += [[rand.choice(angles) for i in range(self.sample_num)]]
            sampled_angles += [angles]
            sampled_angles += [-angles]
            sampled_angles += [2*math.pi - angles]
            orig_angles += [angles]
            cud.pv('angles')
            cud.pv('-angles')
            #sampled_angles += [2*math.pi - angles]

        #return cek.gaussian_kde(sa)
        #return cek.gaussian_kde(sampled_angles, bw_method=0.1)
        k1 = stats.gaussian_kde(orig_angles)
        return stats.gaussian_kde(sampled_angles, bw_method=k1.factor)

    def eval_energy(self, sm, background=True):
        energy = 0
        self.interaction_energies = c.defaultdict(float)
        
        if self.real_data == None:
            col = 0
            self.real_data = self.load_stem_stem_data('fess/stats/stem_stem_orientations.csv')
            self.fake_data = self.load_stem_stem_data('fess/stats/stem_stem_orientations_sampled.csv')
            #self.fake_data = self.load_stem_stem_data('fess/stats/stem_stem_orientations_sampled.csv')

        for (s1,s2) in it.combinations(sm.bg.stems(), r=2):
            orientation = cgg.stem_stem_orientation(sm.bg, s1,s2)
            if orientation[0] < self.max_dist and orientation[4] < self.max_lateral_dist:
                angs = []
                for col in self.cols:
                    sso = cgg.stem_stem_orientation(sm.bg, s1, s2)
                    angs += [sso[col]]

                ang = np.array(angs)
                #ang = min(ang, math.pi - ang)
                real = my_log(self.real_data(ang))
                fake = my_log(self.fake_data(ang))
                #real = my_log( self.real_data(cgg.stem_stem_orientation(sm.bg, s1, s2))[self.col])
                #fake = my_log( self.fake_data(cgg.stem_stem_orientation(sm.bg, s1, s2))[self.col])

                energy += (real - fake)
                #cud.pv('angs, fake, real, real-fake')

                self.interaction_energies[tuple(sorted([s1,s2]))] += (real - fake)

        if energy > 300:
            #pdb.set_trace()
            pass

        return -energy
                
class StemCoverageEnergy(EnergyFunction):
    def __init__(self):
        self.max_dist = 22


    def eval_energy(self, sm, background=True):
        covered = set()
        bg = sm.bg

        for (s1, s2) in it.permutations(bg.stems(), 2):
            if bg.are_any_adjacent_stems(s1, s2):
                continue

            for i in range(bg.stem_length(s1)):
                basis = bg.vbases[s1][i]

                for j in range(bg.stem_length(s2)):
                    np = bg.vposs[s2][j] - bg.vposs[s1][i]
                    new_pos = cuv.change_basis(np, basis, cuv.standard_basis)

                    if abs(new_pos[0]) >= 2.:
                        continue

                    if cuv.magnitude([0, new_pos[1], new_pos[2]]) < self.max_dist:
                        covered.add((s1, i))

        return -math.log(len(covered) + 1)

class CylinderIntersectionEnergy(EnergyFunction):
    def __init__(self):
        self.max_dist = 30.
        self.min_ratio = 0.
        self.max_ratio = 30.

        self.real_data = None
        self.fake_data = None

    def load_data(self, filename):
        import pandas as pa

        t = pa.read_csv(filename, header=None, sep=' ')

        ratios = t[t.columns[1]]
        ratios = list(ratios[~np.isnan(ratios)])

        #new_ratios = [rand.choice(ratios) for i in range(5000)]
        #new_ratios = ratios
        
        #new_ratios += list(np.linspace(self.min_ratio, self.max_ratio, 100))
        new_ratios = ratios

        return (new_ratios, stats.gaussian_kde(new_ratios))

    def calculate_intersection_coverages(self, bg):
        in_cyl_fractions = c.defaultdict(lambda: 0.001)

        for (s1, s2) in it.permutations(bg.stem_like(), 2):
            line = bg.coords[s1]
            cyl = bg.coords[s2]
            extension = 20.

            cyl_vec = cuv.normalize(bg.coords[s2][1] - bg.coords[s2][0])
            cyl = [cyl[0] - extension * cyl_vec,
                   cyl[1] + extension * cyl_vec]
            
            line_len = cuv.magnitude(line[1] - line[0])
            intersects = cuv.cylinder_line_intersection(cyl, line, 
                                                        self.max_dist)
            if len(intersects) > 0 and np.isnan(intersects[0][0]):
                cud.pv('exiting')
                sys.exit(1)

            if len(intersects) == 0:
                in_cyl_len = 0.
            else:
                #in_cyl_len = cuv.magnitude(intersects[1] - intersects[0])
                cyl_basis = cuv.create_orthonormal_basis(cyl_vec)
                intersects_t = cuv.change_basis(intersects.T, cyl_basis, cuv.standard_basis).T
                #in_cyl_len = abs(intersects[1][0] - intersects[0][0])
                in_cyl_len = abs(intersects_t[1][0] - intersects_t[0][0])

            '''
            if s1 == 's4':
                #cud.pv('line')
                #cud.pv('intersects')
                #cud.pv('cyl')
                #cud.pv('in_cyl_len')
                #cud.pv('s2, in_cyl_len / line_len')
            '''

            in_cyl_fractions[s1] += in_cyl_len / line_len
        return in_cyl_fractions

    def eval_energy(self, sm , background=True):
        if self.real_data == None:
            (self.real_data, self.real_kde) = self.load_data('fess/stats/cylinder_intersection_fractions.csv')
            (self.fake_data, self.fake_kde) = self.load_data('fess/stats/cylinder_intersection_fractions_sampled.csv')
            #(self.fake_data, self.fake_kde) = self.load_data('fess/stats/cylinder_intersection_fractions_%s.csv' % (sm.bg.name))

            self.fake_min = min(self.fake_kde(self.fake_data))

            '''
            import matplotlib.pyplot as plt
            fig = plt.figure()

            xs = np.linspace(self.min_ratio, self.max_ratio, 200)

            ax_hist = fig.add_subplot(1,1,1)
            ax_hist.plot(xs, my_log(self.real_data(xs)), 'b')
            ax_hist.plot(xs, my_log(self.fake_data(xs)), 'r')
            
            plt.show()
            '''

        cyl_fractions = self.calculate_intersection_coverages(sm.bg)
        energy = 0.
        for (key,val) in cyl_fractions.items():
            real = my_log(self.real_kde(val))
            fake = my_log(max(self.fake_kde(val), self.fake_min))


            energy += (real - fake)

            #cud.pv('key, val, real, fake, real-fake')

            if np.isnan(energy):
                pdb.set_trace()

        return -energy

class CheatingEnergy(EnergyFunction):
    def __init__(self, real_bg):
        self.real_bg = copy.deepcopy(real_bg)
        self.real_residues = cgg.bg_virtual_residues(self.real_bg)

    def eval_energy(self, sm, background=True):
        new_residues = cgg.bg_virtual_residues(sm.bg)

        return  cbr.centered_rmsd(self.real_residues, new_residues)

class LoopLoopEnergy(EnergyFunction):
    def __init__(self):
        self.real_data = None
        self.fake_data = None

    def load_data(self, filename):
        import pandas as pa

        t = pa.read_csv(filename, header=None, sep=' ')
        t.columns = ['key1', 'type1', 'len1', 'key2', 'type2', 'len2', 'dist', 'seq1', 'seq2', 'longrange', 'angle']

        #loop_loop = t[np.logical_and(t[t.columns[1]] == "l", t[t.columns[4]] == "l")]
        loop_loop = t[np.logical_and(t.type1 == "l", t.type2 == "l")]
        loop_loop_y = loop_loop[loop_loop.longrange == 'Y']
        loop_loop_n = loop_loop[loop_loop.longrange == 'N']

        return (loop_loop.dist.values, cek.gaussian_kde(loop_loop_y.dist), cek.gaussian_kde(loop_loop.dist))

    def eval_energy(self, sm, background=True):
        if self.real_data == None:
            (self.real_data, self.real_d_given_i, self.real_d) = self.load_data('fess/stats/temp.longrange.stats')
            (self.fake_data, self.fake_d_given_i, self.fake_d) = self.load_data('fess/stats/temp.longrange.stats.sampled')

        p_i = 1.

        num = 0
        energy = 0
        for (l1, l2) in it.combinations(sm.bg.loops(), 2):
            (i1,i2) = cuv.line_segment_distance(sm.bg.coords[l1][0], 
                                                sm.bg.coords[l1][1],
                                                sm.bg.coords[l2][0],
                                                sm.bg.coords[l2][1])
            x = cuv.magnitude(i1 - i2)
            if x > 50.:
                x = 50.
            num += 1

            real_i_given_d = self.real_d_given_i(x) / self.real_d(x)
            fake_i_given_d = self.fake_d_given_i(x) / self.fake_d(x)

            #contrib = real_i_given_d - fake_i_given_d
            contrib = my_log(real_i_given_d/fake_i_given_d)
            #cud.pv('l1,l2,contrib,x')

            energy += contrib
        
        if num == 0:
            return 0

        return -energy;
        '''
        x = np.linspace(min(self.real_data), max(self.real_data))

        real_i_given_d = self.real_d_given_i(x) / self.real_d(x)
        fake_i_given_d = self.fake_d_given_i(x) / self.fake_d(x)

        import matplotlib.pyplot as plt

        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        ax.plot(x, my_log(real_i_given_d), 'g')
        ax.plot(x, my_log(fake_i_given_d), 'r')
        ax.plot(x, my_log(real_i_given_d / fake_i_given_d), 'b')
        plt.show()
        '''

class LoopJunctionEnergy(EnergyFunction):
    def __init__(self):
        self.real_data = None
        self.fake_data = None

    def load_data(self, filename):
        import pandas as pa

        t = pa.read_csv(filename, header=None, sep=' ')
        t.columns = ['key1', 'type1', 'len1', 'key2', 'type2', 'len2', 'dist', 'seq1', 'seq2', 'longrange', 'angle']

        #loop_loop = t[np.logical_and(t[t.columns[1]] == "l", t[t.columns[4]] == "l")]
        loop_loop = t[np.logical_and(t.type1 == "l", t.type2 == "m")]
        loop_loop_y = loop_loop[loop_loop.longrange == 'Y']
        loop_loop_n = loop_loop[loop_loop.longrange == 'N']

        return (loop_loop.dist.values, cek.gaussian_kde(loop_loop_y.dist), cek.gaussian_kde(loop_loop.dist))

    def eval_energy(self, sm, background=True):
        if self.real_data == None:
            (self.real_data, self.real_d_given_i, self.real_d) = self.load_data('fess/stats/temp.longrange.stats')
            (self.fake_data, self.fake_d_given_i, self.fake_d) = self.load_data('fess/stats/temp.longrange.stats.sampled')

        p_i = 1.

        num = 0
        energy = 0
        for l1 in sm.bg.loops():
            for l2 in sm.bg.multiloops():
                (i1,i2) = cuv.line_segment_distance(sm.bg.coords[l1][0], 
                                                    sm.bg.coords[l1][1],
                                                    sm.bg.coords[l2][0],
                                                    sm.bg.coords[l2][1])
                x = cuv.magnitude(i1 - i2)
                if x > 50.:
                    x = 50.
                num += 1

                real_i_given_d = self.real_d_given_i(x) / self.real_d(x)
                fake_i_given_d = self.fake_d_given_i(x) / self.fake_d(x)

                #contrib = real_i_given_d - fake_i_given_d
                contrib = my_log(real_i_given_d/fake_i_given_d)
                #cud.pv('l1,l2,contrib,x')

                energy += contrib
            
        if num == 0:
            return 0

        return -energy;
        '''
        x = np.linspace(min(self.real_data), max(self.real_data))

        real_i_given_d = self.real_d_given_i(x) / self.real_d(x)
        fake_i_given_d = self.fake_d_given_i(x) / self.fake_d(x)

        import matplotlib.pyplot as plt

        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        ax.plot(x, my_log(real_i_given_d), 'g')
        ax.plot(x, my_log(fake_i_given_d), 'r')
        ax.plot(x, my_log(real_i_given_d / fake_i_given_d), 'b')
        plt.show()
        '''

class NLoopLoopEnergy(EnergyFunction):
    def __init__(self):
        self.real_dist = None
        self.fake_dist = None

        self.real_struct = cbm.SpatialModel(cgb.BulgeGraph(os.path.join(cbc.Configuration.test_input_dir, "1jj2/graph", "temp.comp")))
        self.fake_struct = cbm.SpatialModel(cgb.BulgeGraph(os.path.join(cbc.Configuration.test_input_dir, "background/1jj2/", "best0.coord")))

        if self.real_dist == None:
            (self.real_dist, self.real_d_given_i, self.real_d, self.real_size, self.real_s_given_i, self.real_s) = self.load_data('fess/stats/temp.longrange.stats')
            (self.fake_dist, self.fake_d_given_i, self.fake_d, self.fake_size, self.fake_s_given_i, self.fake_s) = self.load_data('fess/stats/temp.longrange.stats.sampled')

        # Store the total probabilities for each loop length
        # The loop length is the index and the total probabilites seen
        # are the values
        self.e_reals = c.defaultdict(list)
        self.e_sampleds = c.defaultdict(list)

        e_real = self.all_interaction_probs(self.real_struct)
        e_sampled = self.all_interaction_probs(self.fake_struct)

        for r in e_real:
            self.e_reals[r[0]] += [r[1]]
        for r in e_sampled:
            self.e_sampleds[r[0]] += [r[1]]

        e_reals = self.e_reals
        e_sampleds = self.e_sampleds

        self.ger = c.defaultdict(None)
        self.ges = c.defaultdict(None)

        for (r, l) in e_reals.items():
            if len(e_reals[r]) < 2 or len(e_sampleds[r]) < 2:
                continue

            self.ger[r] = cek.gaussian_kde(self.e_reals[r])
            self.ges[r] = cek.gaussian_kde(self.e_sampleds[r])

        import matplotlib.pyplot as plt
        fig = plt.figure()
        for i in range(16):
            cud.pv('i')
            ax = fig.add_subplot(4,4,i)

            if len(e_reals[i]) < 2 or len(e_sampleds[i]) < 2:
                continue

            cud.pv('e_reals[i]')
            cud.pv('e_sampleds[i]')
            xs = np.linspace(0, max(e_reals[i]), 100)
            ger = self.ger[i]
            ges = self.ges[i]

            ax.plot(xs, my_log(ger(xs)), 'g')
            ax.plot(xs, my_log(ges(xs)), 'r')
            ax.plot(xs, my_log(ger(xs)) - my_log(ges(xs)), 'y')
            ax.set_title(str(i))
        fig.tight_layout()
        plt.show()

    def load_data(self, filename):
        import pandas as pa

        t = pa.read_csv(filename, header=None, sep=' ')
        t.columns = ['key1', 'type1', 'len1', 'key2', 'type2', 'len2', 'dist', 'seq1', 'seq2', 'longrange', 'angle']

        #loop_loop = t[np.logical_and(t[t.columns[1]] == "l", t[t.columns[4]] == "l")]
        loop_loop = t[np.logical_and(t.type1 == "l", t.type2 == "l")]
        loop_loop_y = loop_loop[loop_loop.longrange == 'Y']
        loop_loop_n = loop_loop[loop_loop.longrange == 'N']

        return (loop_loop.dist.values, 
                cek.gaussian_kde(loop_loop_y.dist), 
                cek.gaussian_kde(loop_loop.dist), 
                loop_loop.len1.values, 
                cek.gaussian_kde([float(f) for f in loop_loop_y.len1.values]), 
                cek.gaussian_kde([float(f) for f in loop_loop.len1.values]))

    def interaction_prob(self, sm, l1):
        '''
        Get the total interaction probability of a particular loop.
        '''

        total_p = 0.
        p_i_l1_given_s = (self.real_s_given_i(sm.bg.get_length(l1)) / 
                          self.real_s(sm.bg.get_length(l1)))

        for l2 in sm.bg.loops():
            if l1 == l2:
                continue

            (i1,i2) = cuv.line_segment_distance(sm.bg.coords[l1][0], 
                                                sm.bg.coords[l1][1],
                                                sm.bg.coords[l2][0],
                                                sm.bg.coords[l2][1])
            d = cuv.magnitude(i2 - i1)
            #cud.pv('l1, l2, d')
            if d > 50:
                continue

            p_i_given_d = (self.real_d_given_i(d) / 
                           self.real_d(d))
            p_i_l2_given_s = (self.real_s_given_i(sm.bg.get_length(l2)) / 
                              self.real_s(sm.bg.get_length(l2)))
            contrib = p_i_l1_given_s * p_i_given_d * p_i_l2_given_s
            #total_p += (p_i_given_d * p_i_l2_given_s)
            total_p += contrib[0]

        return total_p
        #total_p *= p_i_l1_given_s

    def all_interaction_probs(self, sm):
        total_ps = []
        for l1 in sm.bg.loops():
            total_p = self.interaction_prob(sm, l1)
            total_ps += [(sm.bg.get_length(l1),total_p)]
            #cud.pv('l1, sm.bg.get_length(l1), total_p')

        return total_ps

    def eval_energy(self, sm, background=True):
        total_ps = self.all_interaction_probs(sm)
        energy = 0.

        for (s, p) in total_ps:
            if len(self.e_reals[s]) < 2 or len(self.e_sampleds[s]) < 2:
                continue

            energy += self.ger[s](p) - self.ges[s](p)

        return -energy
            

class NLoopJunctionEnergy(EnergyFunction):
    def __init__(self):
        self.real_dist = None
        self.fake_dist = None

        self.real_struct = cbm.SpatialModel(cgb.BulgeGraph(os.path.join(cbc.Configuration.test_input_dir, "1jj2/graph", "temp.comp")))
        self.fake_struct = cbm.SpatialModel(cgb.BulgeGraph(os.path.join(cbc.Configuration.test_input_dir, "background/1jj2/", "best0.coord")))

        if self.real_dist == None:
            (self.real_dist, self.real_d_given_i, self.real_d, self.real_size, self.real_s1_given_i, self.real_s2_given_i, self.real_s1, self.real_s2) = self.load_data('fess/stats/temp.longrange.stats')
            (self.fake_dist, self.fake_d_given_i, self.fake_d, self.fake_size, self.fake_s1_given_i, self.fake_s2_given_i, self.real_s1, self.fake_s2) = self.load_data('fess/stats/temp.longrange.stats.sampled')

        # Store the total probabilities for each loop length
        # The loop length is the index and the total probabilites seen
        # are the values
        self.e_reals = c.defaultdict(list)
        self.e_sampleds = c.defaultdict(list)

        e_real = self.all_interaction_probs(self.real_struct)
        e_sampled = self.all_interaction_probs(self.fake_struct)

        for r in e_real:
            self.e_reals[r[0]] += [r[1]]
        for r in e_sampled:
            self.e_sampleds[r[0]] += [r[1]]

        e_reals = self.e_reals
        e_sampleds = self.e_sampleds

        self.ger = c.defaultdict(None)
        self.ges = c.defaultdict(None)

        for (r, l) in e_reals.items():
            if len(e_reals[r]) < 2 or len(e_sampleds[r]) < 2:
                continue

            self.ger[r] = cek.gaussian_kde(self.e_reals[r])
            self.ges[r] = cek.gaussian_kde(self.e_sampleds[r])

        '''
        import matplotlib.pyplot as plt
        fig = plt.figure()
        for i in range(16):
            cud.pv('i')
            ax = fig.add_subplot(4,4,i)

            if len(e_reals[i]) < 2 or len(e_sampleds[i]) < 2:
                continue

            cud.pv('e_reals[i]')
            cud.pv('e_sampleds[i]')
            xs = np.linspace(0, max(e_reals[i]), 100)
            ger = self.ger[i]
            ges = self.ges[i]

            ax.plot(xs, my_log(ger(xs)), 'g')
            ax.plot(xs, my_log(ges(xs)), 'r')
            ax.plot(xs, my_log(ger(xs)) - my_log(ges(xs)), 'y')
            ax.set_title(str(i))
        fig.tight_layout()
        plt.show()
        '''

    def load_data(self, filename):
        import pandas as pa

        t = pa.read_csv(filename, header=None, sep=' ')
        t.columns = ['key1', 'type1', 'len1', 'key2', 'type2', 'len2', 'dist', 'seq1', 'seq2', 'longrange', 'angle']

        #loop_loop = t[np.logical_and(t[t.columns[1]] == "l", t[t.columns[4]] == "l")]
        loop_loop = t[np.logical_and(t.type1 == "m", t.type2 == "l")]
        loop_loop_y = loop_loop[loop_loop.longrange == 'Y']
        loop_loop_n = loop_loop[loop_loop.longrange == 'N']

        return (loop_loop.dist.values, 
                cek.gaussian_kde(loop_loop_y.dist), 
                cek.gaussian_kde(loop_loop.dist), 
                loop_loop.len1.values, 
                cek.gaussian_kde([float(f) for f in loop_loop_y.len1.values]), 
                cek.gaussian_kde([float(f) for f in loop_loop_y.len2.values]), 
                cek.gaussian_kde([float(f) for f in loop_loop.len1.values]),
                cek.gaussian_kde([float(f) for f in loop_loop.len2.values]))

    def interaction_prob(self, sm, l1):
        '''
        Get the total interaction probability of a particular loop.
        '''

        total_p = 0.
        p_i_l1_given_s = (self.real_s1_given_i(sm.bg.get_length(l1)) / 
                          self.real_s1(sm.bg.get_length(l1)))

        for l2 in sm.bg.loops():
            if l1 == l2:
                continue

            (i1,i2) = cuv.line_segment_distance(sm.bg.coords[l1][0], 
                                                sm.bg.coords[l1][1],
                                                sm.bg.coords[l2][0],
                                                sm.bg.coords[l2][1])
            d = cuv.magnitude(i2 - i1)
            #cud.pv('l1, l2, d')
            if d > 50:
                continue

            p_i_given_d = (self.real_d_given_i(d) / 
                           self.real_d(d))
            p_i_l2_given_s = (self.real_s2_given_i(sm.bg.get_length(l2)) / 
                              self.real_s2(sm.bg.get_length(l2)))
            contrib = p_i_l1_given_s * p_i_given_d * p_i_l2_given_s
            #total_p += (p_i_given_d * p_i_l2_given_s)
            total_p += contrib[0]

        return total_p
        #total_p *= p_i_l1_given_s

    def all_interaction_probs(self, sm):
        total_ps = []
        for l1 in sm.bg.multiloops():
            total_p = self.interaction_prob(sm, l1)
            total_ps += [(sm.bg.get_length(l1),total_p)]
            #cud.pv('l1, sm.bg.get_length(l1), total_p')

        #cud.pv('total_ps')
        return total_ps

    def eval_energy(self, sm, background=True):
        total_ps = self.all_interaction_probs(sm)
        energy = 0.

        for (s, p) in total_ps:
            if len(self.e_reals[s]) < 2 or len(self.e_sampleds[s]) < 2:
                continue

            cud.pv('p')
            energy += self.ger[s](p) - self.ges[s](p)

        return -energy

class NLoopStemEnergy(EnergyFunction):
    def __init__(self):
        self.real_dist = None
        self.fake_dist = None

        self.real_struct = cbm.SpatialModel(cgb.BulgeGraph(os.path.join(cbc.Configuration.test_input_dir, "1jj2/graph", "temp.comp")))
        self.fake_struct = cbm.SpatialModel(cgb.BulgeGraph(os.path.join(cbc.Configuration.test_input_dir, "background/1jj2/", "best0.coord")))

        if self.real_dist == None:
            (self.real_dist, self.real_d_given_i, self.real_d, self.real_size, self.real_s1_given_i, self.real_s2_given_i, self.real_s1, self.real_s2, self.real_a_given_i, self.real_a) = self.load_data('fess/stats/temp.longrange.stats')

        # Store the total probabilities for each loop length
        # The loop length is the index and the total probabilites seen
        # are the values
        self.e_reals = c.defaultdict(list)
        self.e_sampleds = c.defaultdict(list)

        e_real = self.all_interaction_probs(self.real_struct)
        e_sampled = self.all_interaction_probs(self.fake_struct)

        for r in e_real:
            self.e_reals[r[0]] += [r[1]]
        for r in e_sampled:
            self.e_sampleds[r[0]] += [r[1]]

        e_reals = self.e_reals
        e_sampleds = self.e_sampleds

        self.ger = c.defaultdict(None)
        self.ges = c.defaultdict(None)

        for (r, l) in e_reals.items():
            if len(e_reals[r]) < 2 or len(e_sampleds[r]) < 2:
                continue

            self.ger[r] = cek.gaussian_kde(self.e_reals[r])
            self.ges[r] = cek.gaussian_kde(self.e_sampleds[r])

        '''
        import matplotlib.pyplot as plt
        fig = plt.figure()
        for i in range(16):
            cud.pv('i')
            ax = fig.add_subplot(4,4,i)

            if len(e_reals[i]) < 2 or len(e_sampleds[i]) < 2:
                continue

            cud.pv('e_reals[i]')
            cud.pv('e_sampleds[i]')
            xs = np.linspace(0, max(e_reals[i]), 100)
            ger = self.ger[i]
            ges = self.ges[i]

            ax.plot(xs, my_log(ger(xs)), 'g')
            ax.plot(xs, my_log(ges(xs)), 'r')
            ax.plot(xs, my_log(ger(xs)) - my_log(ges(xs)), 'y')
            ax.set_title(str(i))
        fig.tight_layout()
        plt.show()
        '''

    def load_data(self, filename):
        import pandas as pa

        t = pa.read_csv(filename, header=None, sep=' ')
        t.columns = ['key1', 'type1', 'len1', 'key2', 'type2', 'len2', 'dist', 'seq1', 'seq2', 'longrange', 'angle']

        #loop_loop = t[np.logical_and(t[t.columns[1]] == "l", t[t.columns[4]] == "l")]
        loop_loop = t[np.logical_and(t.type1 == "l", t.type2 == "s")]
        loop_loop_y = loop_loop[loop_loop.longrange == 'Y']
        loop_loop_n = loop_loop[loop_loop.longrange == 'N']

        cud.pv('loop_loop_y.angle.values')

        return (loop_loop.dist.values, 
                cek.gaussian_kde(loop_loop_y.dist), 
                cek.gaussian_kde(loop_loop.dist), 
                loop_loop.len1.values, 
                cek.gaussian_kde([float(f) for f in loop_loop_y.len1.values]), 
                cek.gaussian_kde([float(f) for f in loop_loop_y.len2.values]), 
                cek.gaussian_kde([float(f) for f in loop_loop.len1.values]),
                cek.gaussian_kde([float(f) for f in loop_loop.len2.values]),
                cek.gaussian_kde([float(f) for f in loop_loop_y.angle.values]),
                cek.gaussian_kde([float(f) for f in loop_loop.angle.values]))

    def interaction_prob(self, sm, l1):
        '''
        Get the total interaction probability of a particular loop.
        '''

        total_p = 0.
        p_i_l1_given_s = (self.real_s1_given_i(sm.bg.get_length(l1)) / 
                          self.real_s1(sm.bg.get_length(l1)))

        for l2 in sm.bg.stems():
            if l1 == l2:
                continue

            (i1,i2) = cuv.line_segment_distance(sm.bg.coords[l1][0], 
                                                sm.bg.coords[l1][1],
                                                sm.bg.coords[l2][0],
                                                sm.bg.coords[l2][1])
            d = cuv.magnitude(i2 - i1)
            #cud.pv('l1, l2, d')
            if d > 50:
                continue

            angle = cgg.receptor_angle(sm.bg, l1, l2)

            p_i_given_d = (self.real_d_given_i(d) / 
                           self.real_d(d))
            p_i_l2_given_s = (self.real_s2_given_i(sm.bg.get_length(l2)) / 
                              self.real_s2(sm.bg.get_length(l2)))
            p_i_given_a = (self.real_a_given_i(angle) / self.real_a(angle))

            contrib = p_i_l1_given_s * p_i_given_d * p_i_l2_given_s * p_i_given_a
            #total_p += (p_i_given_d * p_i_l2_given_s)
            total_p += contrib[0]

        return total_p
        #total_p *= p_i_l1_given_s

    def all_interaction_probs(self, sm):
        total_ps = []
        for l1 in sm.bg.loops():
            total_p = self.interaction_prob(sm, l1)
            total_ps += [(sm.bg.get_length(l1),total_p)]
            #cud.pv('l1, sm.bg.get_length(l1), total_p')

        #cud.pv('total_ps')
        return total_ps

    def eval_energy(self, sm, background=True):
        total_ps = self.all_interaction_probs(sm)
        energy = 0.

        for (s, p) in total_ps:
            if len(self.e_reals[s]) < 2 or len(self.e_sampleds[s]) < 2:
                continue

            energy += self.ger[s](p) - self.ges[s](p)

        return -energy
