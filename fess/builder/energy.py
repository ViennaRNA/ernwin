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
import forgi.threedee.utilities.graph_pdb as cgg
import forgi.threedee.utilities.graph_pdb as ftug
import fess.builder.aminor as fba
import fess.builder.models as cbm
import forgi.utilities.debug as fud
import forgi.threedee.utilities.rmsd as cbr

import scipy.stats as stats
import scipy.stats as ss
import sys


distribution_upper_bound = 1.0
distribution_lower_bound = 1.0

incr = 0.01

def my_log(x):
    return np.log(x + 1e-200)


class MissingTargetException(Exception):
    def __init__(self, message):
        Exception.__init__(self, message)

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


class EnergyFunction(object):
    '''
    The base class for energy functions.
    '''

    def __init__(self):
        self.interaction_energies = None
        #: Used by constraint energies, to store tuples of stems that clash.
        #: Updated every time eval_energy is called.
        self.bad_bulges = []

        self.measures = []
        self.accepted_measures = []

    def accept_last_measure(self):
        if len(self.measures) > 0:
            self.accepted_measures += [self.measures[-1]]

    def reject_last_measure(self):
        if len(self.accepted_measures) > 0:
            self.accepted_measures += [self.accepted_measures[-1]]

    def eval_energy(self, sm, background=True, nodes=None, new_nodes=None):
        '''
        The base energy function simply returns a random number.
        '''
        return rand.random()

    def resample_background_kde(self, struct):
        pass

    def interaction_energy_iter(self, bg, background):
        sm = cbm.SpatialModel(bg)

        for stem in bg.stems():
            ftug.add_virtual_residues(bg, stem)

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

    def get_energy_name(self):
        '''
        Return the name of the energy.
        '''

        return self.__class__.__name__.lower() + ".measures"

    def dump_measures(self, base_directory, iteration=None):
        '''
        Dump all of the coarse grain measures collected so far
        to a file.
        '''
        output_file = op.join(base_directory, self.get_energy_name())

        with open(output_file, 'w') as f:
            f.write(" ".join(map("{:.2f}".format,self.accepted_measures)))
            f.write("\n")

        if iteration is not None:
            with open(output_file + ".%d" % (iteration), 'w') as f:
                f.write(" ".join(map("{:.2f}".format,self.accepted_measures)))
                f.write("\n")

class CoarseGrainEnergy(EnergyFunction):
    def __init__(self, energy_prefactor=10):
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

        pass

    def resample_background_kde(self, struct):
        values = self.accepted_measures

        if len(values) > 100:
            new_kde = self.get_distribution_from_values(values)

            if new_kde is not None:
                self.sampled_kdes[self.measure_category(struct)] = new_kde
            else:
                print >>sys.stderr, "skipping this time..."

            #self.vals[1] = values

    def get_distribution_from_values(self, values):
        '''
        Return a probability distribution from the given values.

        @param values: The values to fit a distribution to.
        @return: A probability distribution fit to the values.
        '''
        floc = -0.1
        fscale =  1.5 * max(values)

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

    def eval_energy(self, sm, background=True, nodes=None, new_nodes=None):
        '''
        A generic function which simply evaluates the energy based on the
        previously calculated probability distributions.
        '''
        '''
        print "-------------------------------"
        for line in traceback.format_stack():
                    print line.strip()
        '''
        cg = sm.bg

        if self.measure_category(cg) not in self.real_kdes.keys():
            (self.real_kdes[self.measure_category(cg)], x) = self.get_distribution_from_file(self.real_stats_fn, 
                                                                                             self.measure_category(cg))
            (self.sampled_kdes[self.measure_category(cg)], self.measures) = self.get_distribution_from_file(self.sampled_stats_fn, self.measure_category(cg))

            self.accepted_measures = self.measures[:]

        kr = self.real_kdes[self.measure_category(cg)]
        ks = self.sampled_kdes[self.measure_category(cg)]

        m = self.get_cg_measure(sm)
        self.measures.append(m)

        if background:
            energy = (np.log(kr(m) + 0.00000001 * ks(m)) - np.log(ks(m)))
            self.prev_energy = energy
            self.prev_cg = m
            #energy = (np.log(kr.integrate_box_1d(0., m) + 0.0001 * ks.integrate_box_1d(0., m)) - np.log(ks.integrate_box_1d(0., m)))
            return -1 * self.energy_prefactor * energy
        else:
            energy = np.log(kr(m))
            return -energy


class ConstantEnergy(EnergyFunction):
    '''
    An energy function that always just returns a constant value (0.).
    '''
    def __init__(self):
        super(ConstantEnergy, self).__init__()

    
    def eval_energy(self, sm, background=True, nodes=None, new_nodes=None):
        return 0.

class RandomEnergy(EnergyFunction):
    '''
    An energy function that always just returns a random value.
    '''
    def __init__(self):
        super(RandomEnergy, self).__init__()

    
    def eval_energy(self, sm, background=True, nodes=None, new_nodes=None):
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
            for j in range(i + 1, len(keys)):
                d1 = keys[i]
                d2 = keys[j]

                point1 = bg.get_point(d1)
                point2 = bg.get_point(d2)

                dist = ftuv.vec_distance(point1, point2)

                #if dist > 6.0 and dist < 25.0:
                if dist > self.min_distance and dist < self.max_distance:
                    yield tuple(sorted([d1, d2]))


class CombinedEnergy:
    def __init__(self, energies=[], uncalibrated_energies=[]):
        self.energies = energies
        self.uncalibrated_energies = uncalibrated_energies
        self.bad_bulges = []

    def save_energy(self, energy, directory):
        if not os.path.exists(directory):
            os.makedirs(directory)

        filename = os.path.join(directory,
                                energy.__class__.__name__ + ".energy")
        print "saving filename:", filename
        pickle.dump(energy, open(filename, 'w'))

    def resample_background_kde(self, struct):
        for e in it.chain(self.energies,
                          self.uncalibrated_energies):
            e.resample_background_kde(struct)
        pass

    def accept_last_measure(self):
        for e in it.chain(self.energies, self.uncalibrated_energies):
            e.accept_last_measure()

    def reject_last_measure(self):
        for e in it.chain(self.energies, self.uncalibrated_energies):
            e.reject_last_measure()

    def dump_measures(self, base_directory, iteration=None):
        for e in it.chain(self.energies,                       
                          self.uncalibrated_energies): 
            e.dump_measures(base_directory, iteration)

    def eval_energy(self, sm, verbose=False, background=True,
                    nodes=None, new_nodes=None):
        total_energy = 0.
        self.bad_bulges = []

        for energy in self.uncalibrated_energies:
            contrib = energy.eval_energy(sm, background,
                                               nodes, new_nodes)
            total_energy += contrib
            self.bad_bulges += energy.bad_bulges
            if verbose:
                print energy.__class__.__name__, contrib

        for energy in self.energies:
            contrib = energy.eval_energy(sm, background=background, nodes=nodes, new_nodes=new_nodes)

            self.bad_bulges += energy.bad_bulges
            total_energy += contrib

            if verbose:
                print energy.__class__.__name__, contrib

        if verbose:
            #import traceback

            #def f():
            #    g()

#            def g():
#                for line in traceback.format_stack():
#                    print line.strip()

            #f()

            print "--------------------------"
            print "total_energy:", total_energy

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

    def eval_energy(self, sm, background=False, nodes=None, new_nodes=None):
        return 0.
        self.last_clashes=[]
        bg = sm.bg
        min_distance = 8.45
        energy = 0.
        #print list(bg.stems())

        if nodes == None:
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
                if ftuv.vec_distance(a1,a2) < 1.8: 
                #if ftuv.magnitude(a1 - a2) < 1.8:
                    clashes += 1

        return clashes

    def eval_energy(self, sm, background=False, nodes = None, new_nodes = None):
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
        self.last_clashes=[]
        bg = sm.bg
        mult = 8
        points = []
        energy = 0.

        if nodes == None:
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

        if new_nodes == None:
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
                clash_pairs += [((s1,i1,a2), (s2,i2,a2))]
        else:
            new_points = []
            indeces = []
            clash_pairs = []
            for d in new_nodes:
                if d[0] == 's':
                    s = d
                    s_len = bg.stem_length(s)
                    #stem_inv = bg.stem_invs[s]

                    for i in range(s_len):
                        (p, v, v_l, v_r) = bg.v3dposs[d][i]

                        new_points += [(p+ mult * v_l, d, i, 1)]
                        new_points += [(p+ mult * v_r, d, i, 0)]

            for i,p in enumerate(points):
                for ni, newp in enumerate(new_points):
                    if p[1] == newp[1]:
                        continue
                    if ftuv.vec_distance(p[0],newp[0]) < 10.:
                        clash_pairs += [(p[1:], newp[1:])]

        potential_clashes = 0
        for (s1, i1, a1), (s2,i2,a2) in clash_pairs:

            if new_nodes != None:
                if s1 not in new_nodes and s2 not in new_nodes:
                    continue

            if s1 == s2:
                continue

            if len(set.intersection(bg.edges[s1], bg.edges[s2])) > 0:
                # the stems are connected
                continue

            potential_clashes += 1
            #fud.pv('s1,s2')

            if (s1,i1,a1) not in self.vras.keys():
                self.vras[(s1,i1,a1)] = cgg.virtual_residue_atoms(bg, s1, i1, a1)
            if (s2,i2,a2) not in self.vras.keys():
                self.vras[(s2,i2,a2)] = cgg.virtual_residue_atoms(bg, s2, i2, a2)

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

    def eval_energy(self, sm, background=True, nodes=None, new_nodes=None):
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
            #dist = cgg.junction_virtual_res_distance(bg, bulge)
            dist = cgg.junction_virtual_atom_distance(bg, bulge)

            #
            #cutoff_distance = (bl) * 5.9 + 13.4
            #cutoff_distance = (bl) * 5.908 + 11.309
            #cutoff_distance = (bl) * 6.4 + 6.4
            cutoff_distance = (bl) * 6.22 + 14.0
            # Note: DOI: 10.1021/jp810014s claims that a typical MeO-P bond is 1.66A long. 
            if (dist > cutoff_distance):
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

    def eval_energy(self, sm, background=True, nodes=None, new_nodes=None):
        closest_distance = self.get_distance(sm)

        if closest_distance < self.distance:
            return 0

        energy = -np.log(self.expon.pdf(closest_distance)) * 10.

        return energy

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
        self.real_residues = cgg.bg_virtual_residues(self.real_bg)

    def eval_energy(self, sm, background=True, nodes=None, new_nodes=None):
	'''
	@param sm: A SpatialModel, which contains a coarse grain model (sm.bg)
	'''
        new_residues = cgg.bg_virtual_residues(sm.bg)

        return  cbr.centered_rmsd(self.real_residues, new_residues)

def get_coords(cg):
    '''
    Return a list of all the coordinates in this structure.
    '''
    coords = []
    for s in cg.sorted_stem_iterator():
        coords += [cg.coords[s][0]]
        coords += [cg.coords[s][1]]
    return coords

import forgi.threedee.utilities.rmsd as ftur

def length_and_rog(cg):
    coords = get_coords(cg)
    rog = ftur.radius_of_gyration(coords)
    total_length = sum([len(list(cg.define_residue_num_iterator(d))) for d in cg.defines.keys()])
    
    return (total_length, rog)

def length_and_rog_from_file(filename):
    cg = ftmc.CoarseGrainRNA(op.expanduser(filename))
    
    return length_and_rog(cg)

class SimpleRadiusOfGyrationEnergy(EnergyFunction):
    def __init__(self):
        super(SimpleRadiusOfGyrationEnergy, self).__init__()

    def eval_energy(self, sm, background=True, nodes=None, new_nodes=None):
        cg = sm.bg
        (length, rog) = length_and_rog(cg)

        return -rog
    
class RadiusOfGyrationEnergy(CoarseGrainEnergy):
    def __init__(self, dist_type="kde", adjustment=1., energy_prefactor=30):
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
        self.adjustment = adjustment # the adjustment is used to enlarge or shrink
                                     # a particular distribution

    def get_name(self):
        return "Radius Of Gyration"

    def get_cg_measure(self, sm):
        (length, rog) = length_and_rog(sm.bg)

        self.measures += [rog]
        return rog

    def get_distribution_from_file(self, filename, length):
        data = np.genfromtxt(load_local_data(filename), delimiter=' ')

        rdata = []

        distribution_upper_bound = 1.0
        distribution_lower_bound = 1.0

        while (len(rdata) < 500):
            distribution_lower_bound -= incr
            distribution_upper_bound += incr

            rdata = data[np.logical_and( data[:,0] > ( distribution_lower_bound ) * length,
                                         data[:,0] < length * ( distribution_upper_bound ))]

        rogs = rdata[:,1]
        return (self.get_distribution_from_values(rogs), list(rogs))


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


class ShortestLoopDistancePerLoop(ShortestLoopDistanceEnergy):
    def __init__(self, loop_name):
        super(ShortestLoopDistancePerLoop, self).__init__()
        self.loop_name = loop_name

        self.real_stats_fn = 'stats/loop_loop3_distances_native.csv'
        self.sampled_stats_fn = 'stats/loop_loop3_distances_sampled.csv'

    def measure_category(self, cg):
        return 1

    def get_distribution_from_file(self, filename, length):
        data = np.genfromtxt(load_local_data(filename), delimiter=' ')

        rogs = data
        return (self.get_distribution_from_values(rogs), list(rogs))

    def get_cg_measure(self, sm):
        #import traceback
        min_dist = 10000.

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

        return min_dist

    def get_energy_name(self):
        '''
        Return the name of the energy.
        '''

        return self.__class__.__name__.lower() + "_" + self.loop_name + ".measures"

    def eval_energy(self, sm, background=True, nodes=None, new_nodes=None):
        '''
        We want to return an energy of 0. if there's less than two hairpin
        loops.
        '''
        if len(list(sm.bg.hloop_iterator())) < 2:
            return 0.
        else:
            return super(ShortestLoopDistancePerLoop, self).eval_energy(sm,
                                                                        background,
                                                                        nodes,
                                                                        new_nodes)

        
class AMinorEnergy(CoarseGrainEnergy):
    def __init__(self, dist_type="kde", adjustment=1., energy_prefactor=30, loop_type='h'):
        import pandas as pd
        super(AMinorEnergy, self).__init__(energy_prefactor=energy_prefactor)

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

    def get_distribution_from_file(self, filename, length):
        data = np.genfromtxt(load_local_data(filename), delimiter=' ')

        rdata = list()

        distribution_lower_bound = 1.
        distribution_upper_bound = 1.

        while (len(rdata) < 500):
            distribution_lower_bound -= incr
            distribution_upper_bound += incr

            rdata = data[np.logical_and( data[:,0] > ( distribution_lower_bound ) * length,
                                         data[:,0] < length * ( distribution_upper_bound ))]
        
        srdata = rdata[rdata[:,1] == self.loop_type]
        rogs = srdata[:,2]

        return (self.get_distribution_from_values(rogs), list(rogs))

    def eval_prob(self, cg, d):
        lt = d[0]
        prob = 0.
        stem_counts = 0
        probs = []

        for s in cg.stem_iterator():
            if s in cg.edges[d]:
                continue

            if ftug.element_distance(cg, d, s) > 30:
                continue

            point = fba.get_relative_orientation(cg, d, s)
            stem_counts += 1
            p = self.prob_funcs[lt](point)
            prob += p
            probs += [p]

        if len(probs) == 0:
            return np.array([0.])

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
        energy = 0

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
                #energy = (np.log(kr.integrate_box_1d(0., m) + 0.0001 * ks.integrate_box_1d(0., m)) - np.log(ks.integrate_box_1d(0., m)))
                energy +=  -1 * self.energy_prefactor * prev_energy
            else:
                energy +=  -np.log(kr(m))
        
        return energy

class SpecificAMinorEnergy(AMinorEnergy):
    def __init__(self, dist_type="kde", adjustment=1., energy_prefactor=30, loop_name='h0'):
        super(SpecificAMinorEnergy, self).__init__(energy_prefactor=energy_prefactor)

        self.real_stats_fn = 'real'
        self.sampled_stats_fn = 'sampled'
        self.loop_name = loop_name

    def get_distribution_from_file(self, filename, category):
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
        energy = 0

        d = self.loop_name

        m = self.eval_prob(sm.bg, d)[0]
        self.measures.append(m)

        if background:
            prev_energy = (np.log(kr(m) + 0.00000001 * ks(m)) - np.log(ks(m)))
            self.prev_energy = energy
            self.prev_cg = m
            #energy = (np.log(kr.integrate_box_1d(0., m) + 0.0001 * ks.integrate_box_1d(0., m)) - np.log(ks.integrate_box_1d(0., m)))
            energy +=  -1 * self.energy_prefactor * prev_energy
        else:
            energy +=  -np.log(kr(m))
        
        return energy
