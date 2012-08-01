#!/usr/bin/python

import pdb
import scipy.stats as ss
import pickle

from corgy.graph.bulge_graph import BulgeGraph
from corgy.utilities.vector import vec_distance
from bobbins_config import ConstructionConfig

from corgy.builder.models import SpatialModel
from corgy.builder.stats import AngleStatsDict, StemStatsDict, LoopStatsDict

from corgy.utilities.data_structures import DefaultDict
from corgy.utilities.statistics import fit_skew, skew

from pylab import plot,show, hist, savefig, clf, ylim

from numpy import mean

from scipy.stats import norm, linregress
from numpy import log, array, sqrt, linspace
from random import random, shuffle, uniform
from scipy.stats import norm, gaussian_kde

from corgy.utilities.statistics import interpolated_kde

from time import sleep
from sys import float_info, stderr

def my_log(x):
    return log(x + 1e-200)


class DistanceIterator:
    '''
    A class for iterating over the elements that are a certain
    distance from each other.

    '''
    def __init__(self, min_distance=0., max_distance=float_info.max):
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

                dist = vec_distance(point1, point2)

                #if dist > 6.0 and dist < 25.0:
                if dist > self.min_distance and dist < self.max_distance:
                    yield tuple(sorted([d1, d2]))

lri_iter = DistanceIterator(6., 25.)

class CombinedEnergy:
    def __init__(self, energies=[]):
        self.energies = energies

    def eval_energy(self, sm, verbose=False, background=True):
        total_energy = 0.

    
        for energy in self.energies:
            contrib = energy.eval_energy(sm, background)

            total_energy += contrib

            if verbose:
                print energy.__class__.__name__, contrib

        return total_energy

class SkewNormalInteractionEnergy:
    '''
    This energy will assume that all elements need to be a certain
    distance from each other to interact.

    This distance is centered at around 15 angstroms.

    The distribution of distances is modeled by a skew-normal-distribution and
    the total energy will be the sum of log probabilities for the interactions.
    '''
    def __init__(self):
        self.fg = None
        self.bgs = dict()

    def get_target_distribution(self, long_range_stats_fn='../fess/stats/temp.longrange.contact'):
        '''
        Get the target distribution of long range interaction
        lengths.
        '''
        f = open(long_range_stats_fn, 'r')
        lengths = []


        length = list(linspace(0, 200, 200))

        for line in f:
            parts = line.strip().split(' ')
            lengths += [float(parts[2])]

        #self.fg = fit_skew(lengths)
        print "len(lengths):", len(lengths)
        lengths = lengths[::len(lengths)/100]
        self.fg = interpolated_kde(lengths)

    def calibrate(self, bg, steps=40):
        '''
        Run a number of simulations and get the background distribution
        for each interaction pair.
        '''
        interaction_distances = DefaultDict([])
        self.get_target_distribution()
        #self.get_target_distributions(angle_stats)

        sm = SpatialModel(bg)
        sm.traverse_and_build()
        for i in range(steps):
            sm.sample_angles()
            sm.sample_stems()
            sm.sample_loops()
            sm.traverse_and_build()

            defines = list(bg.defines.keys())
        
            for j in range(len(defines)):
                for k in range(j+1, len(defines)):
                    if defines[j] not in bg.edges[defines[k]]:
                        interaction = tuple(sorted([defines[j], defines[k]]))

                        distance = vec_distance(bg.get_point(interaction[0]), bg.get_point(interaction[1]))
                        interaction_distances[interaction] += [distance]

        for interaction in interaction_distances.keys():
            interaction_distances[interaction] += list(linspace(0, 200, 100))
            interactions = interaction_distances[interaction][::len(interaction_distances[interaction])/100]
            #print >>stderr, "interactions:", len(interactions)
            self.bgs[interaction] = interpolated_kde(interactions)
            bg = self.bgs[interaction]
            fg = self.fg

            clf()
            distance_range = linspace(0, 100., 300)

            ylim((-20, 10))
            plot(distance_range, log(fg(distance_range)), 'ro')
            plot(distance_range, log(fg(distance_range)) - log(bg(distance_range)), 'go')
            plot(distance_range, log(bg(distance_range)), 'bo')
            figname = 'figs/%s.png' % ("-".join(interaction))
            print >>stderr, "saving: %s..." % (figname)

            savefig(figname, bbox_inches=0)

    def get_energy_contribution(self, bg, interaction, background=True):
        '''
        Get the energy term for an interaction.
        '''

        fg = self.fg
        distance = vec_distance(bg.get_point(interaction[0]), bg.get_point(interaction[1]))

        bgf = self.bgs[interaction]
        bgp = 1.

        if background:
            #print >>stderr, "distance;", distance, "fg:", fg, "interaction:", interaction


            try:
                fgp = fg(distance)
                #print "distance: ", distance, "fgp:", fgp, "interaction:", interaction
            except FloatingPointError as fpe:
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

                    # Ignore interactions between extremely close elements
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

        energy_counts = DefaultDict(0)
        new_energies = dict()

        num_best = 1
        num_worst = 0

        for interaction in sorted_interactions_falling:
            if energy_counts[interaction[0]] < num_best and energy_counts[interaction[1]] < num_best:
                energy_counts[interaction[0]] += 1 
                energy_counts[interaction[1]] += 1
                new_energies[interaction] = energies[interaction]

        energy_counts = DefaultDict(0)
        for interaction in sorted_interactions_rising:
            if energy_counts[interaction[0]] < num_worst and energy_counts[interaction[1]] < num_worst:
                energy_counts[interaction[0]] += 1 
                energy_counts[interaction[1]] += 1
                new_energies[interaction] = energies[interaction]

        new_energies

        return new_energies

    def iterate_over_interaction_energies(self, bg, background):
        sm = SpatialModel(bg)

        self.eval_energy(sm, background)
        for key in self.interaction_energies.keys():
            yield (key, self.interaction_energies[key])
        
    def eval_energy(self, sm, background=True):
        energy_total = 0.
        interactions = 1.
        bg = sm.bg

        energies = dict()

        for (interaction, energy) in self.iterate_over_interactions(bg, background):
            energies[interaction] = energy

        new_energies = self.prune_energies(energies)
        self.interaction_energies = new_energies

        for energy in new_energies.values():
            energy_total += energy
            interactions += 1

        #return -(energy_total / (2. * interactions))
        return -energy_total

class JunctionClosureEnergy:
    def __init__(self):
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
            fit = [mean(stats), 1.0, 0.0]
        else:
            fit = interpolated_kde(stats)

        return fit

    def calibrate(self, bg, steps = 40):
        sm = SpatialModel(bg)
        distances = DefaultDict([])
        
        sm.traverse_and_build()
        all_bulges = set([d for d in bg.defines.keys() if d[0] != 's' and len(bg.edges[d]) == 2])
        closed_bulges = all_bulges.difference(sm.sampled_bulges)

        for i in range(steps):
            '''
            sm.sample_stems()
            sm.sample_angles()
            sm.sample_loops()
            sm.traverse_and_build()
            '''
            print "step:", i

            filename = "training/best%d.coord" % (i)
            bg = BulgeGraph(filename)

            for bulge in closed_bulges:
                bl = abs(bg.defines[bulge][1] - bg.defines[bulge][0])
                distance = vec_distance(bg.coords[bulge][1], bg.coords[bulge][0])
                distances[bulge] += [distance]
        
        
        for bulge in closed_bulges:
            bg_fit = interpolated_kde(distances[bulge])
            fg_fit = self.get_target_distributions(sm.angle_stats, abs(bg.defines[bulge][1] - bg.defines[bulge][0]))

            self.fgs[abs(bg.defines[bulge][1] - bg.defines[bulge][0])] = fg_fit
            self.bgs[abs(bg.defines[bulge][1] - bg.defines[bulge][0])] = bg_fit

            ds = array(distances[bulge])

            fg = fg_fit(ds)
            bg = bg_fit(ds)


            plot(ds, fg, 'bo')
            plot(ds, bg, 'ro')
            plot(ds, fg - bg, 'go')
            #show()

    def eval_energy(self, sm, background=True):
        #bulge = 'b5' 
        bg = sm.bg
        all_bulges = set([d for d in bg.defines.keys() if d[0] != 's' and len(bg.edges[d]) == 2])
        closed_bulges = all_bulges.difference(sm.sampled_bulges)

        energy = array([0.])

        for bulge in closed_bulges:
            bl = abs(bg.defines[bulge][1] - bg.defines[bulge][0])

            fgd = self.fgs[bl]
            bgd = self.bgs[bl]

            dist = vec_distance(bg.coords[bulge][1], bg.coords[bulge][0])

            if background:
                energy += -(fgd(dist) - bgd(dist))
            else:
                energy += -fgd(dist)
        
        #print "energy:", energy
        #print "energy[0]:", energy[0]

        return energy[0]

class LongRangeInteractionCount:
    def __init__(self, di = lri_iter):
        self.distance_iterator = di
        self.name = 'lric'
        self.gamma_fit = None

    def get_target_interactions(self, bg, filename):
        '''
        Calculate the linear regression of interaction counts.

        @param bg: The BulgeGraph of the target structure
        @param filename: The filename of the statistics file
        '''

        f = open(filename, 'r')
        long_range = []
        all_range = []
        for line in f:
            parts = line.strip().split(' ')

            if float(parts[1]) < 400:
                long_range += [float(parts[0])]
                all_range += [sqrt(float(parts[1]))]

        gradient, intercept, r_value, p_value, std_err = linregress(all_range, long_range)

        di = self.distance_iterator
        self.distance_iterator = DistanceIterator()
        total_interactions = self.count_interactions(bg)
        target_interactions = gradient * sqrt(total_interactions) + intercept
        self.distance_iterator = di
        
        return target_interactions

    def calibrate(self, bg, steps = 40):
        filename = '../fess/stats/temp.energy'
        self.target_interactions = self.get_target_interactions(bg, filename)

        sm = SpatialModel(bg)
        energies = []

        ld = len(bg.defines.keys())
        '''
        
        for i in range(ld, ld * ld / 2):
            energies += [i]
        '''

        for i in range(steps):
            '''
            sm.sample_stems()
            sm.sample_angles()
            sm.sample_loops()
            sm.traverse_and_build()
            '''
            filename = "training/best%d.coord" % (i)
            bg = BulgeGraph(filename)

            #energies += [self.count_interactions(sm.bg)]
            energies += [self.count_interactions(bg)]

        #energies += range(10, 250)
        #fit = ss.gamma.fit(energies)
        #fit = ss.gamma.fit(energies)

        #fit = distrib.fit(energies)
        #fit = fit_skew(energies)
        shuffle(energies)
        kde = interpolated_kde([float(energy) for energy in energies[::len(energies)/1000]])
        fit = kde


        #hist(energies, normed=True)
        #plot(energies, ss.gamma.pdf(energies, fit[0], fit[1], fit[2]), 'go')
        #plot(energies, distrib.pdf(energies, fit[0], fit[1], fit[2]), 'go')
        #plot(energies, skew(energies, fit[0], fit[1], fit[2]), 'go')
        #hist(kde(energies), normed=True)
        #plot(energies, kde(energies), 'go')
        #plot(energies, kde(energies), 'go')

        
        print("fit:", fit)
        #show()

        
        energy_range = linspace(min(energies), 250, 1000)

        ylim((-40, 20))

        #plot(energies, log(norm.pdf(energy_range, self.target_interactions, 8.)) - log(kde(energy_range)), 'go')
        plot(energy_range, log(norm.pdf(energy_range, self.target_interactions, 8.)) , 'ro')
        plot(energy_range, log(norm.pdf(energy_range, self.target_interactions, 8.)) - log(kde(energy_range)) , 'go')

        print >>stderr, energies
        print >>stderr, fit, min(energies)

        self.bgf = fit

        for count in range(50, 220,5):
            print "fg - bg:", count, my_log(norm.pdf(float(count), self.target_interactions, 5.)), my_log(self.bgf(float(count))), (my_log(norm.pdf(float(count), self.target_interactions, 5.)) - my_log(self.bgf(float(count))))

        #show()
            
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

        #return float(count)
        contrib = -(my_log(norm.pdf(float(count), self.target_interactions, 8.)) - self.bgf(float(count)))

        return contrib
        #return -(log(norm.pdf(float(count), 89., 8.)) - log(skew(count, self.skew_fit[0], self.skew_fit[1], self.skew_fit[2])))
        

class LongRangeDistanceEnergy:
    def __init__(self):
        self.calibration_size = 1000
        self.name = 'lrde'

        pass

    def calibrate(self, bg, steps = None ):
        '''
        Sample a bunch of structure and see which elements tend to be
        near each by virtue of the sampling technique.
        '''

        if steps == None:
            steps = self.calibration_size

        sm = SpatialModel(bg)
        interactions = DefaultDict(0)

        for i in range(steps):
            sm.sample_stems()
            sm.sample_angles()
            sm.traverse_and_build()

            bg = sm.bg

            for i in lri_iter.iterate_over_interactions(bg):
                interactions[i] += 1
        
        self.energies = dict()

        for key in interactions.keys():
            self.energies[key] = 1 / float(interactions[key])


    def eval_energy(self, sm, background=True):
        '''
        Evaluate the energy of a coarse-grained structure.

        @param bg: The representation of the coarse-grained RNA molecule.
        '''
        bg = sm.bg
        energy = 0.

        for interaction in lri_iter.iterate_over_interactions(bg):
            try:
                energy += self.energies[interaction]
            except KeyError:
                energy += 1 / float(self.calibration_size)

        return -my_log(energy)

    def naive_energy(self, bg):
        keys = bg.defines.keys()
        energy = 0.

        for i in range(len(keys)):
            for j in range(i+1, len(keys)):
                energy += 1

        return -energy

class RandomEnergy():
    '''
    An energy function that always just returns a random value.
    '''
    def __init__(self):
        pass

    def eval_energy(self, sm, background=True):
        return uniform(-5, 3)
