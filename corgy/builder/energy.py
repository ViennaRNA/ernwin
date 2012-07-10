#!/usr/bin/python

import pdb

from corgy.utilities.vector import vec_distance
from bobbins_config import ConstructionConfig

from corgy.builder.models import SpatialModel
from corgy.builder.stats import AngleStatsDict, StemStatsDict, LoopStatsDict

from corgy.utilities.data_structures import DefaultDict
from corgy.utilities.statistics import fit_skew, skew

from pylab import plot,show, hist, savefig, clf

from numpy import mean

from scipy.stats import norm
from numpy import log, array

from time import sleep
from sys import float_info, stderr

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
                    yield (d1, d2)

lri_iter = DistanceIterator(6., 25.)

class CombinedEnergy:
    def __init__(self, energies=[]):
        self.energies = energies

    def eval_energy(self, bg, verbose=False, background=True):
        total_energy = 0.

    
        for energy in self.energies:
            contrib = energy.eval_energy(bg, background)

            #print "name", energy.__class__.__name__, "contrib:", contrib
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

        for line in f:
            parts = line.strip().split(' ')
            lengths += [float(parts[2])]

        self.fg = fit_skew(lengths)

    def calibrate(self, bg, steps=40):
        '''
        Run a number of simulations and get the background distribution
        for each interaction pair.
        '''
        angle_stats = AngleStatsDict(ConstructionConfig.angle_stats_file)
        stem_stats = StemStatsDict(ConstructionConfig.angle_stats_file)
        loop_stats = LoopStatsDict(ConstructionConfig.angle_stats_file)

        interaction_distances = DefaultDict([])

        self.get_target_distribution()
        #self.get_target_distributions(angle_stats)

        sm = SpatialModel(bg, angle_stats, stem_stats, loop_stats)
        sm.traverse_and_build()
        for i in range(steps):
            print >>stderr, "step:", i
            sm.sample_angles()
            sm.sample_stems()
            sm.sample_loops()
            sm.traverse_and_build()

            defines = list(bg.defines.keys())
            for j in range(len(defines)):
                for k in range(j+1, len(defines)):
                    if defines[j] not in bg.edges[defines[k]]:
                        interaction = (defines[j], defines[k])

                        distance = vec_distance(bg.get_point(interaction[0]), bg.get_point(interaction[1]))
                        interaction_distances[interaction] += [distance]

        for interaction in interaction_distances.keys():
            self.bgs[interaction] = fit_skew(interaction_distances[interaction])
            if self.bgs[interaction][1] < 6.:
                print >>stderr, "low standard deviation... changing to 6."
                self.bgs[interaction][1] = 6.

            print >>stderr, "self.bgs[interaction]:", self.bgs[interaction]
            '''
            clf()
            hist(interaction_distances[interaction])
            savefig('figs/%s.png' % ("-".join(interaction)), bbox_inches=0)
            '''


    def get_energy_contribution(self, bg, interaction, background=True):
        '''
        Get the energy term for an interaction.
        '''

        fg = self.fg
        distance = vec_distance(bg.get_point(interaction[0]), bg.get_point(interaction[1]))

        bgf = self.bgs[interaction]
        bgp = 1.

        if distance < 30.:
            if background:
                fgp = skew(distance, fg[0], fg[1], fg[2]) 
                bgp = skew(distance, bgf[0], bgf[1], bgf[2])
            else:
                fgp = skew(distance, fg[0], fg[1], fg[2])
        else:
            distance = 30
            fgp = skew(30, fg[0], fg[1], fg[2])

        if bgp == 0.:
            #print >>stderr, "abnormally low bgp... changing to slightly higher"
            bgp = .0000000000000000000000001

        energy = log(fgp) - log(bgp)

        return energy

    def iterate_over_interactions(self, bg, background=True):
        defines = list(bg.defines.keys())
        #skew30 = skew(30, fg[0], fg[1], fg[2])

        for j in range(len(defines)):
            for k in range(j+1, len(defines)):
                if defines[j] not in bg.edges[defines[k]]:
                    interaction = (defines[j], defines[k])

                    energy = self.get_energy_contribution(bg, interaction, background)
                    #print >>stderr, "interaction:", interaction, "energy:", energy

                    yield (interaction, energy)
        
    def eval_energy(self, bg, background=True):
        energy_total = 0.
        interactions = 1.

        for (interaction, energy) in self.iterate_over_interactions(bg, background):
            energy_total += energy
            interactions += 1.

        #print "energy_total:", energy_total, "interactions:", interactions, "et / int:", energy_total / interactions
        #return -(energy_total / interactions)
        return -(energy_total / interactions)

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
            fit = fit_skew(stats)

        return fit

    def calibrate(self, bg, steps = 40):
        angle_stats = AngleStatsDict(ConstructionConfig.angle_stats_file)
        stem_stats = StemStatsDict(ConstructionConfig.angle_stats_file)
        loop_stats = LoopStatsDict(ConstructionConfig.angle_stats_file)

        #self.get_target_distributions(angle_stats)

        sm = SpatialModel(bg, angle_stats, stem_stats, loop_stats)
        all_bulges = set([d for d in bg.defines.keys() if d[0] != 's' and len(bg.edges[d]) == 2])
        distances = DefaultDict([])
        
        sm.traverse_and_build()
        closed_bulges = all_bulges.difference(sm.sampled_bulges)

        for i in range(steps):
            sm.sample_stems()
            sm.sample_angles()
            sm.traverse_and_build()


            for bulge in closed_bulges:
                bl = abs(bg.defines[bulge][1] - bg.defines[bulge][0])
                distance = vec_distance(bg.coords[bulge][1], bg.coords[bulge][0])
                distances[bulge] += [distance]
                print "bulge", bulge, bl, distance

        
        
        for bulge in closed_bulges:
            bg_fit = fit_skew(distances[bulge])
            fg_fit = self.get_target_distributions(angle_stats, abs(bg.defines[bulge][1] - bg.defines[bulge][0]))

            self.fgs[abs(bg.defines[bulge][1] - bg.defines[bulge][0])] = fg_fit
            self.bgs[abs(bg.defines[bulge][1] - bg.defines[bulge][0])] = bg_fit

            ds = array(distances[bulge])

            fg = log(skew(ds, fg_fit[0], fg_fit[1], fg_fit[2]))
            bg = log(skew(ds, bg_fit[0], bg_fit[1], bg_fit[2]))

            print "bulge", bulge, "bg_fit", bg_fit, "fg_fit", fg_fit

            plot(ds, fg, 'bo')
            plot(ds, bg, 'ro')
            plot(ds, fg - bg, 'go')
            show()

    def eval_energy(self, bg, background=True):
        bulge = 'b5' 
        bl = abs(bg.defines[bulge][1] - bg.defines[bulge][0])

        fgd = self.fgs[bl]
        bgd = self.bgs[bl]

        dist = vec_distance(bg.coords[bulge][1], bg.coords[bulge][0])

        if background:
            return -(log(skew(dist, fgd[0], fgd[1], fgd[2])) - log(skew(dist, bgd[0], bgd[1], bgd[2])))
        else:
            return -(log(skew(dist, fgd[0], fgd[1], fgd[2])))




class LongRangeInteractionCount:
    def __init__(self, di = lri_iter):
        self.distance_iterator = di
        self.name = 'lric'
        self.skew_fit = None

    def calibrate(self, bg, steps = 40):
        angle_stats = AngleStatsDict(ConstructionConfig.angle_stats_file)
        stem_stats = StemStatsDict(ConstructionConfig.angle_stats_file)
        loop_stats = LoopStatsDict(ConstructionConfig.angle_stats_file)

        sm = SpatialModel(bg, angle_stats, stem_stats, loop_stats)
        energies = []

        for i in range(steps):
            sm.sample_stems()
            sm.sample_angles()
            sm.traverse_and_build()

            energies += [self.eval_energy(sm.bg)]

        fit = fit_skew(energies)

        #hist(energies, bins=len(set(energies))-1)
        #plot(energies, len(energies) * skew(energies, fit[0], fit[1], fit[2]), 'o')
        #show()

        print stderr, "fit:", fit
        self.skew_fit = fit
            
    def count_interactions(self, bg):
        '''
        Count the number of long range interactions that occur in the structure.
        '''

        count = 0

        for inter in self.distance_iterator.iterate_over_interactions(bg):
            count += 1

        return count

    def eval_energy(self, bg, background=True):
        count = self.count_interactions(bg)

        #return float(count)
        return -log(norm.pdf(float(count), 89., 8.))
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

        angle_stats = AngleStatsDict(ConstructionConfig.angle_stats_file)
        stem_stats = StemStatsDict(ConstructionConfig.angle_stats_file)
        loop_stats = LoopStatsDict(ConstructionConfig.angle_stats_file)

        sm = SpatialModel(bg, angle_stats, stem_stats, loop_stats)
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


    def eval_energy(self, bg, background=True):
        '''
        Evaluate the energy of a coarse-grained structure.

        @param bg: The representation of the coarse-grained RNA molecule.
        '''
        energy = 0.

        for interaction in lri_iter.iterate_over_interactions(bg):
            try:
                #print >>stderr, "here", self.energies[interaction]
                energy += self.energies[interaction]
            except KeyError:
                #print >>stderr, "not here", 1 / self.calibration_size
                energy += 1 / float(self.calibration_size)

        return -log(energy)

    def naive_energy(self, bg):
        keys = bg.defines.keys()
        energy = 0.

        for i in range(len(keys)):
            for j in range(i+1, len(keys)):
                energy += 1

        return -energy
