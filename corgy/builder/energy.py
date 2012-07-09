#!/usr/bin/python

import pdb

from corgy.utilities.vector import vec_distance
from bobbins_config import ConstructionConfig

from corgy.builder.models import SpatialModel
from corgy.builder.stats import AngleStatsDict, StemStatsDict, LoopStatsDict

from corgy.utilities.data_structures import DefaultDict
from corgy.utilities.statistics import fit_skew, skew

from pylab import plot,show, hist

from scipy.stats import norm
from math import log

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

    def eval_energy(self, bg):
        total_energy = 0.

        for energy in self.energies:
            total_energy += energy.eval_energy(bg)

        return total_energy

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

    def eval_energy(self, bg):
        count = self.count_interactions(bg)

        #return float(count)
        return -log(norm.pdf(float(count), 89., 8.))
        

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


    def eval_energy(self, bg):
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
