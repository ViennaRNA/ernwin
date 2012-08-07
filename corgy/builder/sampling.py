#!/usr/bin/python

import matplotlib.pyplot as plt

from corgy.utilities.data_structures import DefaultDict
from corgy.builder.rmsd import centered_rmsd
from corgy.builder.models import SpatialModel

from sys import stderr
from random import random, sample
from copy import deepcopy

from numpy import exp, array
from scipy.stats import gaussian_kde
import numpy as np

class StatisticsPlotter:
    '''
    Plot a histogram of the rmsd as well as a plot of the energy vs. the rmsd.
    '''
    def __init__(self):
        self.fig = plt.figure(figsize=(9,9))

        self.ax_hist = self.fig.add_subplot(2,1,1)
        self.ax_plot = self.fig.add_subplot(2,1,2)

        self.energies = DefaultDict([])
        self.rmsds = DefaultDict([])
        self.energy_rmsds = []

        self.ax_plot.set_xlabel('rmsd')
        self.ax_hist.set_xlabel('rmsd')
        
        plt.ion()

    def create_contour_plot1(self, m1, m2, ax, xlim, ylim, color):
        new_m1 = []
        new_m2 = []

        for i in range(len(m1)):
            if m1[i] > xlim[0] and m1[i] < xlim[1] and m2[i] > ylim[0] and m2[i] < ylim[1]:
                new_m1 += [m1[i]]
                new_m2 += [m2[i]]

        xmin = max(xlim[0], min(new_m1))
        xmax = min(xlim[1], max(new_m1))

        ymin = max(ylim[0], min(new_m2))
        ymax = min(ylim[1], max(new_m2))

        X, Y = np.mgrid[xmin:xmax:100j, ymin:ymax:100j]
        positions = np.vstack([X.ravel(), Y.ravel()])
        values = np.vstack([m1, m2])
        #print "values:", values
        kernel = gaussian_kde(values)
        Z = np.reshape(kernel(positions).T, X.shape)


        if color == 'b':
            ax.imshow(np.rot90(Z), cmap=plt.cm.Blues, extent=[xmin, xmax, ymin, ymax], alpha=0.5, aspect='auto')
        elif color == 'r':
            ax.imshow(np.rot90(Z), cmap=plt.cm.Reds, extent=[xmin, xmax, ymin, ymax], alpha=0.5, aspect='auto')

        #ax.plot(m1, m2, 'k.', markersize=2)
        #ax.set_xlim([xmin, xmax])
        #ax.set_ylim([ymin, ymax])

    def create_contour_plot(self, m1, m2, ax, xlim, ylim, color):
        new_m1 = []
        new_m2 = []

        for i in range(len(m1)):
            if m1[i] > xlim[0] and m1[i] < xlim[1] and m2[i] > ylim[0] and m2[i] < ylim[1]:
                new_m1 += [m1[i]]
                new_m2 += [m2[i]]

        xmin = max(xlim[0], min(new_m1))
        xmax = min(xlim[1], max(new_m1))

        ymin = max(ylim[0], min(new_m2))
        ymax = min(ylim[1], max(new_m2))

        #positions = np.vstack([X.ravel(), Y.ravel()])
        X, Y = np.mgrid[xlim[0]:xlim[1]:complex(0,len(new_m1)), ylim[0]:ylim[1]:complex(0,len(new_m2))]
        #X,Y = np.meshgrid(new_m1, new_m2)
        positions = np.vstack([X.ravel(), Y.ravel()])
        values = np.vstack([m1, m2])
        #print "values:", values
        kernel = gaussian_kde(values)
        Z = np.reshape(kernel(positions).T, X.shape)

        #print "Z:", Z

        if color == 'b':
            ax.contourf(X, Y, Z, cmap=plt.cm.Blues,alpha=0.5)
        if color == 'r':
            ax.contourf(X, Y, Z, cmap=plt.cm.Reds,alpha=0.5)

    def add_data(self, energy, rmsd, color):
        self.energies[color] += [energy]
        self.rmsds[color] += [rmsd] 

        self.energy_rmsds += [(energy, rmsd, color)]

        #print >>stderr, "energy", "rmsd", color

        all_energies = []
        all_rmsds = []

        for color in self.energies.keys():
            all_energies += list(self.energies[color])
            all_rmsds += list(self.rmsds[color])

        sorted_energy_rmsds = sorted(self.energy_rmsds)

        sorted_energies = sorted(all_energies)
        sorted_rmsds = sorted(all_rmsds)

        self.ax_hist.cla()

        if len(sorted_energies) % 2 == 0:
            self.ax_plot.cla()

            if len(sorted_energies) > 4:
                ylim = (sorted_energies[0] - 5., sorted_energies[3 * len(sorted_energies) / 4] + 5.)
                xlim = (sorted_rmsds[0] - 5., sorted_rmsds[3 * len(sorted_rmsds) / 4] + 5.)

                self.xlim = xlim
                self.ylim = ylim

                self.ax_plot.set_ylim(ylim)
                self.ax_plot.set_xlim(xlim)

            for i in range(min(5, len(sorted_energy_rmsds))):
                self.ax_plot.plot(sorted_energy_rmsds[i][1], sorted_energy_rmsds[i][0], '%so' % (sorted_energy_rmsds[i][2]), alpha=0.5)


            if len(self.energies[color]) > 2. and len(sorted_energies) > 4:
                for color in self.energies.keys():
                    try: 
                        s = sample(sorted_energy_rmsds, min(len(sorted_energy_rmsds), 180))
                        e = [s1[0] for s1 in s if s1[2] == color]
                        r = [s1[1] for s1 in s if s1[2] == color]

                        self.create_contour_plot(array(r), array(e), self.ax_plot, xlim, ylim, color)
                    except Exception as ex:
                        print "exception:", ex, "color:", color

                        continue
                    
            for color in self.energies.keys():
                self.ax_plot.plot(self.rmsds[color], self.energies[color], '%so' % (color), alpha=0.05)

        for color in self.energies.keys():
            self.ax_hist.hist(self.rmsds[color], color=color, alpha=0.5, normed=True)
        plt.draw()

    def diagnose_energy(self, energy_function, bgs):

        energies = [energy_function.eval_energy(SpatialModel(bg), background=True) for bg in bgs]


    def finish(self):
        self.ax_plot.cla()
        sorted_energy_rmsds = sorted(self.energy_rmsds)

        for i in range(min(5, len(sorted_energy_rmsds))):
            self.ax_plot.plot(sorted_energy_rmsds[i][1], sorted_energy_rmsds[i][0], '%so' % (sorted_energy_rmsds[i][2]), alpha=0.5)
                    
        for color in self.energies.keys():
            self.ax_plot.plot(self.rmsds[color], self.energies[color], '%so' % (color), alpha=0.05)

        for color in self.energies.keys():
            self.create_contour_plot(array(self.rmsds[color]), array(self.energies[color]), self.ax_plot, self.xlim, self.ylim, color)

        plt.ioff()
        plt.show()


class SamplingStatistics:
    '''
    Store statistics about a sample.
    '''

    def __init__(self, sm_orig, plotter=None, plot_color=None, silent=False):
        '''
        @param sm_orig: The original Spatial Model against which to collect statistics.
        '''
        self.energy_rmsd_structs = []
        self.centers_orig = sm_orig.bg.get_centers()
        self.counter = 0
        self.plotter = plotter
        self.plot_color = plot_color
        self.silent = silent
        self.verbose = False

    def update_statistics(self, energy_function, sm):
        '''
        Add a newly sampled structure to the set of statistics.

        @param energy_function: The energy_function used to evaluate the structure.
        @param sm: The spatial model that was sampled.
        '''
        self.counter += 1
        sm.traverse_and_build()
        #energy = energy_function.eval_energy(sm.bg, background=True)
        energy = energy_function.eval_energy(sm, background=True)

        centers_new = sm.bg.get_centers()
        r = centered_rmsd(self.centers_orig, centers_new)

        self.energy_rmsd_structs += [(energy, r, deepcopy(sm.bg))]

        sorted_energies = sorted(self.energy_rmsd_structs, key=lambda x: x[0])

        lowest_energy = sorted_energies[0][0]
        lowest_rmsd = sorted_energies[0][1]

        if not self.silent:
            if self.verbose:
                '''
                for energy_func in energy_function.energies:
                    print energy_func.__class__.__name__, energy_func.eval_energy(sm)
                '''

            print "native_energy: %3d %5.2f  %5.2f | min: %5.2f %5.2f" % ( self.counter, energy, r , lowest_energy, lowest_rmsd)

        self.update_plots(energy, r)

        if self.counter % 10 == 0:
            if not self.silent:
                self.save_top(10)

    def save_top(self, n = 100000):
        '''
        Save the top n structures.
        '''
        if n > len(self.energy_rmsd_structs):
            n = len(self.energy_rmsd_structs)

        sorted_energies = sorted(self.energy_rmsd_structs, key=lambda x: x[0])

        for i in range(n):
            sorted_energies[i][2].output('best/best%d.coord' % (i))

    def update_plots(self, energy, rmsd):
        '''
        Maintain plots of all the necessary statistics.
        '''
        if self.plotter != None:
            self.plotter.add_data(energy, rmsd, self.plot_color)
    
    def print_final_stats(self, energy_function):
        sorted_energies = sorted(self.energy_rmsd_structs, key=lambda x: x[0])
        sm = SpatialModel(sorted_energies[0][2])
        sm.get_sampled_bulges()

        print "---------------------------"

        '''
        for energy in energy_function.energies:
            print energy.__class__.__name__, energy.eval_energy(sm)
        '''

        print "-------------------------"
        


class GibbsBGSampler:
    '''
    A Gibbs Sampler for Bulge Graphs.
    '''

    def __init__(self, sm, energy_function, stats):
        '''
        param @sm: SpatialModel that will be used for sampling.
        '''
        self.sm = sm
        self.energy_function = energy_function
        self.stats = stats

        sm.get_sampled_bulges()

        #print >>stderr, "original native_energy:", energy_function.eval_energy(sm, background=True)

    def step(self):
        '''
        Perform another step in the simulation.
        '''

        self.sm.sample_stems()
        self.sm.sample_loops()

        # pick a random bulge to vary
        bulge = self.sm.bg.get_random_bulge()
        dims = self.sm.bg.get_bulge_dimensions(bulge)

        # What are the potential angle statistics for it
        possible_angles = self.sm.angle_stats[dims[0]][dims[1]]

        # only choose 10 possible angles
        if len(possible_angles) > 20:
            possible_angles = sample(possible_angles, 20)

        energies = dict()

        # evaluate the energies of the structure when the original
        # angle is replaced by one of the 10 potential new ones
        for pa in possible_angles:
            self.sm.angle_defs[bulge] = pa
            self.sm.traverse_and_build()
            energy = self.energy_function.eval_energy(self.sm, background=True)
            energies[pa] = energy


        # energy = -log(p(S)) 
        # So... we want to maximize p(S)
        # Therefore, we want to minimize the energy
        max_energy = max(energies.values())
        min_energy = min(energies.values())

        if max_energy - min_energy > 40:
            max_energy = min_energy + 40.

        #print >>stderr, "max_energy:", max_energy
        for pa in possible_angles:
            prev_energy = energies[pa]
            if prev_energy > max_energy:
                prev_energy = max_energy 
            
            prev_energy = prev_energy - min_energy
            energies[pa] = exp(-prev_energy)
            #print >>stderr, "energies[pa]:", energies[pa], "energy:", prev_energy

        # Convert all of the sampled energies into one probability
        total_energy = sum([energies[key] for key in energies.keys()])

        #print >>stderr, "total_energy:", total_energy
        energy_probs = dict()
        for key in energies.keys():
            energy_probs[key] = energies[key] / total_energy

        # sanity check
        total_prob = sum([energy_probs[key] for key in energies.keys()])
        #assert(allclose(total_prob, 1.))

        #pick one new angle to accept given the probabilities of the
        #sampled ones
        prob_remaining = 1.
        for key in energy_probs.keys():
            if random() < energy_probs[key] / prob_remaining:
                self.sm.angle_defs[bulge] = key
                break

            prob_remaining -= energy_probs[key]

        self.stats.update_statistics(self.energy_function, self.sm)

