#!/usr/bin/python

import collections as c
import sys, random, copy
import numpy as np
import scipy.stats as ss
import math, os

#import matplotlib.pyplot as plt

import corgy.builder.config as cbc
import corgy.builder.stats as cbs
import corgy.graph.graph_pdb as cgg
import corgy.builder.rmsd as cbr
import corgy.builder.models as cbm
import corgy.utilities.debug as cud

import numpy as np

class StatisticsPlotter:
    '''
    Plot a histogram of the rmsd as well as a plot of the energy vs. the rmsd.
    '''
    def __init__(self):
        import matplotlib.pyplot as plt
        import pylab as pl
        pl.ion()

        self.fig = plt.figure(figsize=(9, 9))

        self.ax_hist = self.fig.add_subplot(2, 1, 1)
        self.ax_plot = self.fig.add_subplot(2, 1, 2, sharex=self.ax_hist)

        self.energies = c.defaultdict(list)
        self.rmsds = c.defaultdict(list)
        self.energy_rmsds = []

        self.ax_plot.set_xlabel('rmsd')
        self.ax_hist.set_xlabel('rmsd')
        
        plt.ion()

    def create_contour_plot(self, m1, m2, ax, xlim, ylim, color):
        import matplotlib.pyplot as plt
        new_m1 = []
        new_m2 = []

        for i in range(len(m1)):
            if m1[i] > xlim[0] and m1[i] < xlim[1] and m2[i] > ylim[0] and m2[i] < ylim[1]:
                new_m1 += [m1[i]]
                new_m2 += [m2[i]]

        #positions = np.vstack([X.ravel(), Y.ravel()])
        X, Y = np.mgrid[xlim[0]:xlim[1]:complex(0,len(new_m1)), ylim[0]:ylim[1]:complex(0, len(new_m2))]
        #X,Y = np.meshgrid(new_m1, new_m2)
        positions = np.vstack([X.ravel(), Y.ravel()])
        values = np.vstack([m1, m2])
        #print "values:", values
        kernel = ss.gaussian_kde(values)
        Z = np.reshape(kernel(positions).T, X.shape)

        #print "Z:", Z

        if color == 'b':
            ax.contourf(X, Y, Z, cmap=plt.cm.Blues,alpha=0.5)
        if color == 'r':
            ax.contourf(X, Y, Z, cmap=plt.cm.Reds,alpha=0.5)
        if color == 'g':
            ax.contourf(X, Y, Z, cmap=plt.cm.Greens, alpha=0.5)
        if color == 'y':
            ax.contourf(X, Y, Z, cmap=plt.cm.YlOrBr, alpha=0.5)

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
                #xlim = (sorted_rmsds[0] - 5., sorted_rmsds[3 * len(sorted_rmsds) / 4] + 5.)
                xlim = (0, sorted_rmsds[-1] + 0.5)

                self.xlim = xlim
                self.ylim = ylim

                self.ax_plot.set_ylim(ylim)
                #self.ax_plot.set_xlim(xlim)

            for i in range(min(5, len(sorted_energy_rmsds))):
                self.ax_plot.plot(sorted_energy_rmsds[i][1], sorted_energy_rmsds[i][0], '%so' % (sorted_energy_rmsds[i][2]), alpha=0.5)


            if len(self.energies[color]) > 2. and len(sorted_energies) > 4:
                for color in self.energies.keys():
                    try: 
                        #s = random.sample(sorted_energy_rmsds, min(len(sorted_energy_rmsds), 180))
                        #cud.pv('s')
                        s = sorted_energy_rmsds[: 3 * len(sorted_energy_rmsds) / 4]
                        s = random.sample(s, min(len(s), 180))

                        #cud.pv('s2')
                        e = [s1[0] for s1 in s if s1[2] == color]
                        r = [s1[1] for s1 in s if s1[2] == color]

                        self.create_contour_plot(np.array(r), np.array(e), self.ax_plot, xlim, ylim, color)
                    except Exception as ex:
                        print "exception:", ex, "color:", color

                        continue
                    
            for color in self.energies.keys():
                self.ax_plot.plot(self.rmsds[color], self.energies[color], '%so' % (color), alpha=0.05)

        for color in self.energies.keys():
            xlim = (0, sorted_rmsds[-1] + 0.5)
            self.ax_hist.set_xlim(xlim)
            self.ax_hist.hist(self.rmsds[color], color=color, alpha=0.5, normed=True)

        import matplotlib.pyplot as plt
        plt.draw()

    def diagnose_energy(self, energy_function, bgs):
        energies = [energy_function.eval_energy(cbm.SpatialModel(bg), background=True) for bg in bgs]


    def finish(self):
        self.ax_plot.cla()
        sorted_energy_rmsds = sorted(self.energy_rmsds)

        sorted_energies = sorted([s[0] for s in sorted_energy_rmsds])
        sorted_rmsds = sorted([s[1] for s in sorted_energy_rmsds])

        se = sorted_energies[:3 * len(sorted_energies) / 4]
        sr = sorted_rmsds[:3 * len(sorted_rmsds) / 4]

        ylim = (sorted_energies[0] - 5., sorted_energies[3 * len(sorted_energies) / 4] + 5.)
        #xlim = (sorted_rmsds[0] - 5., sorted_rmsds[3 * len(sorted_rmsds) / 4] + 5.)
        xlim = (0, sorted_rmsds[3 * len(sorted_rmsds) / 4] + 5.)

        self.xlim = xlim
        self.ylim = ylim

        self.ax_plot.set_ylim(ylim)
        self.ax_plot.set_xlim(xlim)

        for i in range(min(5, len(sorted_energy_rmsds))):
            self.ax_plot.plot(sorted_energy_rmsds[i][1], sorted_energy_rmsds[i][0], '%so' % (sorted_energy_rmsds[i][2]), alpha=0.5)
                    
        for color in self.energies.keys():
            self.ax_plot.plot(self.rmsds[color], self.energies[color], '%so' % (color), alpha=0.05)

        for color in self.energies.keys():
            self.create_contour_plot(np.array(self.rmsds[color]), np.array(self.energies[color]), self.ax_plot, self.xlim, self.ylim, color)

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
        self.counter = 0
        self.plotter = plotter
        self.plot_color = plot_color
        self.silent = silent
        self.verbose = False

        try:
            self.centers_orig = cgg.bg_virtual_residues(sm_orig.bg)
        except KeyError:
            # if there are no coordinates provided in the original
            # bulge graph file, then don't calculate rmsds
            self.centers_orig = None

    def update_statistics(self, energy_function, sm):
        '''
        Add a newly sampled structure to the set of statistics.

        @param energy_function: The energy_function used to evaluate the structure.
        @param sm: The spatial model that was sampled.
        '''
        self.counter += 1
        #sm.traverse_and_build()

        #energy = energy_function.eval_energy(sm.bg, background=True)
        energy = energy_function.eval_energy(sm, background=True)

        if self.centers_orig != None:
            # no original coordinates provided so we can't calculate rmsds
            centers_new = cgg.bg_virtual_residues(sm.bg)
            r = cbr.centered_rmsd(self.centers_orig, centers_new)
        else:
            r = 0.

        self.energy_rmsd_structs += [(energy, r, copy.deepcopy(sm.bg))]
        #self.energy_rmsd_structs += [(energy, r, sm.bg.copy())]

        sorted_energies = sorted(self.energy_rmsd_structs, key=lambda x: x[0])

        lowest_energy = sorted_energies[0][0]
        lowest_rmsd = sorted_energies[0][1]

        '''
        if energy == lowest_energy:
            for key in sm.angle_defs:
                print >>sys.stderr, key, str(sm.angle_defs[key])
        '''

        if not self.silent:
            if self.verbose:
                '''
                for energy_func in energy_function.energies:
                    print energy_func.__class__.__name__, energy_func.eval_energy(sm)
                '''

            print "native_energy [%s %d]: %3d %5.03g  %5.3f | min: %5.2f %5.2f" % ( sm.bg.name, sm.bg.length, self.counter, energy, r , lowest_energy, lowest_rmsd)
            sys.stdout.flush()

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
            sorted_energies[i][2].output(os.path.join(cbc.Configuration.sampling_output_dir, 'best%d.coord' % (i)))

    def update_plots(self, energy, rmsd):
        '''
        Maintain plots of all the necessary statistics.
        '''
        if self.plotter != None:
            self.plotter.add_data(energy, rmsd, self.plot_color)
    
    def print_final_stats(self, energy_function):
        sorted_energies = sorted(self.energy_rmsd_structs, key=lambda x: x[0])
        sm = cbm.SpatialModel(sorted_energies[0][2])
        sm.get_sampled_bulges()

        print "---------------------------"

        print [e[1] for e in sorted_energies[:10]]

        '''
        for energy in energy_function.energies:
            print energy.__class__.__name__, energy.eval_energy(sm)
        '''

        print "-------------------------"
        

class MCMCSampler:
    '''
    Sample using tradition accept/reject sampling.
    '''
    def __init__(self, sm, energy_function, stats):
        '''
        param @sm: SpatialModel that will be used for sampling.
        '''
        self.sm = sm
        self.energy_function = energy_function
        self.stats = stats
        self.prev_energy = 10000000.
        self.cont_stats = cbs.ContinuousAngleStats(cbs.get_angle_stats())

        sm.get_sampled_bulges()

    def step(self):
        self.sm.sample_stems()
        self.sm.sample_loops()
        self.sm.traverse_and_build()

        # pick a random bulge to vary
        bulge = self.sm.bg.get_random_bulge()
        while bulge not in self.sm.sampled_bulges:
            bulge = self.sm.bg.get_random_bulge()

        dims = self.sm.bg.get_bulge_dimensions(bulge)

        # What are the potential angle statistics for it
        #possible_angles = self.sm.angle_stats[dims[0]][dims[1]]
        #pa = random.choice(possible_angles)

        pa = self.cont_stats.sample_stats(dims)

        prev_angle = self.sm.angle_defs[bulge]
        self.sm.angle_defs[bulge] = pa
        self.sm.traverse_and_build(start = bulge)
        
        energy = self.energy_function.eval_energy(self.sm, background=True)
        if energy > self.prev_energy:
            #cud.pv('self.prev_energy')
            #cud.pv('energy')
            #cud.pv('math.exp(self.prev_energy - energy)')

            if random.random() > math.exp(self.prev_energy - energy):
                self.sm.angle_defs[bulge] = prev_angle
                #print "skipping:", energy, self.prev_energy
            else:
                self.prev_energy = energy
                self.stats.update_statistics(self.energy_function, self.sm)
        else:
            self.prev_energy = energy
            self.stats.update_statistics(self.energy_function, self.sm)

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

        sm.sample_stats()
        sm.get_sampled_bulges()
        #print >>stderr, "original native_energy:", energy_function.eval_energy(sm, background=True)

    def step(self):
        '''
        Perform another step in the simulation.
        '''

        self.sm.sample_stems()
        self.sm.sample_loops()
        self.sm.traverse_and_build()

        # pick a random bulge to vary
        #bulge = self.sm.bg.get_random_bulge()
        (bulge, (s1b, s2b)) = random.choice(self.sm.sampled_bulge_sides)
        dims = self.sm.bg.get_bulge_dimensions(bulge)
        
        if dims == (2,7):
            dims = (3,7)
        # What are the potential angle statistics for it
        possible_angles = self.sm.angle_stats[dims[0]][dims[1]][s1b][s2b]
        #print >>sys.stderr, bulge, dims, len(possible_angles)

        if len(possible_angles) == 0:
            print >>sys.stderr, "No available statistics for bulge %s of size %s" % (bulge, str(dims))
            print >>sys.stderr, "s1b", s1b, "s2b", s2b
        # only choose 10 possible angles
        if len(possible_angles) > 20:
            possible_angles = random.sample(possible_angles, 20)

        energies = dict()

        # evaluate the energies of the structure when the original
        # angle is replaced by one of the 10 potential new ones
        #cud.pv('possible_angles')
        for pa in possible_angles:
            self.sm.angle_defs[bulge][s1b][s2b] = pa
            self.sm.traverse_and_build(start=bulge)
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
            energies[pa] = np.exp(-prev_energy)
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
            if random.random() < energy_probs[key] / prob_remaining:
                self.sm.angle_defs[bulge][s1b][s2b] = key
                break

            prob_remaining -= energy_probs[key]

        self.stats.update_statistics(self.energy_function, self.sm)

