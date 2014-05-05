#!/usr/bin/python

import collections as c
import sys, random, copy
import numpy as np
import math, os
import time

import scipy.stats as ss

#import matplotlib.pyplot as plt

import fess.builder.config as cbc
import forgi.threedee.model.stats as cbs

import forgi.threedee.utilities.graph_pdb as ftug
import forgi.threedee.utilities.vector as ftuv
import forgi.utilities.debug as fud

import forgi.threedee.utilities.rmsd as cbr
import fess.builder.models as cbm

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
        kernel = ss.gaussian_kde(values)
        Z = np.reshape(kernel(positions).T, X.shape)

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
                        s = sorted_energy_rmsds[: 3 * len(sorted_energy_rmsds) / 4]
                        s = random.sample(s, min(len(s), 180))

                        e = [s1[0] for s1 in s if s1[2] == color]
                        r = [s1[1] for s1 in s if s1[2] == color]

                        #self.create_contour_plot(np.array(r), np.array(e), self.ax_plot, xlim, ylim, color)
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

    def __init__(self, sm_orig, plotter=None, plot_color=None, silent=False, output_file=sys.stdout, save_n_best=3, dist1=None, dist2=None, no_rmsd=False, save_iterative_cg_measures=False):
        '''
        @param sm_orig: The original Spatial Model against which to collect statistics.
        '''
        self.energy_rmsd_structs = []
        self.counter = 0
        self.plotter = plotter
        self.plot_color = plot_color
        self.silent = silent
        self.verbose = False
        self.output_file = output_file
        self.save_n_best = save_n_best
        self.sm_orig = sm_orig
        self.energy_orig = None
        self.step_save = 0
        self.save_iterative_cg_measures=save_iterative_cg_measures

        self.dist1 = dist1
        self.dist2 = dist2

        self.highest_rmsd = 0.
        self.lowest_rmsd = 10000000000.
        self.no_rmsd = no_rmsd
        self.creation_time = time.time()

        try:
            self.centers_orig = ftug.bg_virtual_residues(sm_orig.bg)
        except KeyError:
            # if there are no coordinates provided in the original
            # bulge graph file, then don't calculate rmsds
            self.centers_orig = None

    def update_statistics(self, energy_function, sm, prev_energy, tracking_energies = []):
        '''
        Add a newly sampled structure to the set of statistics.

        @param energy_function: The energy_function used to evaluate the structure.
        @param sm: The spatial model that was sampled.
        '''
        self.counter += 1
        #sm.traverse_and_build()

        if self.energy_orig == None:
            self.energy_orig = 0.
            try:
                for s in sm.bg.stem_iterator():
                    ftug.add_virtual_residues(self.sm_orig.bg, s)
                self.energy_orig = energy_function.eval_energy(self.sm_orig)
            except KeyError:
                # most likely no native structure was provided
                pass

        energy = prev_energy #energy_function.eval_energy(sm, background=True)
        #energy = self.sampled_energy
        if self.sampled_energy != energy:
            pass
            #fud.pv('self.sampled_energy, energy')

        if self.centers_orig != None:
            # no original coordinates provided so we can't calculate rmsds
            r = 0.
            if not self.no_rmsd:
                centers_new = ftug.bg_virtual_residues(sm.bg)
                #r = cbr.centered_rmsd(self.centers_orig, centers_new)
                r = cbr.drmsd(self.centers_orig, centers_new)
        else:
            r = 0.

        dist = None
        if self.dist1 and self.dist2:
            if self.dist1 < 0:
                print >> sys.stderr, "The first distance nucleotide number should be greater than or equal to 0."
                sys.exit(1)
            elif self.dist1 > sm.bg.seq_length:
                print >> sys.stderr, "The first distance nucleotide number should be less than the length of the molecule."
                sys.exit(1)
            if self.dist2 < 0:
                print >> sys.stderr, "The second distance nucleotide number should be greater than or equal to 0."
                sys.exit(1)
            elif self.dist2 > sm.bg.seq_length:
                print >> sys.stderr, "The second distance nucleotide number should be less than the length of the molecule."
                sys.exit(1)


            atoms = ftug.virtual_atoms(sm.bg)

            d1 = sm.bg.get_node_from_residue_num(self.dist1)
            d2 = sm.bg.get_node_from_residue_num(self.dist2)

            #fud.pv('d1,d2')
            #fud.pv('atoms[self.dist1]')

            dist = ftuv.vec_distance(atoms[self.dist1]["C1'"],
                                     atoms[self.dist2]["C1'"])


        #self.energy_rmsd_structs += [(energy, r, sm.bg)]
        self.energy_rmsd_structs += [(energy, r, copy.deepcopy(sm.bg))]
        #self.energy_rmsd_structs += [(energy, r, sm.bg.copy())]

        sorted_energies = sorted(self.energy_rmsd_structs, key=lambda x: x[0])
        self.energy_rmsd_structs = sorted_energies[:self.save_n_best]

        if r > self.highest_rmsd:
            self.highest_rmsd = r

        if r < self.lowest_rmsd:
            self.lowest_rmsd = r

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

            output_str = "native_energy [%s %d]: %3d %5.03g  %5.3f | min: %5.2f (%5.2f) %5.2f | extreme_rmsds: %5.2f %5.2f" % ( sm.bg.name, sm.bg.seq_length, self.counter, energy, r , lowest_energy, self.energy_orig, lowest_rmsd, self.lowest_rmsd, self.highest_rmsd)
            output_str += " |"

            '''
            for e in tracking_energies[:1]:
                output_str += " %.2f" % (e.prev_cg)
            '''

            if dist:
                output_str += " | dist %.2f" % (dist)

            output_str += " [time: %.1f]" % (time.time() - self.creation_time)
            output_str += "\n"

            if self.output_file != sys.stdout:
                print output_str.strip()

            if self.output_file != None:
                self.output_file.write(output_str)
                self.output_file.flush()

        self.update_plots(energy, r)

        if self.counter % 10 == 0:
            if not self.silent:
                self.save_top(self.save_n_best)

        if self.step_save > 0:
            if self.counter % self.step_save == 0:
                sm.bg.to_cg_file(os.path.join(cbc.Configuration.sampling_output_dir, 'step%06d.coord' % (self.counter)))
            

    def save_top(self, n = 100000):
        '''
        Save the top n structures.
        '''
        if n > len(self.energy_rmsd_structs):
            n = len(self.energy_rmsd_structs)

        sorted_energies = sorted(self.energy_rmsd_structs, key=lambda x: x[0])
        '''
        if self.save_iterative_cg_measures:
            self.energy_function.dump_measures(cbc.Configuration.sampling_output_dir, self.counter)
        else:
            self.energy_function.dump_measures(cbc.Configuration.sampling_output_dir)
        '''

        #self.energy_function.resample_background_kde(sorted_energies[0][2])

        for i in range(n):
            sorted_energies[i][2].to_cg_file(os.path.join(cbc.Configuration.sampling_output_dir, 'best%d.coord' % (i)))

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
    def __init__(self, sm, energy_function, stats, stats_type='discrete', no_rmsd=False, energies_to_track=[]):
        '''
        param @sm: SpatialModel that will be used for sampling.
        '''
        if stats_type == 'continuous':
            self.cont_stats = cbs.ContinuousAngleStats(cbs.get_angle_stats())
        elif stats_type == 'random':
            self.cont_stats = cbs.RandomAngleStats(cbs.get_angle_stats())
        else:
            self.cont_stats = None
        
        self.step_counter = 0
        self.sm = sm
        self.energy_function = energy_function
        self.stats = stats
        self.stats.energy_function = energy_function
        self.prev_energy = 100000000000.
        self.energies_to_track = energies_to_track
        self.dump_measures = False

        sm.sample_stats()
        constraint_energy = sm.constraint_energy
        if sm.constraint_energy is not None:
            junction_constraint_energy = sm.junction_constraint_energy
            sm.constraint_energy = None
            sm.junction_constraint_energy = None
            sm.traverse_and_build()
            sm.constraint_energy = constraint_energy
            sm.junction_constraint_energy = junction_constraint_energy
            sm.traverse_and_build()
            energy_function.energies += sm.constraint_energy.energies
            energy_function.energies += [sm.junction_constraint_energy]

        sm.constraint_energy = None
        sm.junction_constraint_energy = None
        self.no_rmsd = no_rmsd

        #fud.pv('constraint_energy.eval_energy(sm)')
        #fud.pv('junction_constraint_energy.eval_energy(sm)')
        sm.traverse_and_build()
        #fud.pv('constraint_energy.eval_energy(sm)')
        #fud.pv('junction_constraint_energy.eval_energy(sm)')
        #fud.pv('energy_function.eval_energy(sm)')
        self.prev_energy = energy_function.eval_energy(sm)
        #sys.exit(1)
        sm.get_sampled_bulges()

    def get_random_bulge(self, bg):
        '''
        Return a random interior or multiloop from the structure.
        
        @param bg: A forgi.graph.BulgeGraph data structure.
        @return: The name of one of the bulges in bg selected randomly
        '''
        bulges = []
        for d in bg.defines.keys():
            if d[0] == 'i' or d[0] == 'm':
                bulges += [d]
        return random.choice(bulges)

    def get_random_stem(self, bg):
        '''
        Return a random stem from the bulge graph.
        
        @param bg: A forgi.graph.BulgeGraph data structure.
        @return: The name of one of the stems in bg selected randomly
        '''
        bulges = []
        for d in bg.defines.keys():
            if d[0] == 's':
                bulges += [d]
        return random.choice(bulges)

    def get_random_loop(self, bg):
        '''
        Return a random hairpin loop from the bulge graph.
        
        @param bg: A forgi.graph.BulgeGraph data structure.
        @return: The name of one of the hairpins in bg selected randomly
        '''
        bulges = []
        for d in bg.defines.keys():
            if d[0] == 'h':
                bulges += [d]

        return random.choice(bulges)

    def change_elem(self):
        '''
        Change a random element and accept the change with a probability
        proportional to the energy function.
        '''
        # pick a random element and get a new statistic for it
        d = random.choice(list(self.sm.bg.get_mst()))
        new_stat = random.choice(self.sm.conf_stats.sample_stats(self.sm.bg, d))

        # get the energy before we replace the statistic
        # it's dubious whether this is really necessary since we already
        # store the previous energy in the accept/reject step
        self.prev_energy = self.energy_function.eval_energy(self.sm, 
                                                            background=True)
         
        prev_stat = self.sm.elem_defs[d]

        # replace the statistic, rebuild the struture 
        # and calculate a new energy
        self.sm.elem_defs[d] = new_stat
        self.sm.traverse_and_build(start=d)
        energy = self.energy_function.eval_energy(self.sm, background=True)
        self.stats.sampled_energy = energy

        if energy < self.prev_energy:
            # lower energy means automatic acceptance accordint to the
            # metropolis hastings criterion
            self.prev_energy = energy
            self.energy_function.accept_last_measure()
        else:
            # calculate a probability
            r = random.random()
            if r > math.exp(self.prev_energy - energy):
                # reject the sampled statistic and replace it the old one
                self.sm.elem_defs[d] = prev_stat
                self.sm.traverse_and_build(start=d)
                self.stats.sampled_energy = self.prev_energy
                self.energy_function.reject_last_measure()
            else:
                # accept the new statistic
                self.prev_energy = energy
                self.energy_function.accept_last_measure


    def change_angle(self):
        # pick a random bulge to vary
        bulge = self.get_random_bulge(self.sm.bg)

        if bulge == None:
            # this structure has no bulges
            return

        while bulge not in self.sm.sampled_bulges:
            bulge = self.get_random_bulge(self.sm.bg)

        dims = self.sm.bg.get_bulge_dimensions(bulge)

        # What are the potential angle statistics for it
        #possible_angles = cbs.get_angle_stats()[dims[0]][dims[1]]
        (bulge, ang_type1) = random.choice(self.sm.sampled_bulge_sides)
        dims = self.sm.bg.get_bulge_dimensions(bulge)
        
        #self.sm.traverse_and_build()
        # What are the potential angle statistics for it
        try:
            if self.cont_stats:
                pa = self.cont_stats.sample_stats(tuple(list(dims) + [ang_type1]))
            else:
                possible_angles = cbs.get_angle_stats()[(dims[0], dims[1], ang_type1)]
                pa = random.choice(possible_angles)
        except IndexError:
            (dist, size1, size2, typ1) = cbs.get_angle_stat_dims(dims[0], dims[1], ang_type1)[0]
            if self.cont_stats:
                pa = self.cont_stats.sample_stats(tuple(list(dims) + [ang_type1]))
            else:
                possible_angles = cbs.get_angle_stats()[(size1, size2, ang_type1)]
                pa = random.choice(possible_angles)

        #fud.pv('"pe", self.energy_function.eval_energy(self.sm, background=True)')
        self.prev_energy = self.energy_function.eval_energy(self.sm, background=True)
        prev_angle = self.sm.angle_defs[bulge][ang_type1]
        self.sm.angle_defs[bulge][ang_type1] = pa
        self.sm.traverse_and_build(start=bulge)
        
        energy = self.energy_function.eval_energy(self.sm, background=True)
        #fud.pv('self.prev_energy, energy')
        self.stats.sampled_energy = energy
        if energy > self.prev_energy:
            r = random.random()
            #fud.pv('self.prev_energy, energy')
            #print >>sys.stderr, "checking r:", r
            if r > math.exp(self.prev_energy - energy):
                #fud.pv('self.sm.angle_defs[bulge][ang_type1]')
                #fud.pv('prev_angle')
                self.sm.angle_defs[bulge][ang_type1] = prev_angle
                self.sm.traverse_and_build(start=bulge)
                self.stats.sampled_energy = self.prev_energy
                #fud.pv('self.energy_function.eval_energy(self.sm, background=True)')
                #print >>sys.stderr, "rejecting r:", r
                #fud.pv('self.prev_energy, energy, self.energy_function.uncalibrated_energies[-1].measures[-1:]')
                self.energy_function.reject_last_measure()
            else:
                self.prev_energy = energy
                self.energy_function.accept_last_measure()
        else:
            self.prev_energy = energy
            self.energy_function.accept_last_measure()


    def change_stem(self):
        # pick a random bulge to vary
        stem = self.get_random_stem(self.sm.bg)
        length = self.sm.bg.get_length(stem)

        # What are the potential angle statistics for it
        #possible_angles = cbs.get_angle_stats()[dims[0]][dims[1]]
        
        #self.sm.traverse_and_build()
        # What are the potential angle statistics for it
        possible_stems = cbs.get_stem_stats()[length]

        pa = random.choice(possible_stems)

        #fud.pv('"pe", self.energy_function.eval_energy(self.sm, background=True)')
        self.prev_energy = self.energy_function.eval_energy(self.sm, background=True)
        prev_stem = self.sm.stem_defs[stem]
        #fud.pv('loop,length,self.sm.loop_defs[loop]')
        self.sm.stem_defs[stem] = pa
        self.sm.traverse_and_build(start=stem)

        energy = self.energy_function.eval_energy(self.sm, background=True)
        self.stats.sampled_energy = energy
        if energy > self.prev_energy:
            if random.random() > math.exp(self.prev_energy - energy):
                self.sm.stem_defs[stem] = prev_stem
                self.sm.traverse_and_build()
                self.stats.sampled_energy = self.prev_energy
                #fud.pv('self.energy_function.eval_energy(self.sm, background=True)')
                self.energy_function.reject_last_measure()
            else:
                self.prev_energy = energy
                self.energy_function.accept_last_measure()
        else:
            self.prev_energy = energy
            self.energy_function.accept_last_measure()

    def step(self):
        #self.sm.sample_stems()
        #self.sm.sample_loops()
        #self.sm.sample_angles()
        #self.sm.traverse_and_build()

        self.change_elem()
        '''
        if random.random() < 0.5:
            self.change_angle()
        elif random.random() < 0.5:
            self.change_stem()
        else:
            self.change_loop()
        '''

        if self.dump_measures:
            if self.step_counter % 20 == 0:
                self.energy_function.dump_measures(cbc.Configuration.sampling_output_dir)

        if self.step_counter % 3 == 0:
            self.energy_function.resample_background_kde(self.sm.bg)

        #fud.pv('self.energy_function.uncalibrated_energies[-1].accepted_measures[-1]')
        self.step_counter += 1

        if self.step_counter % 3 == 0:
            self.energy_function.resample_background_kde(self.sm.bg)

        for e in self.energies_to_track:
            e.eval_energy(self.sm)
            e.accepted_measures += [e.measures[-1]]

            if self.step_counter % 20 == 0:
                e.dump_measures(cbc.Configuration.sampling_output_dir)

        self.stats.update_statistics(self.energy_function, self.sm, self.prev_energy, self.energies_to_track)

    def change_loop(self):
        # pick a random bulge to vary
        loop = self.get_random_loop(self.sm.bg)
        length = self.sm.bg.get_length(loop)

        # What are the potential angle statistics for it
        #possible_angles = cbs.get_angle_stats()[dims[0]][dims[1]]
        
        #self.sm.traverse_and_build()
        # What are the potential angle statistics for it
        possible_loops = cbs.get_loop_stats()[length]
        #self.sm.traverse_and_build()


        pa = random.choice(possible_loops)
        #fud.pv('"pe", self.energy_function.eval_energy(self.sm, background=True)')
        self.prev_energy = self.energy_function.eval_energy(self.sm, background=True)
        prev_loop = self.sm.loop_defs[loop]
        #fud.pv('loop,length,self.sm.loop_defs[loop]')
        self.sm.loop_defs[loop] = pa
        self.sm.traverse_and_build(start=loop)

        energy = self.energy_function.eval_energy(self.sm, background=True)
        self.stats.sampled_energy = energy
        if energy > self.prev_energy:
            if random.random() > math.exp(self.prev_energy - energy):
                self.sm.loop_defs[loop] = prev_loop
                self.sm.traverse_and_build(start=loop)
                self.stats.sampled_energy = self.prev_energy
                #fud.pv('self.energy_function.eval_energy(self.sm, background=True)')
                self.energy_function.reject_last_measure()
            else:
                self.prev_energy = energy
                self.energy_function.accept_last_measure()
        else:
            self.prev_energy = energy
            self.energy_function.accept_last_measure()

    '''
    def step(self):
        #self.sm.sample_stems()
        #self.sm.sample_loops()
        #self.sm.sample_angles()
        self.sm.traverse_and_build()

        if random.random() < 0.5:
            self.change_angle()
        else:
            self.change_loop()

        self.stats.update_statistics(self.energy_function, self.sm)
    '''

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
        self.angles_to_sample = 50

        sm.sample_stats()
        sm.get_sampled_bulges()

    def step(self):
        '''
        Perform another step in the simulation.
        '''

        self.sm.sample_stems()
        self.sm.sample_loops()
        self.sm.traverse_and_build()

        # pick a random bulge to vary
        #bulge = self.sm.bg.get_random_bulge()
        (bulge, ang_type1) = random.choice(self.sm.sampled_bulge_sides)
        dims = self.sm.bg.get_bulge_dimensions(bulge)
        
        # What are the potential angle statistics for it
        (dist, size1, size2, type1) = cbs.get_angle_stat_dims(dims[0], dims[1], ang_type1)[0]
        possible_angles = cbs.get_angle_stats()[(size1, size2, ang_type1)]

        if len(possible_angles) == 0:
            print >>sys.stderr, "s1b", s1b, "s2b", s2b
        # only choose 10 possible angles
        if len(possible_angles) > self.angles_to_sample:
            possible_angles = random.sample(possible_angles, self.angles_to_sample)

        #possible_angles += [self.sm.angle_defs[bulge][ang_type1]]
        energies = dict()

        # evaluate the energies of the structure when the original
        # angle is replaced by one of the 10 potential new ones
        for pa in possible_angles:
            self.sm.angle_defs[bulge][ang_type1] = pa
            self.sm.traverse_and_build(start=bulge)
            energy = self.energy_function.eval_energy(self.sm, background=True)
            energies[pa] = energy
        #fud.pv('[v for v in energies.values()]')


        # energy = -log(p(S)) 
        # So... we want to maximize p(S)
        # Therefore, we want to minimize the energy
        max_energy = max(energies.values())
        min_energy = min(energies.values())

        if max_energy - min_energy > 40:
            max_energy = min_energy + 40.

        #fud.pv('max_energy')
        #fud.pv('min_energy')
        for pa in possible_angles:
            prev_energy = energies[pa]
            if prev_energy > max_energy:
                prev_energy = max_energy 
            
            prev_energy = prev_energy - min_energy
            energies[pa] = np.exp(-prev_energy)

        # Convert all of the sampled energies into one probability
        total_energy = sum([energies[key] for key in energies.keys()])

        energy_probs = dict()
        for key in energies.keys():
            energy_probs[key] = energies[key] / total_energy
        #fud.pv('[v for v in energy_probs.values()]')

        # sanity check
        total_prob = sum([energy_probs[key] for key in energies.keys()])
        #assert(allclose(total_prob, 1.))

        #pick one new angle to accept given the probabilities of the
        #sampled ones
        prob_remaining = 1.
        for key in energy_probs.keys():
            if random.random() < energy_probs[key] / prob_remaining:
                #fud.pv('energy_probs[key]')
                #fud.pv('energies[key]')
                self.sm.angle_defs[bulge][ang_type1] = key
                #self.sm.traverse_and_build(start=bulge)
                #fud.pv('self.energy_function.eval_energy(self.sm, background=True)')
                #fud.pv('self.energy_function.eval_energy(self.sm, background=True)')
                break

            prob_remaining -= energy_probs[key]

        self.sm.traverse_and_build(start=bulge)
        self.stats.update_statistics(self.energy_function, self.sm)

