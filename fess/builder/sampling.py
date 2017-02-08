#!/usr/bin/python
from __future__ import absolute_import, division, print_function, unicode_literals
from builtins import (ascii, bytes, chr, dict, filter, hex, input,
                      int, map, next, oct, open, pow, range, round,
                      str, super, zip)
from future.builtins.disabled import (apply, cmp, coerce, execfile,
                             file, long, raw_input, reduce, reload,
                             unicode, xrange, StandardError)


import collections as c
import sys, random, copy
import numpy as np
import math, os
import time
#import psutil
from . import config as conf

import scipy.stats as ss

#import matplotlib.pyplot as plt
from . import samplingStatisticsNew2 as sstats
from . import config as cbc
from . import energy as fbe
import forgi.threedee.model.similarity as ftme
import forgi.threedee.model.stats as ftms

import forgi.threedee.utilities.graph_pdb as ftug
import forgi.threedee.utilities.vector as ftuv
import forgi.utilities.debug as fud

import forgi.threedee.model.descriptors as cbr
import fess.builder.models as cbm
import logging
log = logging.getLogger(__name__)

#from guppy import hpy

def sizeof_fmt(num, suffix='B'):
    #http://stackoverflow.com/a/1094933
    for unit in ['','Ki','Mi','Gi','Ti','Pi','Ei','Zi']:
        if abs(num) < 1024.0:
            return "%3.1f%s%s" % (num, unit, suffix)
        num /= 1024.0
    return "%.1f%s%s" % (num, 'Yi', suffix)

'''def profile(fn):
    #http://stackoverflow.com/a/16624539
    def wrapper(*args, **kwargs):
        process = psutil.Process(os.getpid())
        start_rss, start_vms = process.memory_info()[:2]
        try:
            return fn(*args, **kwargs)
        finally:
            end_rss, end_vms = process.memory_info()[:2]
            print(fn.__doc__[1:25], sizeof_fmt(end_rss - start_rss), sizeof_fmt(end_vms - start_vms), sizeof_fmt(end_rss), sizeof_fmt(end_vms))
    return wrapper'''

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
        try:
          kernel = ss.gaussian_kde(values)
        except: pass
        else:
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
                self.ax_plot.plot(sorted_energy_rmsds[i][1], sorted_energy_rmsds[i][0], 'o', alpha=0.5) #'%so' % (sorted_energy_rmsds[i][2]), alpha=0.5)


            if len(self.energies[color]) > 2. and len(sorted_energies) > 4:
                for color in self.energies.keys():
                    try: 
                        #s = random.sample(sorted_energy_rmsds, min(len(sorted_energy_rmsds), 180))
                        s = sorted_energy_rmsds[: 3 * len(sorted_energy_rmsds) / 4]
                        s = random.sample(s, min(len(s), 180))

                        #self.create_contour_plot(np.array(r), np.array(e), self.ax_plot, xlim, ylim, color)
                    except Exception as ex:
                        print ("exception:", ex, "color:", color)

                        continue
                    
            for color in self.energies.keys():
                self.ax_plot.plot(self.rmsds[color], self.energies[color], 'o', alpha=0.5) #'%so' % (color), alpha=0.05)

        for color in self.energies.keys():
            xlim = (0, sorted_rmsds[-1] + 0.5)
            self.ax_hist.set_xlim(xlim)
            if len(self.rmsds[color])>2:
                self.ax_hist.hist(self.rmsds[color], color=color, alpha=0.5, normed=True)

        import matplotlib.pyplot as plt
        plt.draw()

    def diagnose_energy(self, energy_function, bgs):
        energies = [energy_function.eval_energy(cbm.SpatialModel(bg), background=True) for bg in bgs]

        return energies


    def finish(self):
        self.ax_plot.cla()
        sorted_energy_rmsds = sorted(self.energy_rmsds)

        sorted_energies = sorted([s[0] for s in sorted_energy_rmsds])
        sorted_rmsds = sorted([s[1] for s in sorted_energy_rmsds])

        se = sorted_energies[:3 * len(sorted_energies) // 4]
        sr = sorted_rmsds[:3 * len(sorted_rmsds) // 4]

        ylim = (sorted_energies[0] - 5., sorted_energies[3 * len(sorted_energies) // 4] + 5.)
        #xlim = (sorted_rmsds[0] - 5., sorted_rmsds[3 * len(sorted_rmsds) / 4] + 5.)
        xlim = (0, sorted_rmsds[3 * len(sorted_rmsds) // 4] + 5.)

        self.xlim = xlim
        self.ylim = ylim

        self.ax_plot.set_ylim(ylim)
        self.ax_plot.set_xlim(xlim)

        for i in range(min(5, len(sorted_energy_rmsds))):
            self.ax_plot.plot(sorted_energy_rmsds[i][1], sorted_energy_rmsds[i][0], 'o', alpha=0.5)
                    
        for color in self.energies.keys():
            self.ax_plot.plot(self.rmsds[color], self.energies[color], 'o', alpha=0.05)

        for color in self.energies.keys():
            self.create_contour_plot(np.array(self.rmsds[color]), np.array(self.energies[color]), self.ax_plot, self.xlim, self.ylim, color)

        import matplotlib.pyplot as plt
        plt.ioff()
        plt.show()


class SamplingStatistics:
    '''
    Store statistics about a sample.
    '''

    def __init__(self, sm_orig, plotter=None, plot_color=None, silent=False, output_file=sys.stdout, save_n_best=3, dists=[], no_rmsd=False, save_iterative_cg_measures=False):
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
        self.save_iterative_cg_measures = save_iterative_cg_measures



        self.dists = dists

        self.highest_rmsd = 0.
        self.lowest_rmsd = 10000000000.
        self.no_rmsd = no_rmsd
        self.creation_time = time.time()

        try:
            self.centers_orig = ftug.bg_virtual_residues(sm_orig.bg)        
            self.confusion_matrix_calculator = ftme.ConfusionMatrix(sm_orig.bg)
        except KeyError:
            # if there are no coordinates provided in the original
            # bulge graph file, then don't calculate rmsds
            self.centers_orig = None
            self.confusion_matrix_calculator = None


    def update_statistics(self, energy_function, sm, prev_energy, tracking_energies = None, tracked_energies=None):
        '''
        Add a newly sampled structure to the set of statistics.

        :param energy_function: The energy_function used to evaluate the structure.
        :param sm: The spatial model that was sampled.
        :param prev_energy: The evaluated (accepted) energy of the current step 
        :tracking_energyis: The energy_functions which are calculated for statistics, but not used for sampling.
        :tracked_energies: The energy values of the tracking_energies.
        '''
        self.counter += 1

        if self.energy_orig is None:
            self.energy_orig = 0.
            try:
                self.sm_orig.bg.add_all_virtual_residues()
                self.energy_orig = energy_function.eval_energy(self.sm_orig)
            except KeyError:
                # most likely no native structure was provided
                pass

        energy = prev_energy
        #energy = energy_function.eval_energy(sm, background=True)
        if energy_function.uses_background():
            energy_nobg = energy_function.eval_energy(sm, background=False)
        else:
            energy_nobg=energy

        mcc = None

        if self.centers_orig is not None:
            r = 0.
            if not self.no_rmsd:
                centers_new = ftug.bg_virtual_residues(sm.bg)
                r = cbr.centered_rmsd(self.centers_orig, centers_new)
                #r = cbr.drmsd(self.centers_orig, centers_new)
                cm = self.confusion_matrix_calculator.evaluate(sm.bg)
                mcc = ftme.mcc(cm)
        else:            
            # no original coordinates provided so we can't calculate rmsds
            r = 0.

        dist = None
        dist2 = None

        cg = sm.bg
        dists = []

        for (self.dist1, self.dist2) in self.dists:
            node1 = cg.get_node_from_residue_num(self.dist1)
            node2 = cg.get_node_from_residue_num(self.dist2)

            pos1, len1 = cg.get_position_in_element(self.dist1)
            pos2, len2 = cg.get_position_in_element(self.dist2)

            #fud.pv('node1, node2, pos1, pos2')

            vec1 = cg.coords[node1][1] - cg.coords[node1][0]
            vec2 = cg.coords[node2][1] - cg.coords[node2][0]

            #mid1 = (cg.coords[node1][0] + cg.coords[node1][1]) / 2
            #mid2 = (cg.coords[node2][0] + cg.coords[node2][1]) / 2

            mid1 = cg.coords[node1][0] + pos1 * (vec1 / len1)
            mid2 = cg.coords[node2][0] + pos2 * (vec2 / len2)
            
            #fud.pv('mid1, mid2')

            dists += [ftuv.vec_distance(mid1, mid2)]

        #self.energy_rmsd_structs += [(energy, r, sm.bg)]
        self.energy_rmsd_structs += [(energy_nobg, r, copy.deepcopy(sm.bg))]
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
            _, rog=fbe.length_and_rog(sm.bg)
            #output_str = u"native_energy [{:s} {:d}]: {:3d} {:5.03g} {:5.3f} ROG: {:5.3f} | min:
            output_str = u"native_energy [%s %d]: %3d %5.03g  %5.3f ROG: %5.3f | min: %5.2f (%5.2f) %5.2f | extreme_rmsds: %5.2f %5.2f (%.2f)" % ( sm.bg.name, sm.bg.seq_length, self.counter, energy, r , rog, lowest_energy, self.energy_orig, lowest_rmsd, self.lowest_rmsd, self.highest_rmsd, energy_nobg)
            output_str += " |"

            # assume that the energy function is a combined energy
            if isinstance(self.energy_function, fbe.CombinedEnergy):
                for e in self.energy_function.iterate_energies():
                    if isinstance(e,fbe.DistanceExponentialEnergy):
                        output_str += " [clamp {},{}: {:.1f}]".format(e.from_elem,
                                                                      e.to_elem,
                                                                      e.get_distance(sm))
            if tracked_energies and tracking_energies:
                output_str += " | [tracked Energies]"
                for i,e in enumerate(tracking_energies):
                    sn=e.shortname()
                    if len(sn)>12:
                        sn=sn[:9]+"..."
                    output_str += "  [{}]: ".format(sn)
                    output_str += "%5.03g" % (tracked_energies[i])
            elif tracking_energies:
                output_str += " | [tracked Energies]"
                for e in tracking_energies:
                    sn=e.shortname()
                    if len(sn)>12:
                        sn=sn[:9]+"..."
                    output_str += "  [{}]: ".format(sn)
                    output_str += "%5.03g" % (e.eval_energy(sm))

            if dist:
                output_str += " | dist %.2f" % (dist)

            for dist2 in dists:
                if dist2 is not None:
                    output_str += " | [dist2: %.2f]" % (dist2)

            if mcc is not None:
                output_str += " | [mcc: %.3f]" % (mcc)

            output_str += " [time: %.1f]" % (time.time() - self.creation_time)

            #Print to both STDOUT and the log file.
            if self.output_file != sys.stdout:
                print (output_str.strip())

            if self.output_file != None:
                print(output_str, file=self.output_file)
                self.output_file.flush()

        self.update_plots(energy, r)

        '''
        if self.counter % 1000 == 0:
            import pdb; pdb.set_trace()
        '''

        if self.counter % 10 == 0:
            if not self.silent:
                self.save_top(self.save_n_best, counter=self.counter)

        if self.step_save > 0 and self.counter % self.step_save == 0:
            #If a projection match energy was used, save the optimal projection direction to the file.
            if isinstance(self.energy_function, fbe.CombinedEnergy):
                for e in self.energy_function.iterate_energies():
                    if hasattr(e, "accepted_projDir"):
                        sm.bg.project_from=e.accepted_projDir
            sm.bg.to_cg_file(os.path.join(cbc.Configuration.sampling_output_dir, 'step%06d.coord' % (self.counter)))


    def save_top(self, n = 100000, counter=100, step_save=0):
        '''
        Save the top n structures.
        '''
        # if we don't want to save any structures, then don't save any structures
        if n == 0:
            return

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

        if self.step_save > 0:
            if self.counter % self.step_save == 0:
                sorted_energies[0][2].to_cg_file(os.path.join(cbc.Configuration.sampling_output_dir, 'intermediate_best%d.coord' % (counter)))

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

        print ("---------------------------")

        print ([e[1] for e in sorted_energies[:10]])

        '''
        for energy in energy_function.energies:
            print energy.__class__.__name__, energy.eval_energy(sm)
        '''

        print ("-------------------------")

def load_sampled_elements(sm):
    """
    Try to load the sampled elements from the cg into the sm.
    :returns: True upon success, False upon failure (e.g. if cg was derived from a fasta file)
    """
    try:
        sm.load_sampled_elems()
    except Exception as e:
        log.warning("Cannot use structure from file. Need to resample: {}".format(e), exc_info=True)
        return False
    if not sm.elem_defs:
        return False
    return True

def build_fair(sm, stat_source, target_attempts=None, target_structures=None, randomize_mst = False):
    if target_attempts is None and target_structures is None:
        raise ValueError("Need target_structures or target_attempts")
    log.debug("Target attempts: {}, target structures: {}".format(target_attempts, target_structures))
    attempts = 0
    junction_failures = 0
    structures = []        
    constraint_energy = sm.constraint_energy
    junction_constraint_energy = sm.junction_constraint_energy
    failed_structures = []
    sm.constraint_energy = None
    sm.junction_constraint_energy = None
    while True:
        if sys.stdout.isatty():
            print("\r{} / {}".format(len(structures), attempts), end="")
        attempts += 1
        if target_attempts is not None and attempts>target_attempts:
            break
        log.debug("Fair building: Sampling...")
        sm.sample_stats(stat_source)
        sm.traverse_and_build()
        if constraint_energy.eval_energy(sm)==0 and junction_constraint_energy.eval_energy(sm)==0:
            structures.append(sm.bg.to_cg_string()) #Serializing is cheaper than deepcopy.
            log.debug("... is valid")
        else:
            failed_structures.append(sm.bg.to_cg_string())
        if junction_constraint_energy.eval_energy(sm)>0:
            junction_failures += 1
            log.debug("... junction invalid")
        else:
            log.debug("... clash invalid")

        if target_structures is not None and len(structures)>=target_structures:
            break
        if randomize_mst:
            for ml in sm.bg.find_mlonly_multiloops:
                segment = random.choice(ml)
                sm.set_multiloop_break_segment(segment)
    sm.constraint_energy = constraint_energy
    sm.junction_constraint_energy = junction_constraint_energy
    return structures, failed_structures, attempts, junction_failures

def build_sm(sm, stat_source, verbose = False):
    """
    Build the initial structure that will be used for sampling.
    """
    if sm.constraint_energy is not None or sm.junction_constraint_energy is not None:    
        constraint_energy = sm.constraint_energy
        junction_constraint_energy = sm.junction_constraint_energy
        sm.constraint_energy = None
        sm.junction_constraint_energy = None
        log.info("building without constraint energy...")
        sm.traverse_and_build()
        sm.constraint_energy = constraint_energy
        sm.junction_constraint_energy = junction_constraint_energy
        log.info("building with constraint energy (verbose is '{}')".format(verbose))
        sm.traverse_and_build(verbose=verbose, stat_source=stat_source)
        log.info("finished building")
    sm.traverse_and_build()
    sm.bg.add_all_virtual_residues()


class MCMCSampler(object):
    '''
    Sample using tradition accept/reject sampling.
    '''
    def __init__(self, sm, energy_function, stats, stat_source, dump_measures=False):
        '''
        :param sm: SpatialModel that will be used for sampling.
    
        '''
        #BT: Seems to be not in use
        #if stats_type == 'continuous':
        #    self.cont_stats = ftms.ContinuousAngleStats(ftms.get_angle_stats())
        #elif stats_type == 'random':
        #    self.cont_stats = ftms.RandomAngleStats(ftms.get_angle_stats())
        #else:
        #    self.cont_stats = None
        
        self.step_counter = 0
        self.sm = sm
        self.energy_function = energy_function
        self.stats = stats
        self.stat_source = stat_source
        if not isinstance(stats, sstats.SamplingStatistics):
            self.stats.energy_function = energy_function
        #: Store the energy of the last configuration!
        self.prev_energy = 100000000000.
        self.dump_measures = dump_measures
        self.resampled_energy = True
        self.prev_stats={}
        self.prev_mst = None #Used only in Subclass!
        if sm.constraint_energy is not None:
            self.energy_function.energies += sm.constraint_energy.energies
        if sm.junction_constraint_energy is not None:
            self.energy_function.energies += [sm.junction_constraint_energy]
        self.sm.constraint_energy = None
        self.sm.junction_constraint_energy = None

        #Evaluate the energy so we can accept that measure.
        self.prev_energy = energy_function.eval_energy(sm)
        log.info("Energy was {}".format(self.prev_energy))
        try:
            self.prev_constituing = self.energy_function.constituing_energies
        except AttributeError: 
            pass
        #Accept the measure of the initial structure.
        #This is required so reject_last_measure does not accept the last measure from the 
        #file a second time. And for the use_accepted_measure flag of eval_energy
        self.energy_function.accept_last_measure()
        self.sm.get_sampled_bulges() #Store in sm which bulges are sampled (vs broken ml-segments)
        if isinstance(stats, sstats.SamplingStatistics):
            self.stats.print_header()

    #@profile
    def change_elem(self):
        '''
        Change a random element and accept the change with a probability
        proportional to the energy function.
        '''

        # pick a random element and get a new statistic for it
        possible_elements=list(self.sm.bg.get_mst())
        pe=set(possible_elements)
        d = random.choice(possible_elements)
        movestring =  self.change_one_element(d)
        ms, accepted = self.accept_reject()
        return movestring + ms, accepted

    def accept_reject(self):
        movestring=[]
        energy = self.energy_function.eval_energy(self.sm, background=True)
        movestring.append("{:.3f}".format(self.prev_energy))
        movestring.append("->")
        movestring.append("{:.3f};".format(energy))
        if energy <= self.prev_energy:
            movestring.append("A")
            # lower energy means automatic acceptance accordint to the
            # metropolis hastings criterion
            self.accept(energy)
            accepted = True
        else:
            # calculate a probability
            r = random.random()
            if r > math.exp(self.prev_energy - energy):
                movestring.append("R")
                # reject the sampled statistic and replace it the old one
                cg_stri = self.sm.bg.to_cg_string()
                with open(os.path.join(conf.Configuration.sampling_output_dir, 'before_reset.coord'), "w") as f:
                    f.write(cg_stri)

                self.reject()
                self.sm.new_traverse_and_build(start='start')
                accepted = False
                if not self.energy_function.uses_background():
                    rec_prev_energy = self.energy_function.eval_energy(self.sm)
                    log.debug("Energy after resetting is {}".format(rec_prev_energy))
                    if rec_prev_energy != self.prev_energy:
                        cg_stri = self.sm.bg.to_cg_string()
                        with open(os.path.join(conf.Configuration.sampling_output_dir, 'after_reset.coord'), "w") as f:
                            f.write(cg_stri)
                    assert rec_prev_energy == self.prev_energy, "{}!={}. Energy changed after resetting".format(rec_prev_energy, self.prev_energy)
            else:            
                movestring.append("A")
                self.accept(energy)
                accepted = True
                #print ("...still accepting")            
        return "".join(movestring), accepted

    def accept(self, energy):
        """
        :param energy: The energy of the accepted state. 
                       This is to avoid expensive recalculation of the energy
        """
        # accept the new statistic
        self.prev_energy = energy
        try:
            self.prev_constituing =  self.energy_function.constituing_energies
        except AttributeError: pass
        self.energy_function.accept_last_measure()
        for d in self.prev_stats:
            if d[0] =="m" and d not in self.sm.bg.mst:
                if d in self.sm.bg.sampled: del self.sm.bg.sampled[d]
        self.prev_mst = None
        self.prev_stats = {}
        if isinstance(self.energy_function, fbe.CombinedEnergy):
            for e in self.energy_function.iterate_energies():
                if hasattr(e, "accepted_projDir"):
                    self.sm.bg.project_from=e.accepted_projDir

    def reject(self):
        if self.prev_mst is not None:
            oldonly = self.prev_mst - self.sm.bg.mst
            self.sm.change_mst(self.prev_mst)
            for m in oldonly:
                if m[0]!="m": continue
                if m in self.sm.bg.sampled:
                    del self.sm.bg.sampled[m]        
            self.prev_mst = None
        for d, stats in self.prev_stats.items():
            log.debug("reject: Resetting {} to {}".format(d, stats.pdb_name))
            self.sm.elem_defs[d] = stats
        self.energy_function.reject_last_measure()
        self.prev_stats = {}
        
    def change_one_element(self, d):
        """
        Change the stats for the selected define and accept or 
        reject the new structure based on the energy.

        :param d: The define to change.
        """
        
        movestring=[]
        movestring.append(d+":")
        possible_stats=self.stat_source.get_possible_stats(self.sm.bg, d)
        new_stat = random.choice(possible_stats)

        # we have to replace the energy because we've probably re-calibrated
        # the energy function
        if self.resampled_energy and self.energy_function.uses_background():
            self.prev_energy = self.energy_function.eval_energy(self.sm, background=True, use_accepted_measure=True)
            try:
                self.prev_constituing =  self.energy_function.constituing_energies
            except AttributeError: pass
            self.resampled_energy = False

        self.prev_stats={d: self.sm.elem_defs[d]}
        movestring.append(self.prev_stats[d].pdb_name)
        movestring.append("->")
        movestring.append(new_stat.pdb_name)
        movestring.append(";")
        #if isinstance(new_stat, ftms.AngleStat) and isinstance(self.sm.conf_stats.angle_stats, ftms.ClusteredAngleStats):
        #    movestring.append(str(self.sm.conf_stats.angle_stats.cluster_of(self.prev_stats[d])))
        #    movestring.append("->")
        #    movestring.append(str(self.sm.conf_stats.angle_stats.cluster_of(new_stat)))
        #    movestring.append(";")
        self.sm.elem_defs[d] = new_stat
        self.sm.new_traverse_and_build(start=d, include_start=True)
        return "".join(movestring)

    def step(self):
        movestring, accepted =self.change_elem()
        
        if self.dump_measures:
            if self.step_counter % 20 == 0:
                self.energy_function.dump_measures(cbc.Configuration.sampling_output_dir, 
                                                   self.step_counter)

        if self.step_counter % 3 == 0:
            self.resampled_energy = True
            self.energy_function.resample_background_kde(self.sm.bg)

        self.step_counter += 1
        if isinstance(self.stats, sstats.SamplingStatistics):
            self.stats.update_statistics( self.sm, self.prev_energy, self.prev_constituing, movestring )

        self.energy_function.update_adjustment(self.step_counter, self.sm.bg)
        
        return accepted
class ImprovedMultiloopMCMC(MCMCSampler):
    """
    This Sampler is like the MCMCSampler except for multiloops.

    If a multiloop segment is picked, a different move step is used:
    After changing the stats for the selected multiloop segment, 
    the loop closure energy is evaluated.
    If it is >0, OTHER multiloop segments are changed, 
    until a certain number of tries was performed or the multiloop fulfills the constraints.
    
    Stats that lead to a lot of clashes are picked less often for the first multiloop segment, to
    avoid introducing a statistical bias.
    """
    def __init__(self, sm, energy_function, stats, stat_source, dump_measures=False):        
        self.junction_energy = sm.junction_constraint_energy
        super(ImprovedMultiloopMCMC, self).__init__(sm, energy_function, stats, stat_source, dump_measures)
        #: A dict of the form `{ define : { stat : [#choosen, #rejects ] }}`
        self.stats_weights = c.defaultdict(lambda : c.defaultdict(lambda: [0,0]))

        self.prev_stats={}
        self.prev_mst=None
    def change_elem(self):
        # pick a random element and get a new statistic for it
        possible_elements=list(self.sm.bg.mst)
        pe=set(possible_elements)
        d = random.choice(possible_elements)
        self.prev_stats={}
        self.prev_mst=None
        if d[0]!="m":
            movestring=self.change_one_element(d) #Use the superclass method.
            ms, accepted = self.accept_reject()
            return movestring + ms, accepted

        junction_nodes = set( x for x in self.sm.bg.find_bulge_loop(d, 200) if x[0]=="m" )
        missing_nodes = junction_nodes - pe
        defined_junction_nodes = junction_nodes &  pe

        for node in defined_junction_nodes: #self.prev_stats was emptied earlier in this function
            self.prev_stats[node]= self.sm.elem_defs[node]
        if not defined_junction_nodes: #Open multiloop (at terminus), not a cycle
            self.prev_stats[d]= self.sm.elem_defs[d]
        assert d in self.prev_stats

        # Break another multiloop segment and sample stats for this segment!
        if len(missing_nodes)>0:
            self.prev_mst = copy.copy(self.sm.bg.mst)
            #print("Breaking {}".format(d))
            d = self.sm.set_multiloop_break_segment(d) #This is done AFTER prev_stats are saved! d now is the previousely broken segment
            #print("Closed {}".format(d))

        #Update defined_junction_nodes
        defined_junction_nodes = junction_nodes & self.sm.bg.mst

        movestring=[]
        movestring.append(d+":") #The previousely broken but now closed segment.
        possible_stats=self.stat_source.get_possible_stats(self.sm.bg, d)
        random.shuffle(list(possible_stats))

        searching=True
        #For this multiloop segment, pick stats that have fewer rejects more often!
        xth_try=0
        while searching:
            xth_try+=1
            for newstat in possible_stats:
                r = random.random()
                if self.stats_weights[d][newstat][0] == 0:
                    weight=1
                else:
                    weight=1-self.stats_weights[d][newstat][1]/(self.stats_weights[d][newstat][0])
                if r < weight:
                    searching=False
                    break
            if xth_try>10000:
                break
        self.stats_weights[d][newstat][0]+=1
        movestring.append("ST{};".format(xth_try))

        self.sm.elem_defs[d] = newstat
        self.sm.traverse_and_build(start=d)
        if self.junction_energy is not None:
            energy = self.junction_energy.eval_energy(self.sm, nodes=junction_nodes)
        else:
            energy = 0
        num_tries=0
        if len(defined_junction_nodes)>1: #If pseudoknots are allowed, more than one multiloop segment can be broken.
            #print ("... coaxial stacks now: {}".format(rcs.report_all_stacks(self.sm.bg)), file=sys.stderr)
            while energy>0:
                #The previous attempt is failed!
                self.stats_weights[d][newstat][1]+=1
                num_tries+=1
                if num_tries>10000:
                    #There might not be any conformation possible for this stat.
                    #Just return a wrong conformation and let the energy function take care of it.
                    break 
                #We give it another try. Next attempt
                self.stats_weights[d][newstat][0]+=1
                other_d = random.choice(list(set(defined_junction_nodes)-set([d])))
                #print ("... now changing {}".format(other_d), file=sys.stderr)
                assert other_d in self.prev_stats
                possible_other_stats = self.stat_source.get_possible_stats(self.sm.bg, other_d)
                self.sm.elem_defs[other_d] = random.choice(possible_other_stats)
                self.sm.traverse_and_build(start=other_d)
                #print ("... coaxial stacks now: {}".format(rcs.report_all_stacks(self.sm.bg)), file=sys.stderr)
                assert self.junction_energy is not None
                energy = self.junction_energy.eval_energy(self.sm, nodes=junction_nodes)
        movestring.append("TRIES{};".format(num_tries))
        ms, accepted = self.accept_reject()
        movestring.append(ms)
        return "".join(movestring), accepted

class ExhaustiveExplorer(MCMCSampler):
    def __init__(self, sm, energy_function, stats, stat_source, loop_of_interest):
        super(ExhaustiveExplorer, self).__init__(sm, energy_function, stats, stat_source)
        self.possible_stats = self.stat_source.get_possible_stats(sm.bg, loop_of_interest)
        self.loop_of_interest=loop_of_interest
    def next_choice(self):
        while True:
            random.shuffle(self.possible_stats)
            for stat in self.possible_stats:
                yield stat
    def change_elem(self):
        new_stat=next(self.next_choice())
        self.sm.elem_defs[self.loop_of_interest] = new_stat
        self.sm.traverse_and_build(start=self.loop_of_interest)
        self.sm.bg.infos["sampledFromExhaustive"]=["{} {}".format(self.loop_of_interest,c)]
        #Calculate the enrgy, but always accept
        energy = self.energy_function.eval_energy(self.sm, background=True)
        self.prev_energy = energy
        try:
            self.prev_constituing =  self.energy_function.constituing_energies
        except AttributeError: pass
        self.energy_function.accept_last_measure()
        return str(c), True

"""
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
        (dist, size1, size2, type1) = ftms.get_angle_stat_dims(dims[0], dims[1], ang_type1)[0]
        possible_angles = ftms.get_angle_stats()[(size1, size2, ang_type1)]

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


        # energy = -log(p(S)) 
        # So... we want to maximize p(S)
        # Therefore, we want to minimize the energy
        max_energy = max(energies.values())
        min_energy = min(energies.values())

        if max_energy - min_energy > 40:
            max_energy = min_energy + 40.

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

        # sanity check
        #assert(allclose(total_prob, 1.))

        #pick one new angle to accept given the probabilities of the
        #sampled ones
        prob_remaining = 1.
        for key in energy_probs.keys():
            if random.random() < energy_probs[key] / prob_remaining:
                self.sm.angle_defs[bulge][ang_type1] = key
                #self.sm.traverse_and_build(start=bulge)
                break

            prob_remaining -= energy_probs[key]

        self.sm.traverse_and_build(start=bulge)
        self.stats.update_statistics(self.energy_function, self.sm)
"""
