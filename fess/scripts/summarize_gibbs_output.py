#!/usr/bin/python

import math
import sys, re
import os.path as op
import collections as c
import numpy as np
import scipy.stats as ss
import math as m

import borgy.utilities.debug as cud

from optparse import OptionParser

def density_visit_dir(rmsds, dirname, names):
    if 'log.txt' not in names:
        return

    filename = op.join(dirname, 'log.txt')
    trial_id = op.split(dirname)[1]

    with open(filename) as f:
        lines = f.readlines()
        for line in lines:
            if line.find("native_energy") == 0:
                parts =  line.split()

                pdb_name = parts[1].strip('][')
                pdb_size = int(parts[2].strip(']:'))

                # keep track of the energy and rmsd of 
                # the sampled structures
                rmsds[(pdb_name, pdb_size)] += [(float(parts[4]), float(parts[5]), int(parts[3]))]

def visit_dir(rmsds, dirname, names):
    if 'log.txt' not in names:
        return

    filename = op.join(dirname, 'log.txt')
    trial_id = op.split(dirname)[1]

    with open(filename) as f:
        file_size = op.getsize(filename) 
        if file_size < 2048:
            return
        f.seek(-1024, 2)
        line = f.readlines()[-1]

        if line.find("native_energy") == 0:
            m = re.search('\[(.*) (.*)\].*([0-9]+).*min:[ ]+(.*)[ ]+\((.*)\)+[ ](.*)[ ]\|', line) 
            m1 = re.search('\[(.*) (.*)\].*extreme_rmsds:[ ]+(.*)[ ]+(.*)', line) 

            # group(1) = pdb id
            # group(2) = length in nucleotides
            # group(3) = iteration
            # group(4) = energy
            # group(5) = rmsd
            rmsds[(int(m.group(2)), m.group(1))] += [(float(m.group(4)), float(m.group(6)), float(m.group(5)), filename, trial_id, float(m1.group(3)), float(m1.group(4)), int(m.group(3)))]

def summarize_rmsds(rmsds, compact=False, base_dir='', nth=0, minimal=False):
    keys = rmsds.keys()
    #keys.sort(key=lambda x: min([k[1] for k in rmsds[x]]))
    keys.sort()

    min_rmsds = []
    ratios = []
    best_rmsds = []
    worst_rmsds = []


    for key1,key2 in keys:
        # key1 is the length of the sequence and key2 is the pdb name
        key = (key1,key2)
        rmsds[key].sort(key=lambda x: x[0])

        if minimal:
            print "ernwin", key2, rmsds[key][nth][1]
        elif compact:
            print base_dir, key2, rmsds[key][nth][4]
        else:
            pdb_name = key2 
            pdb_len = key1
            trial_num = rmsds[key][nth][4]
            rmsd = rmsds[key][nth][1]
            energy = rmsds[key][nth][0]
            native_energy = rmsds[key][nth][2]

            best_rmsd = rmsds[key][nth][5]
            worst_rmsd = rmsds[key][nth][6]

            ratio = (rmsd - best_rmsd) / (worst_rmsd - best_rmsd)

            print "[",pdb_name + "/" + trial_num,"|", pdb_len,"]", "[", str(rmsd), "::", str(energy) , "(" , str(native_energy), ")" , "]", "<", "%.2f" % (ratio), ">", str(best_rmsd), str(worst_rmsd)

            min_rmsds += [rmsds[key][nth][1]]
            ratios += [ratio]
            best_rmsds += [best_rmsd]
            worst_rmsds += [worst_rmsd]

    min_rmsds.sort()

    if not compact and not minimal:
        print "average: %.2f median %.2f avg_ratio: %.2f avg_best_rmsd: %.2f avg_worst_rmsd: %.2f" % (np.mean(np.array(min_rmsds)), min_rmsds[len(min_rmsds) / 2], np.mean(np.array(ratios)), np.mean(np.array(best_rmsds)), np.mean(np.array(worst_rmsds)))

def summarize_rmsds_by_density(all_rmsds, plot=False):
    if plot:
        import matplotlib.pyplot as plt
        plt.ioff()

    counter = 0
    num_figs = len(all_rmsds[0].keys())
    rows = int(math.floor(math.sqrt(num_figs)))
    cols = int(m.ceil(num_figs / float(rows)))
    densest_rmsds = c.defaultdict(list)

    if plot:
        fig = plt.figure(figsize=(rows*5,rows*4))

        # Adjust the numbering so that the histograms are under the plots
        plot_nums =[i for i in range(1, rows * cols + 1)]
        #print "plot_nums:", plot_nums, len(plot_nums)
        plot_nums = np.array(plot_nums)
        #print rows, cols
        plot_nums = plot_nums.reshape((rows, cols))
        #print plot_nums
        #print plot_nums.T
        plot_nums =  np.array(plot_nums.T).reshape((rows * cols,))

    colors = ['#d7191c','#fdae61','#ffffbf','#abdda4', '#2b83ba']

    for (current_color, rmsds) in zip(colors, all_rmsds):
        counter = 0
        keys = rmsds.keys()
        keys.sort(key=lambda x: x[1])
        for (pdb_name, pdb_size) in keys:
            energy_rmsds = np.array(rmsds[(pdb_name, pdb_size)])
            max_energy = 10000000
            low_energy_rmsds = energy_rmsds[energy_rmsds[:,0] < max_energy]
            xs = None

            try:
                if xs == None:
                    xs = np.linspace(0., max(low_energy_rmsds[:,1]), num=200)
                kde = ss.gaussian_kde(low_energy_rmsds[:,1])
            except ValueError as ve:
                # Usually happens when there aren't enough low energies
                # to create a kde
                continue

            kxs = kde(xs)

            if plot:
                row = counter / rows
                col = counter - rows * (counter / rows)

                #print "row:", row, "col:", col

                #ax_plot = fig.add_subplot(rows, cols, plot_nums[counter])
                ax_hist = fig.add_subplot(rows, cols, plot_nums[counter])

                #ax_plot.set_title(pdb_name)
                #ax_plot.set_xlim(min(xs), max(xs))

                ax_hist.set_title(pdb_name)

                #ax_plot.plot(low_energy_rmsds[:,1], low_energy_rmsds[:,0], 'o', alpha=0.2)
                ax_hist.set_xlim(min(xs), max(xs))
                ys = kde(xs)
                ax_hist.plot(xs, list(kde(xs)[:-1]) + [0.], alpha=0.5, color=current_color, linewidth=3.)

            counter += 1
            densest_rmsds[pdb_name] = [xs[kxs == max(kxs)][0]]
            print "%s %d %.1f %.1f" % (pdb_name, pdb_size, xs[kxs == max(kxs)][0], np.average(xs, weights=kxs))

    if plot:
        fig.tight_layout()
        plt.savefig('current_plot.png', bbox_inches='tight')
        plt.show()

    sum_rmsds = [0. for i in densest_rmsds.items()[0][1]]
    total = 0.
    for (key, value) in densest_rmsds.items():
        for i,val in enumerate(value):
            sum_rmsds[i] += val
        total += 1
    for s in sum_rmsds:
        print "---------------"
        print "avg:", s / total

def calc_stepwise_divergence(all_rmsds):
    for rmsds in all_rmsds:
        for pdb_name in sorted(rmsds.keys()):
            energy_rmsds = np.array(rmsds[pdb_name])

            kde = ss.gaussian_kde(energy_rmsds[:,1])
            xs = np.linspace(0., max(energy_rmsds[:,1]), num=200)
            ys = kde(xs)

            max_iterations = int(max(energy_rmsds[:,2]))

            for i in range(100, max_iterations, 100):
                curr_energy_rmsds = energy_rmsds[energy_rmsds[:,2] < i]
                new_kde = ss.gaussian_kde(curr_energy_rmsds[:,2])

                ys_new = new_kde(xs)
                q = ys
                p = ys_new
                div = p / q
                ldiv = np.log(div)
                kd = np.dot(ys, ldiv)

                print i, kd

def main():
    usage = './summarize_gibbs_output.py base_dir'
    parser = OptionParser()

    parser.add_option('-n', '--nth-best', dest='nth_best', default=0, help="Display the n-th best energy structure (0-based)", type='int')
    parser.add_option('-c', '--compact', dest='compact', default=False, action='store_true', help='Display a compact representation of the structures consisting of only the name and number.')
    parser.add_option('-d', '--density', dest='density', default=False, action='store_true', help='Display the rmsd of the structure with the greatest probability density in terms of the sampled structures.')
    parser.add_option('-p', '--plot', dest='plot', default=False, action='store_true', help='Plot the densities.')
    parser.add_option('-r', '--random', dest='random', default=None, help='Location of the random energies')
    parser.add_option('-m', '--minimal', dest='minimal', default=False, action='store_true', help='Display a minimal output consisting of the pdb name and the rmsd')
    parser.add_option('-s', '--steps', dest='steps', default=False, action='store_true', help='Display the KL divergence between the true distribution and the sampled one at each step of the simulation')

    (options, args) = parser.parse_args()

    if len(args) < 1:
        parser.print_help()
        sys.exit(1)

    # summarize by the lowest energy
    rmsds = c.defaultdict(list)
    for arg in args:
        op.walk(arg, visit_dir, rmsds)
    summarize_rmsds(rmsds, options.compact, args[0], nth=options.nth_best, minimal=options.minimal)

    # summarize by the structure with the greatest density
    rmsds = c.defaultdict(list)
    if options.density:
        all_rmsds = []

        for arg in args:
            these_rmsds = c.defaultdict(list)
            op.walk(arg, density_visit_dir, these_rmsds)

            all_rmsds += [these_rmsds]

        if options.steps:
            calc_stepwise_divergence(all_rmsds)
        else:
            summarize_rmsds_by_density(all_rmsds, options.plot)



if __name__ == '__main__':
    main()

