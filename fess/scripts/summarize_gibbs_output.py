#!/usr/bin/python

import sys, re
import os.path as op
import collections as c
import numpy as np
import scipy.stats as ss
import math as m

import corgy.utilities.debug as cud

from optparse import OptionParser

def density_visit_dir(rmsds, dirname, names):
    if 'out.txt' not in names:
        return

    filename = op.join(dirname, 'out.txt')
    trial_id = op.split(dirname)[1]

    with open(filename) as f:
        lines = f.readlines()
        for line in lines:
            if line.find("native_energy") == 0:
                parts =  line.split()

                pdb_name = parts[1].strip('][')

                # keep track of the energy and rmsd of 
                # the sampled structures
                rmsds[pdb_name] += [(float(parts[4]), float(parts[5]))]

def visit_dir(rmsds, dirname, names):
    if 'out.txt' not in names:
        return

    filename = op.join(dirname, 'out.txt')
    trial_id = op.split(dirname)[1]

    with open(filename) as f:
        file_size = op.getsize(filename) 
        if file_size < 2048:
            return
        f.seek(-1024, 2)
        line = f.readlines()[-1]

        if line.find("native_energy") == 0:
            m = re.search('\[(.*) (.*)\].*min:[ ]+(.*)[ ]+\((.*)\)+[ ](.*)', line) 

            # group(1) = pdb id
            # group(2) = length in nucleotides
            # group(3) = energy
            # group(4) = rmsd
            rmsds[(int(m.group(2)), m.group(1))] += [(float(m.group(3)), float(m.group(5)), float(m.group(4)), filename, trial_id)]

def summarize_rmsds(rmsds, compact=False, base_dir='', nth=0):
    keys = rmsds.keys()
    #keys.sort(key=lambda x: min([k[1] for k in rmsds[x]]))
    keys.sort()

    min_rmsds = []
    for key1,key2 in keys:
        key = (key1,key2)
        rmsds[key].sort(key=lambda x: x[0])
        if compact:
            print base_dir, key2, rmsds[key][nth][4]
        else:
            print "[",key2 + "/" + rmsds[key][nth][4],"|", key1,"]", "[", rmsds[key][nth][1], "::", str(rmsds[key][nth][0]) , "(" , str(rmsds[key][nth][2]), ")" , "]", " ".join([str(k[1]) for k in rmsds[key][1:5]])

        min_rmsds += [rmsds[key][nth][1]]
    min_rmsds.sort()

    if not compact:
        print "average: %.2f median %.2f" % (np.mean(np.array(min_rmsds)), min_rmsds[len(min_rmsds) / 2])

def summarize_rmsds_by_density(rmsds, plot=False, random_rmsds = None):
    if plot:
        import matplotlib.pyplot as plt
        plt.ioff()

    counter = 0
    num_figs = len(rmsds.keys())
    rows = 3
    cols = int(m.ceil(num_figs / float(rows)))
    densest_rmsds = c.defaultdict(list)

    if plot:
        fig = plt.figure(figsize=(20,16))

        # Adjust the numbering so that the histograms are under the plots
        plot_nums =[i for i in range(1, rows * cols * 2 + 1)]
        print "plot_nums:", plot_nums, len(plot_nums)
        plot_nums = np.array(plot_nums)
        print rows, cols
        plot_nums = plot_nums.reshape((rows * 2, cols))
        print plot_nums
        print plot_nums.T
        plot_nums =  np.array(plot_nums.T).reshape((rows * cols * 2,))
        plot_nums

    for pdb_name in rmsds.keys():
        energy_rmsds = np.array(rmsds[pdb_name])
        max_energy = 5.
        low_energy_rmsds = energy_rmsds[energy_rmsds[:,0] < max_energy]
        xs = None

        if random_rmsds != None and pdb_name in random_rmsds.keys():
            rr = np.array(random_rmsds[pdb_name])
            low_random_rmsds = rr[rr[:,0] < 30]

            if len(low_random_rmsds[:,1]) < 5:
                # this can occur if we've sampled a lot of high energy
                # structures
                #
                # Just ignore it and move onto the next structure
                continue
            xs = np.linspace(0., max(low_random_rmsds[:,1]), num=200)


        try:
            if xs == None:
                xs = np.linspace(0., max(low_energy_rmsds[:,1]), num=200)
            kde = ss.gaussian_kde(low_energy_rmsds[:,1])
        except ValueError as ve:
            # Usually happens when there aren't enough low energies
            # to create a kde
            continue

        kxs = kde(xs)

        if random_rmsds != None:
            random_kde = ss.gaussian_kde(low_random_rmsds[:,1])
            rand_kxs = random_kde(xs)

        if plot:


            row = counter / 3
            col = counter - 3 * (counter / 3)

            print "row:", row, "col:", col

            ax_plot = fig.add_subplot(rows*2, cols, plot_nums[counter])
            ax_hist = fig.add_subplot(rows*2, cols, plot_nums[counter+1])

            ax_plot.set_title(pdb_name)
            ax_plot.set_xlim(min(xs), max(xs))

            ax_hist.set_title(pdb_name)

            ax_plot.plot(low_energy_rmsds[:,1], low_energy_rmsds[:,0], 'bo', alpha=0.2)
            ax_hist.set_xlim(min(xs), max(xs))
            ax_hist.plot(xs, kde(xs), 'b')

            if random_rmsds != None:
                ax_hist.plot(xs, random_kde(xs), 'r')
                #ax_plot.plot(low_random_rmsds[:,1], low_random_rmsds[:,0], 'ro')

        counter += 2

        if random_rmsds != None:
            print pdb_name, xs[kxs == max(kxs)], xs[rand_kxs == max(rand_kxs)]
            densest_rmsds[pdb_name] = [xs[kxs == max(kxs)][0], xs[rand_kxs == max(rand_kxs)][0]]
        else:
            print pdb_name, xs[kxs == max(kxs)][0]
            densest_rmsds[pdb_name] = [xs[kxs == max(kxs)][0]]


    if plot:
        fig.tight_layout()
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


def main():
    usage = './summarize_gibbs_output.py base_dir'
    parser = OptionParser()

    parser.add_option('-n', '--nth-best', dest='nth_best', default=0, help="Display the n-th best energy structure (0-based)", type='int')
    parser.add_option('-c', '--compact', dest='compact', default=False, action='store_true', help='Display a compact representation of the structures consisting of only the name and number.')
    parser.add_option('-d', '--density', dest='density', default=False, action='store_true', help='Display the rmsd of the structure with the greatest probability density in terms of the sampled structures.')
    parser.add_option('-p', '--plot', dest='plot', default=False, action='store_true', help='Plot the densities.')
    parser.add_option('-r', '--random', dest='random', default=None, help='Location of the random energies')

    (options, args) = parser.parse_args()

    if len(args) < 1:
        parser.print_help()
        sys.exit(1)

    rmsds = c.defaultdict(list)
    if options.density:
        random_rmsds = None
        if options.random != None:
            random_rmsds = c.defaultdict(list)
            op.walk(options.random, density_visit_dir, random_rmsds)

        op.walk(args[0], density_visit_dir, rmsds)
        summarize_rmsds_by_density(rmsds, options.plot, random_rmsds=random_rmsds)
    else:
        for arg in args:
            op.walk(arg, visit_dir, rmsds)

        summarize_rmsds(rmsds, options.compact, args[0], nth=options.nth_best)

if __name__ == '__main__':
    main()

