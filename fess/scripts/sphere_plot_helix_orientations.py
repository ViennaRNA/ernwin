#!/usr/bin/python

import random
import numpy as np
import math as m
from matplotlib import cm
from matplotlib.ticker import LinearLocator
import sys, pandas as pa
import scipy.stats as ss
from optparse import OptionParser
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d, Axes3D

from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import proj3d

class Arrow3D(FancyArrowPatch):
    '''
    From:
        
    http://stackoverflow.com/questions/11140163/python-matplotlib-plotting-a-3d-cube-a-sphere-and-a-vector
    '''

    def __init__(self, xs, ys, zs, *args, **kwargs):
        FancyArrowPatch.__init__(self, (0,0), (0,0), *args, **kwargs)
        self._verts3d = xs, ys, zs

    def draw(self, renderer):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, renderer.M)
        self.set_positions((xs[0],ys[0]),(xs[1],ys[1]))
        FancyArrowPatch.draw(self, renderer)


def main():
    usage = """
    ./sphere_plot_helix_orientations.py 
    """
    num_args=0
    parser = OptionParser()

    parser.add_option('-r', '--resolution', dest='resolution', default=10, help="The resolution of the resulting plot", type='int')
    parser.add_option('-a', '--angle', dest='angle', default=0, help="The angle of the camera", type='float')
    parser.add_option('-f', '--fig-name', dest='fig_name', default='', help="The name of the file to save the figure to. If it is not specified, the figure will not be saved", type='str')

    #parser.add_option('-u', '--useless', dest='uselesss', default=False, action='store_true', help='Another useless option')

    (options, args) = parser.parse_args()

    if len(args) < num_args:
        parser.print_help()
        sys.exit(1)

    column_names = ['type', 'pdb', 's1', 's2', 'u', 'v', 't', 'r', 'u1', 'v1', 'atype', 'something1', 'something2', 'sth3', 'sth4']
    real_stats = pa.read_csv('fess/stats/temp.real.stats', header=None, sep=' ', names=column_names, engine='python')
    sampled_stats = pa.read_csv('fess/stats/temp.sampled.stats', header=None, sep=' ', names=column_names, engine='python')

    real_stats = real_stats[real_stats["type"] == "angle"]
    real_us = real_stats[['u', 'v']].as_matrix()
    sampled_us = sampled_stats[['u','v']].as_matrix()

    real_us_orig = np.copy(real_us)
    sampled_us_orig = np.copy(sampled_us)

    real_us = np.vstack([real_us+[0,-2*m.pi], 
                         real_us+[0,0], 
                         real_us+[0,2*m.pi], 
                         real_us+[m.pi,-2*m.pi], 
                         real_us+[m.pi,0], 
                         real_us+[m.pi, 2*m.pi],
                         real_us+[-m.pi,-2*m.pi],
                         real_us+[-m.pi,0], 
                         real_us+[-m.pi, 2*m.pi]] )
    sampled_us = np.vstack([sampled_us+[0,-2*m.pi], 
                         sampled_us+[0,0], 
                         sampled_us+[0,2*m.pi], 
                         sampled_us+[m.pi,-2*m.pi], 
                         sampled_us+[m.pi,0], 
                         sampled_us+[m.pi, 2*m.pi],
                         sampled_us+[-m.pi,-2*m.pi],
                         sampled_us+[-m.pi,0], 
                         sampled_us+[-m.pi, 2*m.pi] ])

    #real_us = real_us[random.sample(range(len(real_us)), len(real_us)/3)]
    #sampled_us = sampled_us[random.sample(range(len(sampled_us)), len(sampled_us)/3)]
    
    print "len(real_us):", len(real_us)
    k_ro = ss.gaussian_kde(real_us_orig.T, bw_method="silverman")
    k_so = ss.gaussian_kde(sampled_us_orig.T, bw_method="silverman")

    print len(real_us), len(sampled_us)

    k_r = ss.gaussian_kde(real_us.T, bw_method=k_ro.factor/4.)
    k_s = ss.gaussian_kde(sampled_us.T, bw_method=k_so.factor/4.)

    '''
    k_r = k_ro
    k_s = k_so
    '''

    #U,V = np.mgrid[0:m.pi:complex(0,options.resolution), -m.pi:m.pi:complex(0, options.resolution)]
    U,V = np.mgrid[0:m.pi:complex(0,options.resolution), -m.pi:m.pi:complex(0, options.resolution)]
    positions = np.vstack([U.ravel(), V.ravel()])
    vals_r = np.log(k_r(positions))
    vals_s = np.log(k_s(positions))
    vals = 0.7 * (vals_r - vals_s)
    #vals = 3. * (np.exp(vals_r) - np.exp(vals_s))
    #print "vals_r:", vals_r
    #print "vals_s:", vals_s
    #print "vals:", vals, 
    print len(vals), len(vals[vals > 0]), np.mean(vals)

    fig = plt.figure(figsize=(10,10))
    ax = Axes3D(fig)

    ulen = len(U)
    vlen = len(V)

    X = np.sin(U) * np.cos(V)
    Y = np.sin(U) * np.sin(V)
    Z = np.cos(U)


    W = np.reshape(vals.T, U.shape)
    colors = cm.jet(W)

    '''
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111)
    xmin = 0 - m.pi
    xmax = m.pi + m.pi
    ymin = -m.pi - 2*m.pi
    ymax = m.pi + 2*m.pi

    ax1.imshow(np.rot90(W), cmap=plt.cm.gist_earth_r, extent=[xmin, xmax, ymin, ymax])
    ax1.plot(real_us[:,0], real_us[:,1], 'k.', markersize=2)
    ax1.set_xlim([xmin, xmax])
    ax1.set_ylim([ymin, ymax])
    plt.show()
    sys.exit(1)
    '''


    a = Arrow3D([-1.3,1.3],[0,0],[0,0], mutation_scale=20, lw=5, arrowstyle="-|>", color="g")
    ax.add_artist(a)
    surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, facecolors=colors,
            linewidth=0, antialiased=False)


    ax._axis3don=False
    ax.set_zlim3d(-1, 1)
    ax.w_zaxis.set_major_locator(LinearLocator(6))
    ax.view_init(0, options.angle)

    '''
    plt.subplots_adjust(left=0.4, right=0.9, top=0.9, bottom=0.1)

    for i in xrange(0, 360, 40):
        savefig("fig%d.png", (i))
    '''

    '''
    sm = cm.ScalarMappable(cmap=cm.jet)
    sm.set_array(W)
    fig.colorbar(sm)
    '''

    if options.fig_name != "":
        plt.savefig(options.fig_name, bbox_inches='tight')
    else:
        plt.show()

if __name__ == '__main__':
    main()

