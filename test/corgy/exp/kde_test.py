import unittest
import time
import timeit

import corgy.exp.kde as cek
import corgy.utilities.debug as cud

import numpy as np
import numpy.random as nr
import scipy.stats as ss
import scipy.ndimage as sn

import matplotlib.pyplot as plt

class TestKDE(unittest.TestCase):
    def setUp(self):
        print

    def test_3d(self):
        source_data = nr.rand(1000,3)

        xmin  = 0
        xmax = 1
        ymin = 0
        ymax = 1
        zmin = 0
        zmax = 0
        resolution = 60 

        X, Y, Z = np.mgrid[xmin:xmax:np.complex(0, resolution), ymin:ymax:np.complex(0, resolution), zmin:zmax:np.complex(0, resolution)]
        positions = np.vstack([X.ravel(), Y.ravel(), Z.ravel()])

        #cud.pv('source_data.shape')
        #cud.pv('source_data')

        ss_kde = ss.gaussian_kde(source_data.T)
        my_kde = cek.gaussian_kde(source_data.T)

        cud.pv('len(positions[0])')

        t = time.time()
        ss_kde(positions)
        print "ss_kde time:", time.time() - t

        t = time.time()
        my_kde(positions)
        print "my_kde time:", time.time() - t

    def test_2d(self):
        source_data = nr.rand(10,2)

        xmin  = 0
        xmax = 1
        ymin = 0
        ymax = 1
        resolution = 10 

        X, Y = np.mgrid[xmin:xmax:np.complex(0, resolution), ymin:ymax:np.complex(0, resolution)]

        positions = np.vstack([X.ravel(), Y.ravel()])

        ss_kde = ss.gaussian_kde(source_data.T)
        my_kde = cek.gaussian_kde(source_data.T)

        t = time.time()
        vals_ss = ss_kde(positions)
        print "ss_kde time:", time.time() - t

        t = time.time()
        vals_k = my_kde(positions)
        print "ss_kde time:", time.time() - t

        fig = plt.figure()
        ax = fig.add_subplot(121)

        ax.imshow(np.rot90(vals_ss.reshape(X.shape)), extent=(xmin,xmax,ymin,ymax), origin='lower')
        ax.plot(source_data.T[0], source_data.T[1], 'k.', markersize=2)
        ax.set_xlim([xmin, xmax])
        ax.set_ylim([ymin, ymax])

        ax = fig.add_subplot(122)

        ax.imshow(np.rot90(vals_k.reshape(X.shape)), extent=(xmin,xmax,ymin,ymax), origin='lower')
        ax.plot(source_data.T[0], source_data.T[1], 'k.', markersize=2)
        ax.set_xlim([xmin, xmax])
        ax.set_ylim([ymin, ymax])

        #plt.show()

    def test_ndimage(self):
        source_data = nr.rand(200,2)

        xmin  = 0
        xmax = 1
        ymin = 0
        ymax = 1
        resolution = 30 

        X, Y = np.mgrid[xmin:xmax:np.complex(0, resolution), ymin:ymax:np.complex(0, resolution)]
        positions = np.vstack([X.ravel(), Y.ravel()])

        ss_kde = ss.gaussian_kde(source_data.T)
        vals_ss = ss_kde(positions)

        img = np.zeros((resolution, resolution))

        xincr = (xmax - xmin) / float(resolution)
        yincr = (ymax - ymin) / float(resolution)

        for d in source_data:
            img[int(d[0] / xincr), int(d[1] / yincr)] += 1
        img = sn.gaussian_filter(img, (1,1))

        fig = plt.figure()
        ax = fig.add_subplot(121)
        ax.imshow(img.T, extent=(xmin, xmax, ymin, ymax), origin='lower')
        ax.plot(source_data.T[0], source_data.T[1], 'k.', markersize=3)
        ax.set_xlim([xmin, xmax])
        ax.set_xlim([ymin, ymax])

        ax = fig.add_subplot(122)
        ax.imshow(np.rot90(vals_ss.reshape(X.shape)), extent=(xmin,xmax,ymin,ymax), origin='lower')
        ax.plot(source_data.T[0], source_data.T[1], 'k.', markersize=2)
        ax.set_xlim([xmin, xmax])
        ax.set_ylim([ymin, ymax])

        #plt.show()

        
        
