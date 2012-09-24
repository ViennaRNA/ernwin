import unittest

from corgy.utilities.statistics import interpolated_kde
from scipy.stats import norm
from scipy.stats import gaussian_kde
import matplotlib.pylab as pl
import numpy as np
from numpy import log, allclose

class TestInterpolatedKde(unittest.TestCase):

    def test_estimate(self):
        
        samples = norm.rvs(size=100)
        bsamples = np.linspace(-4, 4)

        kde = gaussian_kde(samples)
        ik = interpolated_kde(samples, (min(bsamples), max(bsamples)))

        kdata = log(kde(samples))
        idata = ik(samples)

        self.assertTrue(allclose(kdata, idata))

        '''
        pl.hist(samples, normed=True)
        #pl.plot(samples, kde(samples), 'o')
        pl.plot(bsamples, kde(bsamples), 'o')

        pl.plot(bsamples, ik(bsamples))
        pl.plot(bsamples, norm.pdf(bsamples))

        pl.show()
        '''
