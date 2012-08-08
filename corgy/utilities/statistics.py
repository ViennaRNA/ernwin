from scipy import linspace
from scipy import pi, sqrt, exp, mean, std
from scipy.special import erf
from scipy.stats import norm, gaussian_kde
from scipy import interpolate
from scipy.interpolate import interp1d

from scipy.optimize import leastsq

from numpy import log

def my_log(x):
    return log(x + 1e-200)


def pdf(x):
    return 1/sqrt(2 * pi) * exp(-x ** 2 / 2)

def cdf(x):
    return (1 + erf(x / sqrt(2))) / 2

def skew(x,e=0,w=1,a=0):
    t = (x - e) / w
    return 2 / w * pdf(t) * cdf(a * t)
    # You can of course use the scipy.stats.norm versions
    # return 2 * norm.pdf(t) * norm.cdf(a*t)

def fit_skew(x):
    '''
    Fit a skew-normal distribution to a set of data.

    @param x: The data set.
    @return: The parameters (e, w, a)
    '''

    #Get a kernel density estimate of the probability density
    kde = gaussian_kde(x)
    fzz = kde(x)

    #normalize the data
    optm = lambda l, x: skew(x, l[0], l[1], l[2]) - fzz 
    fit = leastsq(optm, [mean(x), std(x), 0.5], (x,))[0]

    '''
    hist(x)
    plot(x, len(x) * kde(x), 'o')
    plot(x, len(x) * skew(x,fit[0], fit[1], fit[2]), 'ro', linewidth=2)
    #show()
    '''

    return fit


def test_fit_skew():
    '''
    Test the fit_skew function by creating an artificial data set
    and fitting the distribution to it.
    '''
    n = 4000


    e = 0.0 # location
    w = 3.0 # scale
    a = 4.0

    dataset = []
    for x in linspace(-10, 10, n):
        mult = int(100 * skew(x, e, w, a))
        dataset += [x + norm.rvs(0, 0.1) for i in range(mult)]

    fit_skew(dataset)

#test_fit_skew()

class interpolated_kde:
    '''
    An interpolated kernel density estimate.
    '''

    def __init__(self, data, limits=None, n=5000):
        '''
        Create a kde from the data and then fit a cubic spline to it.
        '''
        if limits == None:
            limits = (min(0, min(data)), max(500, max(data)))

        x = linspace(limits[0], limits[1], n)
        self.kde = gaussian_kde(data)
        y = my_log(self.kde(x))

        # the data needs to be sorted for the interpolation to succeed
        d = zip(x, y)
        d.sort()
        x,y = zip(*d)

        self.inter = interpolate.splrep(x, y)

    def evaluate(self, data):
        '''
        Get the probability density of the data as if we were
        using the kde.

        @param data: An array-like set of data.
        '''
        return interpolate.splev(data, self.inter)

    def __call__(self, data):
        return self.evaluate(data)

