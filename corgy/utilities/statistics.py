import scipy.stats as ss
import scipy.interpolate as si
import scipy.optimize as so

import numpy as np
import scipy.special as sspec

import numpy as np

def my_log(x):
    return np.log(x + 1e-200)


def pdf(x):
    return 1/np.sqrt(2 * np.pi) * np.exp(-x ** 2 / 2)

def cdf(x):
    return (1 + sspec.erf(x / np.sqrt(2))) / 2

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
    kde = ss.gaussian_kde(x)
    fzz = kde(x)

    #normalize the data
    optm = lambda l, x: skew(x, l[0], l[1], l[2]) - fzz 
    fit = so.leastsq(optm, [np.mean(x), np.std(x), 0.5], (x,))[0]

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
    for x in np.linspace(-10, 10, n):
        mult = int(100 * skew(x, e, w, a))
        dataset += [x + ss.norm.rvs(0, 0.1) for i in range(mult)]

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
        self.kde = ss.gaussian_kde(data)
        y = my_log(self.kde(x))

        # the data needs to be sorted for the interpolation to succeed
        d = zip(x, y)
        d.sort()
        x,y = zip(*d)

        self.inter = si.splrep(x, y)

    def evaluate(self, data):
        '''
        Get the probability density of the data as if we were
        using the kde.

        @param data: An array-like set of data.
        '''
        return si.splev(data, self.inter)

    def __call__(self, data):
        return self.evaluate(data)

def meshgrid2(*arrs):
    '''
    Snatched from:
    
    http://stackoverflow.com/questions/1827489/numpy-meshgrid-in-3d.
    '''
    arrs = tuple(reversed(arrs))  #edit
    lens = map(len, arrs)
    dim = len(arrs)

    sz = 1
    for s in lens:
        sz*=s

    ans = []    
    for i, arr in enumerate(arrs):
        slc = [1]*dim
        slc[i] = lens[i]
        arr2 = np.asarray(arr).reshape(slc)
        for j, sz in enumerate(lens):
            if j!=i:
                arr2 = arr2.repeat(sz, axis=j) 
        ans.append(arr2)

    return tuple(ans[::-1])

class InterpolatedMultiKDE:
    '''
    An interpolated kernel density estimate.
    '''

    def __init__(self, data, grain=None):
        '''
        Create a kde from the data and then fit a cubic spline to it.

        @param grain: Use n points to create the grid over which to interpolate.
                      If it's not passed in, then the number of points is equal
                      to the length of the data set.
        '''
        '''
        grid_arrays = []

        if grain == None:
            skip_points = 1
        else:
            skip_points = len(data) / grain

        if len(data.shape) == 1:
            # 1-d array
            grid_arrays += [data[0::skip_points]]
        else:
            for i in xrange(data.shape[1]):
                grid_arrays += [data[:,i][0::skip_points]]

        grid = meshgrid2(*grid_arrays)
        grid_points = np.array(zip(*(x.flat for x in grid)))
        '''

        kde = ss.gaussian_kde(data.T)

        #self.interp = si.LinearNDInterpolator(grid_points, kde(grid_points.T))
        self.interp = si.LinearNDInterpolator(data, kde(data.T))


    def evaluate(self, data):
        '''
        Get the probability density of the data as if we were
        using the kde.

        @param data: An array-like set of data.
        '''
        return self.interp([data])

    def __call__(self, data):
        return self.evaluate(data)


