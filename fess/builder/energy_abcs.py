#!/usr/bin/python
from __future__ import absolute_import, division, print_function, unicode_literals
from builtins import (ascii, bytes, chr, dict, filter, hex, input, #pip install future
                      int, map, next, oct, open, pow, range, round,
                      str, super, zip)
from future.builtins.disabled import (apply, cmp, coerce, execfile,
                             file, long, raw_input, reduce, reload,
                             unicode, xrange, StandardError)
                             
from abc import ABCMeta, abstractmethod, abstractproperty
import numpy as np
import scipy.stats

import logging
log = logging.getLogger(__name__)


DEFAULT_ENERGY_PREFACTOR = 30

class EnergyFunction(object):
    '''
    The base class for energy functions.
    '''
    __metaclass__ = ABCMeta
    
    def __init__(self, prefactor = DEFAULT_ENERGY_PREFACTOR, adjustment = 1):
        
        #: Used by constraint energies, to store tuples of stems that clash.
        #: Updated every time eval_energy is called.
        self.bad_bulges = []

        #: The measure encountered during last energy evaluation
        self._last_measures = None
        #: The reference distribution. 
        #: In the case of EnergyFunctions that are not CoarseGrainEnergy instances, 
        #: this is only used to dump the measures to a file.
        self.accepted_measures = []
        
        #: The energy function can be adjusted with a prefactor (weight) 
        #: and an adjustment (offset from the target value)
        self.prefactor, self._pf_stepwidth, self._pf_update_freq = self._parse_prefactor(prefactor)
        self.adjustment, self._adj_stepwidth, self._adj_update_freq = self._parse_prefactor(adjustment)
        
        #: How many sampling steps have been performed 
        #: (Used for reference ratio method and simulated annealing)
        self.step=0
        
        #: Name and shortname of the energy 
        self.name = self.__class__.__name__.lower()
        
    def accept_last_measure(self):
        """
        The last measure should contribute to the new reference distribution.
        
        The structure evaluated with the last call to eval_energy was accepted
        during sampling. Thus the corresponding measure should be included in 
        the reference distribution. (Reference ratio method)
        """
        if self._last_measure is not None:        
            self.accepted_measures.append(self._last_measure)
        self._step_complete()
   
    def reject_last_measure(self):
        """
        The last measure of the reference distribution should be duplicated.
        
        The structure evaluated with the last call to eval_energy was rejected
        during sampling. Thus the last accepted measure should be adefaultdictdded a second 
        time to the reference distribution. (Reference ratio method)
        """        
        if len(self.accepted_measures) > 0:
            self.accepted_measures.append(self.accepted_measures[-1])
        self._step_complete()

    @abstractproperty #!Note: Can be overwritten by a simple class-level variable, does not have to be a property.
    def _shortname(self):
        raise NotImplementedError

    @property
    def shortname(self):
        if self.prefactor==DEFAULT_ENERGY_PREFACTOR:
            pre=""
        else:
            pre=self.prefactor
        if self.adjustment==1.0:
            adj=""
        else:
            adj=self.adjustment
        return "{}{}{}".format(pre, self._shortname, adj)

        return 

    @abstractmethod
    def eval_energy(self, cg, background=True, nodes=None):
        raise NotImplementedError
    
    def dump_measures(self, base_directory, iteration=None):
        '''
        Dump all of the accepted measures collected so far
        to a file.
        '''
        output_file = op.join(base_directory, self.name+"_"+str(hex(id(self)))+".measures")
        with open(output_file, 'w') as f:
            f.write(" ".join(map("{:.4f}".format,self.accepted_measures)))
            f.write("\n")

        if iteration is not None:
            with open(output_file + ".%d" % (iteration), 'w') as f:
                f.write(" ".join(map("{:.4f}".format,self.accepted_measures)))
                f.write("\n")
    
    def _parse_prefactor(self, prefactor):
        """
        Creates a tuple from a scalar or returns a given tuple
        
        :param prefactor: A number or a 3-tuple
        :returns: A tuple prefactor, stepwidth, update_frequency
        """
        if isinstance(prefactor, tuple):
            return prefactor
        else:
            return prefactor, None, float("inf")

    def _step_complete(self):
        """
        Bookkeeping after every performed sampling step (NOT energy evaluation).
        
        Called by accept/reject last measure. 
        If required, update "temperature" in simulated annealing simulations.
        """
        self.step+=1
        if self.step % self._pf_update_freq == 0 and self.step>0:
            self._update_pf()
        if self.step % self._adj_update_freq == 0 and self.step>0:
            self._update_adj()
    def _update_pf(self):
        self.prefactor += self._pf_stepwidth
    def _update_adj(self):
        self.adjustment += self._adj_stepwidth

class CoarseGrainEnergy(EnergyFunction):
    """
    A base-class for Energy functions that use a background distribution.
    """
    def __init__(self, rna_length, prefactor, adjustment):
        super(CoarseGrainEnergy, self).__init__(prefactor, adjustment)
        
        #: Change this to anything but "kde" to use a beta distribution (UNTESTED).
        self.dist_type = "kde"

        self.reference_distribution, self.accepted_measures = self._get_distribution_from_file(self.sampled_stats_fn, rna_length)
        self.target_distribution, self.target_values =  self._get_distribution_from_file(self.real_stats_fn, rna_length)
        if adjustment!=1:
            self._adjust_target_distribution()

        #: The previous evaluated energy
        self.prev_energy = None
        
        #: Resample the reference distribution every n steps
        self.kde_resampling_frequency = 3
        
    def _step_complete(self):
        """
        Call superclass _step_complete and resample background_kde every n steps.
        """
        super(CoarseGrainEnergy, self)._step_complete()
        if self.step % self.kde_resampling_frequency == 0:
            self.resample_background_kde()
    
    def resample_background_kde(self):
        """
        Update the reference distribution based on the accepted values
        """
        values = self.accepted_measures
        if len(values) > 100:
            new_kde = self._get_distribution_from_values(values)
            if new_kde is not None:
                self.reference_distribution = new_kde
                
    @abstractmethod
    def _get_distribution_from_file(self):
        raise NotImplementedError
    def _get_distribution_from_values(self, values):
        '''
        Return a probability distribution from the given values.

        :param values: A list of values to fit a distribution to.
        :return: A probability distribution fit to the values.
        '''

        log.debug("Getting distribution from values of shape {}".format(np.shape(values)))
        if self.dist_type == "kde":
            try:
                k = scipy.stats.gaussian_kde(values)
            except np.linalg.linalg.LinAlgError:
                log.exception()
                return None
        else:        
            floc = -0.1
            fscale =  1.5 * max(values)
            f = scipy.stats.beta.fit(values, floc=floc, fscale=fscale)
            k = lambda x: scipy.stats.beta.pdf(x, f[0], f[1], f[2], f[3])
        return k
    @abstractmethod
    def _get_cg_measure(self, cg):
        raise NotImplementedError
    
    def eval_energy(self, cg, background=True, nodes=None, use_accepted_measure=False, plot_debug=False):
        '''
        A generic function which simply evaluates the energy based on the
        previously calculated probability distributions.
        
        :param use_accepted_measure: Do not calculate the measure from the cg. Instead use the
                                     last accepted measure.
        :param plot_debug: Plot the distributions for debugging
        '''
        if use_accepted_measure:
            m = self.accepted_measures[-1]
        else:
            m = self._get_cg_measure(cg)

        if plot_debug: #For debuging
            import matplotlib.pyplot as plt
            xs=np.linspace(0, 200, 2000)
            fig,ax1 = plt.subplots()
            ax2 = ax1.twinx()
            ax1.plot(xs, self.reference_distribution(xs), label="sampled")
            ax1.plot(xs, self.target_distribution(xs), label="reference")
            ax2.plot(xs, -(np.log(self.target_distribution(xs) + 0.00000001 * self.reference_distribution(xs)) - np.log(self.reference_distribution(xs))), label="energy", color="red")
            plt.title(self.shortname())
            ax1.legend(loc="lower left")
            ax2.legend()
            plt.show()

        log.debug("Measure is {:1.4f}".format(m))
        
        self._last_measure = m
        
        if background:
            ref_val = self.reference_distribution(m)
            tar_val = self.target_distribution(m)
            energy, = (np.log( tar_val + 0.00000001 * ref_val) - np.log(ref_val))
            self.prev_energy = energy
            return -1 * self.prefactor * energy
        else:
            l = np.log(self.target_distribution(m))
            log.debug("Energy, = {}".format(l))
            energy, = l
            return -energy * self.prefactor

    def _update_adj(self):
        super(RadiusOfGyrationEnergy, self)._update_adjustment()
        self._adjust_target_distribution()
        
    def _adjust_target_distribution(self):
        self.target_distribution = self._get_distribution_from_values(self.target_values*self.adjustment)
