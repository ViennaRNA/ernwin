#!/usr/bin/python
from __future__ import absolute_import, division, print_function, unicode_literals
from builtins import (ascii, bytes, chr, dict, filter, hex, input, #pip install future
                      int, map, next, oct, open, pow, range, round,
                      str, super, zip)
from future.utils import with_metaclass
from abc import ABCMeta, abstractmethod, abstractproperty
import numpy as np
import scipy.stats
import warnings
import logging
import itertools as it
import time
import math
import os

from commandline_parsable import parsable_base
from logging_exceptions import log_to_exception

from ..utils import get_version_string
from six.moves import map
import six


log = logging.getLogger(__name__)


DEFAULT_ENERGY_PREFACTOR = 30
INCR = 0.01

@parsable_base(False, required_kwargs=["cg"], factory_function="from_cg",
               name_attr="_shortname", helptext_sep="\n", help_attr="HELPTEXT",
               allow_pre_and_post_number=True, help_intro_list_sep="\n")
class EnergyFunction(with_metaclass(ABCMeta, object)):
    '''
    The base class for energy functions.
    '''
    __metaclass__ = ABCMeta

    @classmethod
    def from_cg(cls, prefactor, adjustment, cg, **kwargs):
        """
        Factory function. Return this energy for the given cg.

        :returns: An instance if this class or a CombinedEnergy containing instances of this class.
        """
        return cls(prefactor, adjustment)

    def __init__(self, prefactor=None, adjustment=None):
        self.log = logging.getLogger(self.__class__.__module__+"."+self.__class__.__name__)
        if prefactor is None:
            prefactor = DEFAULT_ENERGY_PREFACTOR
        if adjustment is None:
            adjustment = 1.
        #: Used by constraint energies, to store tuples of stems that clash.
        #: Updated every time eval_energy is called.
        self.bad_bulges = []

        #: The measure encountered during last energy evaluation
        self._last_measure = None
        #: The reference distribution.
        #: In the case of EnergyFunctions that are not CoarseGrainEnergy instances,
        #: this is only used to dump the measures to a file.
        self.accepted_measures = []

        #: The energy function can be adjusted with a prefactor (weight)
        #: and an adjustment (offset from the target value)
        self.prefactor, self._pf_stepwidth, self._pf_update_freq = self._parse_prefactor(prefactor, DEFAULT_ENERGY_PREFACTOR)
        self.adjustment, self._adj_stepwidth, self._adj_update_freq = self._parse_prefactor(adjustment, 1.0)

        #: How many sampling steps have been performed
        #: (Used for reference ratio method and simulated annealing)
        self.step=0

        #: Name and shortname of the energy
        # We need to check the class, not the instance, because implementation of name as a property
        # in subclasses may raise an error on incompletely initialized instanzes.
        if not hasattr(type(self), "name"):
            self.name = self.__class__.__name__.lower()

    @property
    def last_accepted_measure(self):
        return self.accepted_measures[-1]

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
        during sampling. Thus the last accepted measure should be added a second
        time to the reference distribution. (Reference ratio method)
        """
        if len(self.accepted_measures) > 0:
            self.accepted_measures.append(self.accepted_measures[-1])
        self._step_complete()

    @abstractproperty #!Note: Can be overwritten by a simple class-level variable, does not have to be a property.
    def _shortname(self):
        """
        A shortcut for the energy name. Ideally 3 to 4 letters, all-caps.

        This should be implemented by a simple class-level variable.

        It is used to generate energies from commandline-options
        (via fess.builder.energy.energies_from_string).
        Further more, the `shortname` property `_shortname` it for the description of the energy.
        """
        raise NotImplementedError

    @abstractproperty
    def HELPTEXT(self):
        """
        A preformatted helptext describing the energy.

        This is used by get_argparse_help and intended for
        output on the commandline via argparse's `--help` option.
        """
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
    def eval_energy(self, cg, background=True, nodes=None, **kwargs):
        raise NotImplementedError

    def dump_measures(self, base_directory, iteration=None):
        '''
        Dump all of the accepted measures collected so far
        to a file.
        '''
        out = " ".join(map("{:.4f}".format,self.accepted_measures))+"\n"
        output_file = op.join(base_directory, self.name+"_"+str(hex(id(self)))+".measures")
        with open(output_file, 'w') as f:
            f.write(out)

        if iteration is not None:
            with open(output_file + ".{:d}".format(iteration), 'w') as f:
                f.write(out)

    def _parse_prefactor(self, value, default):
        """
        Convert prefactor and adjustment to a 3-tuple

        :param value: A number, tuple or a string with either a number or 3 '_'-separated numbers.
        :returns: A tuple prefactor, stepwidth, update_frequency
        """
        if not value:
            return default, 0, float("inf")
        try:
            if isinstance(value[0], (float, int)): # A tuple of numbers:
                return value # Return the tuple unmodified. Unpacking the tuple will check its length
        except TypeError:
            return value,  0, float("inf") # Value is a number
        # Now value is a string or not a legal value.
        if "_" in value:
            p, s, f = value.split("_")
            return float(p), float(s), int(f)
        else:
            return float(value),  0, float("inf")

    def _step_complete(self):
        """
        Bookkeeping after every performed sampling step (NOT energy evaluation).

        Called by accept/reject last measure.
        If required, update "temperature" in simulated annealing simulations.
        """
        log.debug("step complete called on %s instance at %s.",type(self).__name__, hex(id(self)))
        self.step+=1
        if self.step % self._pf_update_freq <1 and self.step>0:
            self._update_pf()
        log.debug("step: %s  %% _adj_update_freq: %s = %s", self.step, self._adj_update_freq, self.step % self._adj_update_freq)
        if self.step % self._adj_update_freq < 1 and self.step>0:
            self._update_adj()
    def _update_pf(self):
        self.prefactor += self._pf_stepwidth
    def _update_adj(self):
        log.debug("Updating adjustment from %s to %s", self.adjustment, self.adjustment+self._adj_stepwidth)
        self.adjustment += self._adj_stepwidth

    def __bool__(self):
        """
        Returns always True, because a normal energy is never empty.
        This is implemented here to provide compatibility of the
        interface of Energy with CombinedEnergy
        """
        return True

    __nonzero__=__bool__


class InteractionEnergy(EnergyFunction):
    """
    A base-class for A-Minor and SLD energy. Uses the reference ratio method
    """
    def __init__(self, rna_length, num_loops, prefactor, adjustment):
        self.num_loops = num_loops
        super(InteractionEnergy, self).__init__(prefactor, adjustment)
        if num_loops == 0:
            log.warning("The number of loops is 0 for %s", self.shortname)
        self.reset_distributions(rna_length)

    @classmethod
    def qualifying_loops(cls, cg, loop_iterator):
        logger = logging.getLogger(cls.__module__+"."+cls.__name__)
        for l in loop_iterator:
            if l in cg.incomplete_elements or l in cg.interacting_elements:
                logger.debug("Loop %s does not qualify because it is interacting or incomplete.", l)
                logger.debug("Interacting %s.", cg.interacting_elements)
                logger.debug("Incomplete %s.", cg.incomplete_elements)
                continue
            yield l


    @abstractmethod
    def reset_distributions(self, rna_length):
        """
        This should set the following attributes:

        self.reference_interactions: A list with the fraction of loops that interact per entry.
                                     Each entry must correspond to self.num_loops loops.
        self.target_interactions: The fraction of loops that interact in the PDB
        """
        raise NotImplementedError

    @abstractmethod
    def _get_cg_measure(self, cg):
        """
        Return the fraction of loops that interact, i.e. the interaction-count/self.num_loops
        """
        raise NotImplementedError

    def eval_energy(self, cg, background=True, nodes=None, use_accepted_measure=False,
                    plot_debug=False, **kwargs):
        if plot_debug or nodes is not None:
            raise NotImplementedError("'plot_debug' and 'nodes' args are not implemented"
                                      " for InteractionEnergies")
        if use_accepted_measure:
            m = self.accepted_measures[-1]
        else:
            m = self._get_cg_measure(cg)
        self._last_measure = m

        reference_perc = sum(self.reference_interactions)/len(self.reference_interactions)
        target_perc    = self.target_interactions
        if background:
            e_i = -np.log( target_perc) + np.log(reference_perc)
            e_noi = -np.log( 1-target_perc) + np.log(1-reference_perc)
        else:
            e_i=-np.log(target_perc)
            e_noi=-np.log(1-target_perc)
        self.log.debug("%s , %s", e_i, e_noi)
        energy_contrib = e_i-e_noi

        num_interactions = int(round(m*self.num_loops))
        self.log.debug("Energy contribution of %s = %s, "
                  "with %s interactions (background=%s)",
                  self.shortname,
                  energy_contrib, num_interactions, background)
        energy=energy_contrib*num_interactions
        return self.prefactor*energy


    def _set_target_distribution(self):
        """
        Not needed, since the adjustment is multiplied to the
        target value in every eval_energy call
        """
        pass

    def _resample_background_kde(self):
        self.reference_interactions = self.accepted_measures[:]


class CoarseGrainEnergy(EnergyFunction):
    """
    A base-class for Energy functions that use a background distribution.
    """
    #: Change this to anything but "kde" to use a beta distribution (UNTESTED).
    dist_type = "kde"

    @classmethod
    def from_cg(cls, prefactor, adjustment, cg, **kwargs):
        """
        Factory function. Return this energy for the given cg.

        :returns: An instance if this class or a CombinedEnergy which is empty or contains instances of this class.
        """
        return cls(rna_length = cg.seq_length, prefactor=prefactor, adjustment=adjustment)

    def __init__(self, rna_length, prefactor=None, adjustment=None):
        super(CoarseGrainEnergy, self).__init__(prefactor, adjustment)

        self.reset_distributions(rna_length)

        #: The previous evaluated energy
        self.prev_energy = None

        #: Resample the reference distribution every n steps
        self.kde_resampling_frequency = 3


    @abstractproperty
    def sampled_stats_fn(self):
        """
        Can be overridden by class-level variable.

        A filename containing the target_distribution or None
        """
        raise NotImplementedError
    @abstractproperty
    def HELPTEXT(self):
        """
        Can be overridden by class-level variable.

        A filename containing the reference_distribution
        """
        raise NotImplementedError


    @classmethod
    @abstractmethod
    def generate_target_distribution(cls, *args, **kwargs):
        """
        Provide an implementation how the target distribution (and all other required files)
        can be generated. Document this method well to ensure reproducibility.

        This must be implemented as a classmethod, so changes to the class will
        be reflected in the target distribution.
        """
        raise NotImplementedError

    @staticmethod
    def _print_file_header(file_, cg_filenames):
        """
        Print a header to files generated with `generate_target_distribution`.

        This should allow for easy reproducibility where possible.

        :param file_: A file opened for writing
        :param cg_filenames: A list of strings (filenames)
        """
        print("# Generated on "+ time.strftime("%d %b %Y %H:%M:%S %Z"), file=file_)
        print("# Generated from the following files: ", file=file_)
        for fn in cg_filenames:
            print("#   {}".format(fn), file=file_)
        print("# Working directory: {}".format(os.getcwd()), file=file_)
        print("# Version: {}".format(get_version_string().strip()), file=file_)

    def plot_distributions(self, from_=None, to_=None, val = None):
        """
        :param from_: Minimal value of x-axis
        :param to_:   Maximal value of x-axis
        """
        import matplotlib.pyplot as plt

        if from_ is None:
            try:
                from_ = min(it.chain(self.accepted_measures, self.target_values))
            except AttributeError:
                from_ = min(self.accepted_measures)*0.8
            except ValueError:
                from_ = 0
        if to_ is None:
            try:
                to_ = max(it.chain(self.accepted_measures, self.target_values))
            except AttributeError:
                to_ = max(self.accepted_measures)*1.2
            except ValueError:
                to_ = 10
        xs=np.linspace(from_, to_, 2000)
        fig,ax1 = plt.subplots()
        ax2 = ax1.twinx()
        ax1.plot(xs, self.reference_distribution(xs), label="sampled")
        ax1.plot(xs, self.target_distribution(xs), label="target")
        ax2.plot(xs, -(np.log(self.target_distribution(xs) + 0.00000001 *
                     self.reference_distribution(xs)) - np.log(self.reference_distribution(xs))),
                 label="energy", color="red")
        plt.title(self.shortname)
        ax1.legend(loc="lower left")
        ax2.legend()
        if val is not None:
            ax1.plot([val,val], [0, max(self.target_distribution(xs))], "g-")

        plt.show(block = False)

    def reset_distributions(self, rna_length):
        """
        Reset the reference and target distribution to the values from files.

        :param rna_length: The length of the RNA in nucleotides.
        """
        if self.sampled_stats_fn is not None:
            log.debug("Loading sapmled measures into accepted_measures")
            self.accepted_measures = list(self._get_values_from_file(self.sampled_stats_fn, rna_length))
        #If sampled_stats_fn is None, we assume accepted_measures is given in the constructor
        if self.accepted_measures:
            self.reference_distribution = self._get_distribution_from_values(self.accepted_measures)
        else:
            raise ValueError("Either sampled_stats_fn or accepted_measures has to be set "
                             "before calling CoarseGrainEnergy.__init__ or "
                             "CoarseGrainEnergy.reset_distribution "
                             "for {}.".format(type(self).__name__))
        if self.real_stats_fn is not None:
            self.target_values =  self._get_values_from_file(self.real_stats_fn, rna_length)
        self._set_target_distribution()

    def _step_complete(self):
        """
        Call superclass _step_complete and resample background_kde every n steps.
        """
        log.info("CoarseGrainedEnergy step complete")
        super(CoarseGrainEnergy, self)._step_complete()
        if self.step % self.kde_resampling_frequency == 0:
            self._resample_background_kde()

    def _resample_background_kde(self):
        """
        Update the reference distribution based on the accepted values
        """
        log.debug("Resampling background KDE for %s. Now %d accepted measures", type(self).__name__, len(self.accepted_measures))
        values = self.accepted_measures
        new_kde = self._get_distribution_from_values(values)
        if new_kde is not None:
            self.reference_distribution = new_kde
            log.debug("Density of ref AFTER resampling = %s", self.reference_distribution(self.accepted_measures[-1]))
        else:
            log.warning("Distribution is None. Cannot change background_kde")

    @classmethod
    @abstractmethod
    def _get_values_from_file(cls, filename, nt_length):
        raise NotImplementedError

    @staticmethod
    def _values_within_nt_range(data, length, target_col, length_col="nt_length", target_len=500):
        rdata = []

        distribution_lower_bound = 1.
        distribution_upper_bound = 1.

        while (len(rdata) < target_len and len(rdata)<len(data)):
            try:
                distribution_lower_bound -= INCR
                distribution_upper_bound += INCR
                rdata = data[data[length_col] > distribution_lower_bound * length]
                rdata = rdata[rdata[length_col] < length * distribution_upper_bound]
            except KeyboardInterrupt:
                log.error("len(rdata) is {}, len(data)={}, bound= {}...{}".format(len(rdata),len(data),
                     distribution_lower_bound ) * length, length * ( distribution_upper_bound ))
                raise
        if len(rdata)==0:
            raise ValueError("No data found for distribution")
        log.info("%d datapoints", len(rdata))
        return rdata[target_col]

    @classmethod
    def _get_distribution_from_values(cls, values):
        '''
        Return a probability distribution from the given values.

        :param values: A list of values to fit a distribution to.
        :return: A probability distribution fit to the values.
        '''

        log.debug("Getting distribution from values of shape {}".format(np.shape(values)))
        if cls.dist_type == "kde":
            try:
                k = scipy.stats.gaussian_kde(values)
            except np.linalg.linalg.LinAlgError:
                log.exception("Setting KDE for %s to None because of", values)
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

    def eval_energy(self, cg, background=True, nodes=None, use_accepted_measure=False, plot_debug=False, **kwargs):
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
            self.plot_distributions(val=m)
        if self.log.isEnabledFor(logging.DEBUG):
            try:
                self.log.debug("Measure is {:1.4f}".format(m))
            except ValueError:
                self.log.debug("Measure is {}".format(m))


        self._last_measure = m

        if background:
            ref_val = self.reference_distribution(m)
            tar_val = self.target_distribution(m)
            if len(ref_val)>1:
                self.log.debug("ref_values: %s", ref_val)
                self.log.debug("tar_values: %s", tar_val)
                ref_val=np.maximum(ref_val,0)
                tar_val=np.maximum(tar_val,0)
                energies = np.log( tar_val ) - np.log(ref_val)
                energies = np.minimum(energies, 10**8)
                energies = np.maximum(energies, -10**8)
                energy = sum(energies)
                if np.isinf(energy):
                    log.warning(energies)
                    for i,r in enumerate(ref_val):
                        if r<=0 or tar_val[i]<=0:
                            log.warning("%sth measure %s has prob. %s in reference distribution and %s in taret distribution", i, m[i],r, tar_val[i])
            else:
                energy, = (np.log( tar_val ) - np.log(ref_val))
            self.log.debug("Energy (not yet scaled) = {}".format(energy))
            self.prev_energy = energy
            return -1 * self.prefactor * energy
        else:
            l = np.log(self.target_distribution(m))
            self.log.debug("Energy, = {}".format(l))
            if len(l)>1:
                energy = sum(l)
            else:
                energy, = l
            return -energy

    def _update_adj(self):
        super(CoarseGrainEnergy, self)._update_adj()
        self._set_target_distribution()

    def _set_target_distribution(self):
        log.info("Adjusting target distribution (base class)")
        scaled_vals = np.asarray(self.target_values)*self.adjustment
        self.target_distribution = self._get_distribution_from_values(scaled_vals)
