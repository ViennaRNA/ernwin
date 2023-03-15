from __future__ import absolute_import, division, print_function, unicode_literals
from builtins import (ascii, bytes, chr, dict, filter, hex, input,
                      int, map, next, oct, open, pow, range, round,
                      str, super, zip)
# from future.builtins.disabled import (apply, cmp, coerce, execfile,
#                             file, long, raw_input, reduce, reload,
#                             unicode, xrange, StandardError)

import sys
import math
import time
import copy
import os.path
import warnings
import logging
from collections import defaultdict
import contextlib
import itertools as it

import numpy as np
import pandas as pd

from logging_exceptions import log_to_exception, log_exception

import forgi.threedee.model.descriptors as ftur
import forgi.threedee.model.similarity as ftme
import forgi.threedee.utilities.vector as ftuv
import forgi.threedee.utilities.pdb as ftup
import forgi.utilities.commandline_utils as fuc
import fess.builder.reconstructor as fbr

from . import config as conf
from . import energy as fbe
from ..SortedCollection import SortedCollection
from six.moves import map
from six.moves import range

log = logging.getLogger(__name__)
__metaclass__ = type


class StatisticsCollector(object):
    """
    A base class for Spatial-Model Statistic objects.
    """
    header = []  # This is used for outfile-parsing.

    def __init__(self, *args, **kwargs):
        """
        :var self.header: A list of strings to describe the fields of history collected.
        :var self.history: A list of lists.
        :var self.silent: This is set to true, if a Collector cannot initialize itself properly.
        :var self.display: If this is True, the collector will write to stdout.
        """
        self.history = [[]]
        self.silent = False
        self.display = True

    @property
    def header_str(self):
        """
        Return the header as a tab-separated string.
        """
        return "\t".join(self.header)

    def update(self, sm, step):
        """
        Add statistics for the given spatial model to the history
        and return them as a string ready for printing to the screen.
        """
        raise NotImplementedError("Needs to be subclassed")

    @staticmethod
    def parse_value(stri):
        return stri


class CombinedStatistics(StatisticsCollector):
    def __init__(self, collectors, separator=""):
        """
        Used as a convenient way to combine several
        statisticsCollector instances into one object.

        :param collectors: A list of StatisticsCollector instances.
        :param separator: If given, this character will be used to delimit the individual entries.
        """
        self._members = collectors
        if separator:
            self._joiner = "\t"+separator+"\t"
        else:
            self._joiner = "\t"

    @property
    def header_str(self):
        return self._joiner.join(member.header_str for member in self._members if not member.silent)

    @property
    def header(self):
        """
        This only includes headers for which a history is written.
        """
        return [header for member in self._members for header in member.header_list
                if member.history is not None and not member.silent]

    @property
    def history(self):
        hist = []
        for member in self._members:
            if member.history is None or member.silent:
                continue
            hist += member.history
        return hist
        # return [ history for member in self._members for history in member.history if not member.silent and member.history is not None ]

    def update(self, sm, step):
        line = []
        for member in self._members:
            if not member.silent:
                log.debug("Updating {}".format(type(member)))
                fields = member.update(sm, step)
                if member.display:
                    line.append(fields)
        return self._joiner.join(line)

    def to_file(self):
        history = self.history
        length = len(history[0])
        j = None
        with open(os.path.join(conf.Configuration.sampling_output_dir, "monitor.txt"), "w") as f:
            print("#"+self._joiner.join(h for member in self._members for h in member.header if not member.silent and member.history is not None), file=f)
            for i in range(length):
                try:
                    line = []
                    for h in history:
                        if isinstance(h[i], list):
                            if self._joiner == ",":
                                j = " "
                            else:
                                j = ","
                            line.append(j.join(map(str, h[i])))
                        else:
                            line.append(str(h[i]))
                except Exception as e:
                    log.error("Members are: {}".format(
                        ", ".join(map(str, self._members))))
                    log.exception("Could not dump complete history to file 'monitor.txt': {} occurred for history column {} line {}".format(
                        type(e).__name__, j, i))
                    return
                else:
                    print(self._joiner.join(line), file=f)
        for member in self._members:
            if hasattr(member, 'to_file'):
                member.to_file()


class ROGStatistics(StatisticsCollector):
    """
    Store and print the Radius of Gyration.
    """
    header = ["ROG"]

    def update(self, sm, step):
        rog = ftur.radius_of_gyration(sm.bg.get_ordered_stem_poss())
        self.history[0].append(rog)
        return "{:6.2f} A".format(rog)

    @staticmethod
    def parse_value(stri):
        stri = stri.split()
        assert stri[1] == "A", stri
        return float(stri[0])


class AsphericityStatistics(StatisticsCollector):
    """
    Store and print the Asphericity.
    """
    header = ["Asphericity"]

    def update(self, sm, step):
        rog = ftur.asphericity(sm.bg.get_ordered_stem_poss())
        self.history[0].append(rog)
        return "{:6.2f}".format(rog)

    @staticmethod
    def parse_value(stri):
        return float(stri)


class AnisotropyStatistics(StatisticsCollector):
    """
    Store and print the Anisotropy.
    """
    header = ["Anisotropy"]

    def update(self, sm, step):
        rog = ftur.anisotropy(sm.bg.get_ordered_stem_poss())
        self.history[0].append(rog)
        return "{:6.2f}".format(rog)

    @staticmethod
    def parse_value(stri):
        return float(stri)


class ACCStatistics(StatisticsCollector):
    """
    Store and print the Adjacency Correlation Coefficient
    """
    header = ["ACC"]

    def __init__(self, reference_cg):
        super(ACCStatistics, self).__init__()
        try:
            self._cm_calc = ftme.AdjacencyCorrelation(reference_cg)
        except Exception as e:
            log.exception(e)
            warnings.warn(
                "Cannot report ACC. {} in reference SM: {}".format(type(e), str(e)))
            self.silent = True
            self.history = None

    def update(self, sm, step):
        try:
            if self.silent:
                return
            else:
                acc = ftme.mcc(self._cm_calc.evaluate(sm.bg))
                self.history[0].append(acc)
            return "{:5.3f}".format(acc)
        except ZeroDivisionError:
            self.history[0].append(float("nan"))
            return "{:5.3f}".format(float("nan"))

    @staticmethod
    def parse_value(stri):
        return float(stri)

class StemRMSD(StatisticsCollector):
    """
    Store and print stem-RMSD to multiple references
    """
    header = ["stem_rmsd"]

    def __init__(self, ref_fn, name, args):
        """
        :param ref_fn: A file name of an RNA.
        """
        super(StemRMSD, self).__init__()
        try:
            args.rna = [ref_fn]
            reference_cg, = fuc.cgs_from_args(args, rna_type="cg", enable_logging=False)
        except Exception as e:
            log.exception("Cannot create StemRMSD.")
            self.silent=True
            return
        self._reference = reference_cg
        if name:
            self.header = ["stem_rmsd_of_"+name]
        self.has_logged=False

    def update(self, sm, step):
        try:
            srmsd = ftme.cg_stem_rmsd(sm.bg, self._reference)
        except Exception as e:
            if not self.has_logged:
                log.exception("Cannot calculate StemRMSD.")
                self.has_logged=True
            srmsd=float("nan")
        self.history[0].append(srmsd)
        return "{:6.3f}".format(srmsd)

    @staticmethod
    def parse_value(stri):
        return float(stri)

class RMSDStatistics(StatisticsCollector):
    """
    Store and print the Root Mean Square Deviation from the reference SM.
    """
    header = ["RMSD", "dRMSD", "maxRMSD",
              "maxdRMSD"]  # All allowed headers. __init__ sets self.header, which takes priority.

    def __init__(self, reference_cg, show_min_max=True, save_n_best=0, mode="RMSD"):
        """
        :param reference_sm: The reference spatial model, against which to collect statistics.
        :param show_min_max: Bool. Whether to show the minimal and maximal RMSD so far.
        :param mode: "RMSD" or "dRMSD"
        """
        super(RMSDStatistics, self).__init__()
        self.best_cgs = SortedCollection(
            key=lambda x: x[1], maxlen=save_n_best)
        self.mode = mode
        self.reference_cg = reference_cg
        self.res_ids=None
        try:
            self._reference = reference_cg.get_ordered_virtual_residue_poss()
        except:
            self.silent = True
            self.history = None
        else:
            self._showMinMax = show_min_max
            if show_min_max == "max":
                self.header = [mode, "max"+mode]
                self.history = [[], []]
                self._maxRMSD = float("-inf")
            elif show_min_max:
                self.header = [mode, "min"+mode, "max"+mode]
                self.history = [[], [], []]
                self._minRMSD = float("inf")
                self._maxRMSD = float("-inf")
            else:
                self.header = [mode]

    def update(self, sm, step):
        if not self.silent:
            if self.res_ids:
                curr_vress=[]
                for res in self.res_ids:
                    curr_vress.append(sm.bg.get_virtual_residue(res, allow_single_stranded = True))
                curr_vress = np.array(curr_vress)
            else:
                curr_vress = sm.bg.get_ordered_virtual_residue_poss()
            if self.mode == "RMSD":
                try:
                    rmsd = ftme.rmsd(self._reference, curr_vress)
                except Exception as e:
                    self.res_ids = set(self.reference_cg.seq._seqids) & set(sm.bg.seq._seqids)
                    self._reference = []
                    curr_vress=[]
                    for res in self.res_ids:
                        curr_vress.append(sm.bg.get_virtual_residue(res, allow_single_stranded = True))
                        self._reference.append(self.reference_cg.get_virtual_residue(res, allow_single_stranded = True))
                    curr_vress = np.array(curr_vress)
                    self._reference = np.array(self._reference)
                    log.error("Calculating rmsd between %s and %s", self._reference, curr_vress)
                    rmsd = ftme.rmsd(self._reference, curr_vress)
            elif self.mode == "dRMSD":
                rmsd = ftme.drmsd(self._reference, curr_vress)
            if self.best_cgs.can_insert_right((None, rmsd)):
                self.best_cgs.insert_right((sm.bg.to_cg_string(), rmsd))
            if step % 10 == 0:
                for i, bg in enumerate(self.best_cgs):
                    with open(os.path.join(conf.Configuration.sampling_output_dir,
                                           'best_rmsd{:d}.coord'.format(i)), 'w') as f:
                        f.write(bg[0])
            if self._showMinMax:
                if rmsd > self._maxRMSD:
                    self._maxRMSD = rmsd
                if self._showMinMax == True:
                    if rmsd < self._minRMSD:
                        self._minRMSD = rmsd
                    self.history[1].append(self._minRMSD)
                self.history[0].append(rmsd)
                self.history[-1].append(self._maxRMSD)
                if self._showMinMax == "max":
                    return "{:6.3f} A\t{:6.3f} A".format(rmsd, self._maxRMSD)
                else:
                    return "{:6.3f} A\t{:6.3f} A\t{:6.3f} A".format(rmsd, self._minRMSD, self._maxRMSD)
            else:
                self.history[0].append(rmsd)
                return "{:6.3f} A".format(rmsd)
        else:
            return

    @staticmethod
    def parse_value(stri):
        stri = stri.split()
        assert stri[1] == "A"
        return float(stri[0])


class EnergyTracking(StatisticsCollector):
    """
    After every step, evaluate an energy not used for sampling
    """
    header = ["Tracked Energy"]

    def __init__(self, energy_function, background=False):
        """
        :param energy_function: The energy function which will be evaluated in every update step.
        :param background: Bool. If true, call eval_energy with background=True and
              call accept_last_measure().
             .. warning:: If background is set to True, this class will call accept_last_measure()
                          of the energy function, which might have side-effects if a reference to
                          the energy_function is used elsewhere.
        """
        super(EnergyTracking, self).__init__()
        self._energy_function = energy_function
        self._background = background
        self.header = ["Energy-Name", "Tracked-Energy"]
        self.history = [[], []]

    def update(self, sm, step):
        if self._background:
            log.debug("EnergyTracking calling eval_energy with background")
            energy = self._energy_function.eval_energy(sm.bg, background=True, sampled_stats=sm.elem_defs)
            self._energy_function.accept_last_measure()
        else:
            log.debug("EnergyTracking calling eval_energy without background")
            energy = self._energy_function.eval_energy(sm.bg, background=False, sampled_stats=sm.elem_defs)
        self.history[0].append(self._energy_function.shortname)
        self.history[1].append(energy)
        if isinstance(energy, np.ndarray) and len(energy) == 1:
            energy = "{:10.3f}".format(energy[0])
        elif isinstance(energy, float):
            energy = "{:10.3f}".format(energy)
        sn = self._energy_function.shortname
        if sn == "":
            sn = self._energy_function.get_name()
        return "{}: {}".format(sn, energy)

    @property
    def header_str(self):
        return "Tracked Energy"

    @staticmethod
    def parse_value(stri):
        stri = stri.split()
        return stri[0][:-1], float(stri[1])


class EnergyMeasure(StatisticsCollector):
    """
    After every step, log the last accepted measure of the energy. Does not call eval_energy!
    """
    header = ["measure_of_"]

    def __init__(self, energy_function):
        """
        :param energy_function: The energy function from which the last measure will be logged.
        """
        super(EnergyMeasure, self).__init__()
        self._energy_function = energy_function
        if isinstance(energy_function, fbe.AMinorEnergy):
            header = []
            for loop in self._energy_function.LOOPS:
                header.append("measure_of_{}({})".format(
                    self._energy_function.shortname, loop))
            self.header = header
            self.history = [[] for _ in range(len(header))]
        else:
            self.header = ["measure_of_"+self._energy_function.shortname]
            self.history = [[]]

    def update(self, sm, step):
        measure = self._energy_function.last_accepted_measure
        if isinstance(self._energy_function, fbe.AMinorEnergy):
            out = []
            for i, loop in enumerate(self._energy_function.LOOPS):
                self.history[i].append(measure[loop])
                out.append("{} {:d}".format(loop, measure[loop]))
            return "\t".join(out)
        else:
            self.history[0].append(measure)
            return "{:10.3f}".format(measure)

    @staticmethod
    def parse_value(stri):
        if stri[0] in "ih":
            print(stri)
            return int(stri.split(":")[1])
        return float(stri)


class ShowTime(StatisticsCollector):
    """
    After every step, show the elapsed time (since start_time)
    """
    header = ["time"]

    def __init__(self, start_time):
        """
        :param start_time: A numeric value or the string "now".
                           Elapsed time since start_time will be printed, collected.
        """
        super(ShowTime, self).__init__()
        if start_time == "now":
            self._start_time = time.time()
        else:
            self._start_time = start_time
        self.history = None  # This does not save its history.

    def update(self, sm, step):
        elapsed = time.time()-self._start_time
        if elapsed < 1000:
            return "{:5.1f} sec".format(elapsed)
        elif elapsed/60.0 < 1000:
            return "{:5.1f} min".format(elapsed/60.0)
        else:
            return "{:5.1f} h".format(elapsed/3600.0)

    @staticmethod
    def parse_value(stri):
        time, unit = stri.split()
        time = float(time)
        if unit == "h":
            time *= 3600
        elif unit == "min":
            time *= 60
        elif unit == "sec":
            pass
        else:
            assert False
        return time


class Delimitor(StatisticsCollector):
    """
    A dummy collector that always prints a delimitor and does not save any Statistics
    """
    header = []  # Empty header, because no need to parse the dilimiter in file-parsing mode.

    def __init__(self, delimitor="|"):
        super(Delimitor, self).__init__()
        self._delimitor = delimitor
        self.history = None
        self.header = ["|"]

    def update(self, sm, step):
        return self._delimitor


class Distance(StatisticsCollector):
    """
    The distance between two nucleotides
    """
    header = ["Distance_"]

    def __init__(self, nuc1, nuc2):
        super(Distance, self).__init__()
        self._nuc1 = int(nuc1)
        self._nuc2 = int(nuc2)
        self.header = ["Distance_{}-{}".format(self._nuc1, self._nuc2)]

    def update(self, sm, step):
        dist = ftuv.vec_distance(sm.bg.get_virtual_residue(self._nuc1, True),
                                 sm.bg.get_virtual_residue(self._nuc2, True))
        self.history[0].append(dist)
        return "{:6.2f} A".format(dist)

    @staticmethod
    def parse_value(stri):
        return float(stri.split()[0])


class AngleTracker(StatisticsCollector):
    """
    The distance between two nucleotides
    """
    header = ["Angle_"]

    def __init__(self, stem1, stem2):
        super(AngleTracker, self).__init__()
        self._elem1 = stem1
        self._elem2 = stem2
        self.header = ["Angle_{}-{}".format(self._elem1, self._elem2)]

    def update(self, sm, step):
        angle = ftuv.vec_angle(sm.bg.coords.get_direction(self._elem1),
                               sm.bg.coords.get_direction(self._elem2))
        self.history[0].append(angle)
        return "{:6.2f}".format(math.degrees(angle))

    @staticmethod
    def parse_value(stri):
        return math.radians(float(stri))

class LocalFragmentCoverage(StatisticsCollector):
    """
    What fraction of different stats (fragments) were accepted.

    Related: 10.1002/prot.24987
    """
    header_str = "Local-Coverage"
    header = [header_str]

    def __init__(self, stat_source, cg):
        super(LocalFragmentCoverage, self).__init__()
        self.history = None
        if stat_source is None:
            self.silent = True
            self.history = None
            return
        #: A fess.builder.stat_container.StatStorage instance
        self._stat_source = stat_source
        #: Keep track of the coverage for all elements, so it only has to be updated if elements changed.
        self._coverage =  defaultdict(lambda: 0.)
        #: For each element, store the pdb_names of stats that have been accepted at least once
        self._known_fragments = defaultdict(set)
        #: The elem_defs as {elem : pdb_name} from the previouse iteration, to avoid recalculation.
        self._prev_elem_defs = defaultdict(lambda: None)

    def update(self, sm, step):
        if not self.silent:
            if self.history is None:
                self.history = [[] for x in range(len(sm.bg.defines))]
            elem_defs = {elem: stat.pdb_name for elem,
                         stat in sm.elem_defs.items()}
            for elem, stat_name in elem_defs.items():
                if stat_name != self._prev_elem_defs[elem]:
                    # We have a new stat for this elem. Update coverage
                    self._known_fragments[elem].add(stat_name)
                    new_cov = self._stat_source.coverage_for(
                        self._known_fragments[elem], sm.bg, elem)
                    assert new_cov >= self._coverage[elem]
                    self._coverage[elem] = new_cov

            coverage = sum(self._coverage.values())/len(self._coverage)
            for i, elem in enumerate(sorted(self._coverage)):
                self.history[i].append(self._coverage[elem])
            return "{:.1%}".format(coverage)

    @staticmethod
    def parseValue(stri):
        assert stri[-1] == "%"
        return float(stri[:-1])/100


class MlClosingOptions(StatisticsCollector):
    """
    How many stats from the stat container would fit the broken multiloops?
    """
    header_str = "Fitting_ml_stats"
    header = [header_str]

    def __init__(self, stat_source, cutoff=3):
        super(MlClosingOptions, self).__init__()
        self._stat_source = stat_source
        self._cutoff = cutoff
        self.complete_multiloops = []

    def update(self, sm, step):
        if not self.silent:
            self.history[0].append([])
            count = 0
            for ml in sm.bg.mloop_iterator():
                if ml not in sm.bg.mst:
                    count += self._find_stats_for(ml, sm)
            return '{:3d}'.format(count)

    def _find_stats_for(self, ml, sm):
        count = 0
        for stat in self._stat_source.iterate_stats_for(sm.bg, ml):
            diff = sm.ml_stat_deviation(ml, stat)
            if diff < self._cutoff:
                self.history[0][-1].append('{}:{}'.format(stat.pdb_name, diff))
                all_mls = {}
                for m in sm.bg.mloop_iterator():
                    if m in sm.bg.mst:
                        all_mls[m] = sm.elem_defs[m].pdb_name
                    elif m == ml:
                        all_mls[m] = stat.pdb_name
                    else:
                        all_mls[m] = None
                self.complete_multiloops.append(all_mls)
                count += 1
        return count

    @staticmethod
    def parseValue(stri):
        return int(stri)

    def to_file(self):
        df = pd.DataFrame(self.complete_multiloops)
        df.to_csv(os.path.join(
            conf.Configuration.sampling_output_dir, "ml_closing_options.txt"))


###################################################################################################
##########
###################################################################################################
_statisticsDefaultOptions = {
    "showtime": "now",
    "rmsd": True,
    "drmsd": False,
    "acc": True,
    "rog": True,
    "name": True,
    "length": True,
    "step": True,
    "local_coverage": True,
    "extreme_energies": "lh",
    "reference_energy": None,
    "extreme_rmsd": True,
    "silent": False,
    "constituing_energies": True,
    "history": "all",
    "step_save": 0,
    "save_n_best": 0,
    "save_min_rmsd": 0,
    "measure": [],
    "distance": [],
    "angles":False,
    "asphericity": True,
    "anisotropy": True,
    "ml_closing": False
}

def remove_common_pre_and_postfixes(fns):
    for i in range(len(fns[0])):
        try:
            if not all(fn[i]==fns[0][i] for fn in fns):
                break
        except IndexError:
            break
    else:
        return [""]*len(fns)
    fns=[fn[i:] for fn in fns]
    for i in range(1, len(fns[0])):
        i=-i
        try:
            if not all(fn[i]==fns[0][i] for fn in fns):
                break
        except IndexError:
            break
    else:
        assert False
    fns=[fn[:-i+1] for fn in fns]
    return fns


class SamplingStatistics:
    def __init__(self, cg, energy_functions=[], stat_source=None,
                 options="all", output_directory=None, args=None):
        """
        :param cg: The CoarseGrainRNA against which to collect statistics.
                        .. warning:: This should not be modified.
                                     (Its best to create a new SpatialModel only for this
                                     constructor with a deepcopy of the BulgeGraph used elsewhere.)
        :param enery_functions: A list of fbe.CombinedEnergy or fbe.EnergyFunction instances.
                      These energies will be evaluated with "background=False" every time
                      `update_statistics` is called.
                      If the energy used for sampling uses a background KDE, it makes sense to
                      pass a reference to this energy here as well.
        :param options: What statistics to calculate and print, what files to save.
                      A string like "all" or a dictionary, which will be used to update
                      the default option dictionary.
                      Key-Value pairs could be:

                      * `"rog": True|False`:
                            Radius of Gyration
                      * `"mcc": True|False`:
                            Matthews correlation coefficient to sm_orig
                            If sm_orig is missing 3D coordinates, this will always be false.
                      * `"rmsd": True|False`:
                            Root Mean Square Deviation to sm_orig
                            If sm_orig is missing 3D coordinates,
                            this will always be false.
                      * `"extreme_rmsd": TRUE | FALSE | "max"`
                            Whether to show the minimal and maximal RMSD.

                      * `"name": True|False`:
                            The name of the spatial model.
                      * `"length": True|False`
                            The length of the RNA
                      * `"step": True|False`
                            Show the sampling step
                      * `"extreme_energies": STRING`
                            "lh": show lowest_energy and highest_energy
                            "l": lowest energy only, "h": highest only
                            The order is not relevant, unrecognized characters are ignored.

                      * `"constituing_energies": TRUE|FALSE|"no_clash"`
                            Show the energies that generated the sampling energies
                            (If the sampling energy is a CombinedEnergy)
                            "no_clash": Show constituing energies except clash energies.
                      * `"showtime": starttime|"now"|False`:
                            Show the elapsed time at each iteration step.
                            If `starttime` is given, show elapsed time since start time, else
                            since the creation time of this object.
                      * `"step_save": INT`:
                            Save the CoarseGraiNRNA every n steps. If this is 0, don't save it.
                      * `"save_n_best": INT`:
                            Save the n structures with lowest energy
                      * `"save_min_rmsd": INT`:
                            Save the n structures with lowest RMSD.
                      * `"measure": [energy_function, ...]
                            For the energy functions in the list, print the measure (not the enrgy).
                            This is useful for energies with a background.
                      * `"distance": A list of tuples of ints.
                            Display the distance between these two nucleotides.
        """
        self.collector = None
        self.step = 0
        self._outfile = None

        collectors = []
        if options == "all":
            self.options = _statisticsDefaultOptions
        elif isinstance(options, dict):
            self.options = copy.copy(_statisticsDefaultOptions)
            self.options.update(options)
        else:
            raise ValueError("Options '{}' not recognized. Please use a  "
                             "dict or the string 'all'".format(options))

        if self.options["rog"]:
            collectors.append(ROGStatistics())
        if self.options["acc"]:
            collectors.append(ACCStatistics(cg))
        if self.options["rmsd"]:
            r_col = RMSDStatistics(cg, show_min_max=self.options["extreme_rmsd"],
                                   save_n_best=self.options["save_min_rmsd"])
            if not r_col.silent:
                collectors.append(Delimitor())
                collectors.append(r_col)

        if self.options["drmsd"]:
            dr_col = RMSDStatistics(cg, show_min_max=False,
                                    mode="dRMSD")
            if not dr_col.silent:
                collectors.append(Delimitor())
                collectors.append(dr_col)

        if self.options["stem_rmsds"]:
            names = remove_common_pre_and_postfixes(self.options["stem_rmsds"])
            for i, fn in enumerate(self.options["stem_rmsds"]):
                collectors.append(StemRMSD(fn, names[i], args))
            collectors.append(Delimitor())




        if self.options["asphericity"]:
            collectors.append(Delimitor())
            collectors.append(AsphericityStatistics())
        if self.options["anisotropy"]:
            collectors.append(AnisotropyStatistics())


        if self.options["distance"]:
            collectors.append(Delimitor())
            for n1, n2 in self.options["distance"]:
                collectors.append(Distance(n1, n2))

        collectors.append(Delimitor())

        if self.options["angles"]:
            for elem in cg.defines:
                if elem[0] in "im":
                    stem1, stem2 = cg.connections(elem)
                    collectors.append(AngleTracker(stem1, stem2))
            collectors.append(Delimitor())


        if self.options["local_coverage"]:
            collectors.append(LocalFragmentCoverage(stat_source, cg))
            collectors.append(Delimitor())

        for ef in energy_functions:
            log.debug("Add Energy tracking for %s", ef.shortname)
            collectors.append(EnergyTracking(ef))
        if energy_functions:
            collectors.append(Delimitor())

        for m in self.options["measure"]:  # This has to be AFTER tracking energies!
            collectors.append(EnergyMeasure(m))

        if self.options["ml_closing"] is not False:
            collectors.append(Delimitor())
            collectors.append(MlClosingOptions(stat_source))

        if self.options["showtime"] is not False:
            collectors.append(ShowTime(self.options["showtime"]))

        # Note: SortedCollection(maxlen=0) is valid!
        # should contain tuples: `(bg, energy)`
        self.best_cgs = SortedCollection(
            key=lambda x: x[1], maxlen=self.options["save_n_best"])

        self.collector = CombinedStatistics(collectors)

        if output_directory is None:
            self.out_dir = conf.Configuration.sampling_output_dir
        else:
            self.out_dir = output_directory
            if not os.path.exists(self.out_dir):
                os.makedirs(self.out_dir)
                log.info("Directory {} created.".format(self.out_dir))

        if self.options["reconstruct_n"] > 0 and self.options["reconstruct_cg"] and self.options["reconstruct_pdb"]:
            self.reconstructor = fbr.Reconstructor(self.options["reconstruct_pdb"],
                                                   self.options["reconstruct_cg"],
                                                   False,
                                                   self.options["reconstruct_cache_dir"])
        else:
            self.reconstructor = None

    @contextlib.contextmanager
    def open_outfile(self):
        try:
            with open(os.path.join(self.out_dir, "out.log"), "w") as f:
                self._outfile = f
                self.print_header()
                log.debug("Opened file %s", os.path.join(
                    self.out_dir, "out.log"))
                yield f
        finally:
            self._outfile = None

    def printline(self, line):
        """Print to both STDOUT and the log file."""
        print(line)
        if self._outfile is None:
            log.warning("No output-file opened!")
        else:
            print(line, file=self._outfile)
            self._outfile.flush()

    def print_header(self):
        """
        Print the header line
        """
        contribs = ["Step\tSampling_Energy"]
        if self.options["constituing_energies"]:
            contribs.append("( Constituing_Energies )")
        contribs.append(self.collector.header_str)
        contribs.append("Sampling Move")
        contribs.append("Rej.Clashes")
        contribs.append("Rej.BadMls")

        self.printline("\t".join(contribs))

    def update_statistics(self, sm, energy, member_energies=[], change="", clashes=[], bad_mls=[]):
        """
        Add a new structure to the statistics.

        :param sm: The spatial model.
        :param energy: The energy used for sampling of the structure, or None.
        :param member_energies: A list of tuples `(energy_shortname, value)`
        """
        self.step += 1
        line = ["{:6d}\t{:10.3f}".format(self.step, energy)]
        if self.options["constituing_energies"] == "no_clash":
            ignore_names = [fbe.RoughJunctionClosureEnergy().shortname,
                            fbe.StemVirtualResClashEnergy().shortname]
        else:
            ignore_names = []
        line.append("( "+" ".join("{} {:10.3f}".format(*x) for x in member_energies
                                  if x[0] not in ignore_names)+" )")
        line.append(self.collector.update(sm, self.step))
        line.append(change)
        if isinstance(clashes, str):
            line.append(clashes)
        else:
            line.append(",".join(["{}+{}".format(*c) for c in clashes]) or "noclash" )
        if isinstance(bad_mls, str):
            line.append(bad_mls)
        else:
            line.append(",".join(bad_mls) or "noclash")

        self.printline("\t".join(line))

        if self.best_cgs.can_insert_right((None, energy)):
            self.best_cgs.insert_right((sm.bg.to_cg_string(), energy))

        if self.step % 10 == 0:
            for i, cg_stri in enumerate(self.best_cgs):
                with open(os.path.join(self.out_dir,
                                       'best{:d}.coord'.format(i)), 'w') as f:
                    f.write(cg_stri[0])

        if self.options["step_save"] > 0 and self.step % self.options["step_save"] == 0:
            cg_stri = sm.bg.to_cg_string()
            with open(os.path.join(self.out_dir,
                                   'step{:06d}.coord'.format(self.step)), "w") as f:
                # print("Opened file", os.path.join(self.out_dir,
                #                  'step{:06d}.coord'.format(self.step)))
                print(cg_stri, file=f)

        if self.reconstructor is not None and self.step % self.options["reconstruct_n"] == 0:
            try:
                chains = self.reconstructor.reconstruct(sm)
                ftup.output_multiple_chains(list(chains.values()),
                                            os.path.join(self.out_dir,
                                                         'step{:06d}.reconstr.pdb'.format(self.step)))
            except:
                log.exception(
                    "Could not reconstruct all-atom PDB, because of the following error:")


####################################################################################
# Handling of commandline arguments
####################################################################################
def update_parser(parser):
    monitor_options = parser.add_argument_group("Controlling output",
                                                description="These options control what statistics "
                                                "about the sampling trajectory ERNWIN records and what output files ERNWIN writes.")
    monitor_options.add_argument('--save-n-best', default=3,
                                 help='Save the best (lowest energy) n structures.', type=int)
    monitor_options.add_argument('--save-min-rmsd', default=3,
                                 help='Save the best (lowest rmsd) n structures.', type=int)
    monitor_options.add_argument('--step-save', default=100, help="Save the structure at every n'th step.",
                                 type=int)
    monitor_options.add_argument('--reconstruct-every-n', default=500,
                                 help="Reconstruct every n'th structure to an all-atom pdb.\n"
                                 "Use 0 to disable reconstruction",
                                 type=int)

    monitor_options.add_argument("--source-pdb-dir", type=str,
                                 help="A directory where all the pdb files, from which fragments "
                                 "can be extracted, are located.")
    monitor_options.add_argument("--source-cg-dir", type=str,
                                 help="A directory where the cg-files corresponding to the pdb files in "
                                 "source-pdb-dir are located.")
    monitor_options.add_argument("--reconstruction-cache-dir", type=str,
                                 help="A directory where ernwin can store fragments \n"
                                      "used for reconstruction.\n"
                                      "If this is given, ernwin stores fragments it extracts from PDB files\n"
                                      "in this directory, and uses them from there when it needs the same \n"
                                      "fragment again. This speeds up reconstruction compared to always loading\n"
                                      "the full PDB file.")

    monitor_options.add_argument('--dump-energies', default=False, action='store_true',
                                 help='Dump the measures used for energy calculation to file')  # UNUSED OPTION. REMOVE
    monitor_options.add_argument('--no-rmsd', default=False,
                                 help='Refrain from trying to calculate the rmsd.', action='store_true')
    monitor_options.add_argument('--dist', type=str,
                                 help="One or more pairs of nucleotide positions.\n"
                                 "The distance between these nucleotifdes will be \n"
                                 "calculated.\n"
                                 "Example: '1,15:3,20..' will print the distance between\n"
                                 "          nucleoitide 1 and 15 and the distance \n"
                                 "          between nucleotide 3 and 20.")
    monitor_options.add_argument('--angles', action="store_true",
                                 help="Track angles between all adjacent stems")
    monitor_options.add_argument('--track-energies', action='store', type=str, default="",
                                 help="A list of energies.\n"
                                 "Each energy is given in the format specified\n"
                                 "for the --energy option.\n"
                                 "These energies are not used for sampling, \n"
                                 "but only calculated after each step.")
    monitor_options.add_argument('--stem-rmsd-to', type=str,
                                 help="A comma-seperated list of rna filenames.\n"
                                      "track the stem-rmsd to these structures.")
    monitor_options.add_argument('--fewer-statistics', default=False,
                                 help='Do not output so many statistics every step.', action='store_true')


def from_args(args, cg, energy, stat_source, out_dir, show_min_rmsd=True):
    options = {}
    if args.fewer_statistics:
        options.update({
            "acc": False,
            "name": False,
            "length": False,
            "local_coverage": False,
            "extreme_rmsd": False,
            "asphericity": False,
            "anisotropy": False,
            "ml_closing": False})
    if args.no_rmsd:
        options["rmsd"] = False
        options["acc"] = False
    options["step_save"] = args.step_save
    options["reconstruct_n"] = args.reconstruct_every_n
    options["reconstruct_pdb"] = args.source_pdb_dir
    options["reconstruct_cg"] = args.source_cg_dir
    options["reconstruct_cache_dir"] = args.reconstruction_cache_dir
    options["save_n_best"] = args.save_n_best
    options["save_min_rmsd"] = args.save_min_rmsd
    options["measure"] = []
    if not show_min_rmsd:
        options["extreme_rmsd"] = "max"
    if args.stem_rmsd_to:
        fns = args.stem_rmsd_to.split(",")
        options["stem_rmsds"]=fns
    else:
        options["stem_rmsds"]=[]

    if args.dump_energies:
        for e in energy.energies:
            options["measure"].append(e)
    if args.dist:
        options["distance"] = list(map(str.split, args.dist.split(
            ':'), it.repeat(",")))  # map is from future!
    else:
        options["distance"] = []

    if args.angles:
        options["angles"]=True
    else:
        options["angles"]=False
    track_energies = fbe.EnergyFunction.from_string(
        args.track_energies, cg=cg, iterations=None, stat_source=stat_source,
        pdd_target=args.pdd_file)
    assert isinstance(energy, fbe.CombinedEnergy)
    track_energies += [e for e in energy.energies
                       if isinstance(e, (fbe.CombinedEnergy, fbe.CoarseGrainEnergy, fbe.InteractionEnergy))]

    return SamplingStatistics(copy.deepcopy(cg), track_energies, stat_source,
                              output_directory=out_dir, options=options, args=args)
####################################################################################
# Parsing files generated by this module
####################################################################################


class OutfileParser(object):
    def __init__(self):
        self._collectors = None

    def _get_lookup_table(self):
        """
        Return a dictionary with headers as keys and
        StatisticCollector classes (not instances) as values
        """
        return {h: cls for cls in StatisticsCollector.__subclasses__()
                if not cls == CombinedStatistics
                for h in cls.header}

    def _collector_from_header(self, header):
        """
        Return the StatisticsCollector subclass which supports the given header.

        Note that the header ((in the parsed file) might contain additional
        characters to the right that are not present in the classes header list.
        Thus, if key-based lookup fails, this method falls back to startswith()
        """
        try:
            return self._lookup_table[header]
        except KeyError:
            for key, cls in self._lookup_table.items():
                if header.startswith(key):
                    return cls
        if header=="Sampling Move":
            return header
        log.error("Cannot Parse header '%s'", header)
        return None

    def _init_collector_lookup(self, headers):
        self._lookup_table = self._get_lookup_table()
        self._collectors = []
        for header in headers:
            self._collectors.append(self._collector_from_header(header))

    @staticmethod
    def parse(filepath):
        parser = OutfileParser()
        return parser._parse(filepath)

    def _parse(self, filepath):
        meta = {}
        with open(filepath) as file:
            headers = None
            data = None
            for line_no, line in enumerate(file):
                try:
                    line = line.strip()
                    if not line:
                        continue
                    elif line.startswith("# Random Seed:"):
                        meta["seed"] = int(line.split()[-1])
                    elif line.startswith("# Command"):
                        meta["command"] = line.split('`')[1]
                    elif line.startswith("# Version"):
                        fields = line.split()
                        meta["ernwin_version"] = fields[3].rstrip(",")
                        meta["forgi_version"] = fields[5]
                    elif line.startswith("#"):
                        continue
                    elif headers is None:
                        headers = line.split("\t")
                        self._init_collector_lookup(headers)
                        data = []
                        for i in range(len(headers)):
                            data.append([])
                    else:
                        fields = line.split('\t')
                        for i, field in enumerate(fields):
                            if i == 0:  # Step
                                data[i].append(int(field))
                            elif i == 1:  # Sampling_Energy
                                data[i].append(float(field))

                            cls = self._collectors[i]
                            if cls=="Sampling Move":
                                data[i].append(field)
                            elif cls is not None:
                                data[i].append(cls.parse_value(field))


                except Exception as e:
                    with log_to_exception(log, e):
                        log.error(
                            "Exception occurred during parsing of line %d '%s'", line_no, line)
                    raise
            data_dic = {}
            for i, header in enumerate(headers):
                if data[i]:
                    if isinstance(data[i][0], tuple):
                        data_dic["{}_{}".format(header, data[i][0][0])] = [
                            x[1] for x in data[i]]
                    else:
                        data_dic[header] = data[i]
            data_dic["move_type"]=[]
            data_dic["accepted"]=[]
            data_dic["delta_E"]=[]
            data_dic["stats_moved"]=[]

            for d in data[-1]:
                field, _, accepted = d.rpartition(";")
                typ, _, field=field.partition(":")
                data_dic["accepted"].append(accepted)
                if typ=="RE":
                    data_dic["delta_E"].append(float("nan"))
                    data_dic["move_type"].append("RE")
                    data_dic["stats_moved"].append(float("nan"))
                else:
                    self.update_data_move(data_dic, typ, field)
        return data_dic


    def update_data_move(self, data, loop, field):
        fields=field.split(";")
        e0,_,e1 = fields[-1].partition("->")
        data["delta_E"].append(float(e1)-float(e0))
        data["stats_moved"].append(len(fields)-1)
        data["move_type"].append(loop[0])
