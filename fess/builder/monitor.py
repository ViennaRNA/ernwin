from __future__ import absolute_import, division, print_function, unicode_literals
from builtins import (ascii, bytes, chr, dict, filter, hex, input,
                      int, map, next, oct, open, pow, range, round,
                      str, super, zip)
from future.builtins.disabled import (apply, cmp, coerce, execfile,
                             file, long, raw_input, reduce, reload,
                             unicode, xrange, StandardError)

import forgi.threedee.model.descriptors as ftur
import forgi.threedee.model.similarity as ftme
import forgi.threedee.utilities.vector as ftuv
from . import config as conf
from . import energy as fbe
import sys
import time
import copy
import os.path
import numpy as np
from ..aux.SortedCollection import SortedCollection
import warnings
import logging
from collections import defaultdict
log = logging.getLogger(__name__)
__metaclass__=type

class StatisticsCollector(object):
    """
    A base class for Spatial-Model Statistic objects.
    """       
    header=[] #This is used for outfile-parsing.
    def __init__(self, *args, **kwargs):
        """
        :var self.header: A list of strings to describe the fields of history collected.
        :var self.history: A list of lists.
        :var self.silent: This is set to true, if a Collector cannot initialize itself properly.
        :var self.display: If this is True, the collector will write to stdout.
        """
        self.history=[[]]
        self.silent = False
        self.display=True
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
        self._members=collectors
        if separator:
          self._joiner="\t"+separator+"\t"
        else:
          self._joiner="\t"

    @property
    def header_str(self):
        return self._joiner.join(member.header_str for member in self._members if not member.silent)
    @property
    def header(self):
        """
        This only includes headers for which a history is written.
        """
        return [ header for member in self._members for header in member.header_list 
                 if member.history is not None and not member.silent ]
    @property
    def history(self):
        hist = []
        for member in self._members:
            if member.history is None or member.silent:
                continue
            hist+=member.history
        return hist
        #return [ history for member in self._members for history in member.history if not member.silent and member.history is not None ]
    def update(self, sm, step):
        line=[]
        for member in self._members:
            if not member.silent:
                log.debug("Updating {}".format(type(member)))
                fields = member.update(sm, step)
                if member.display:
                    line.append(fields)
        return self._joiner.join(line)

    def to_file(self):
        history=self.history
        length = len(history[0])
        with open(os.path.join(conf.Configuration.sampling_output_dir,"monitor.txt"), "w") as f:
            print("#"+self._joiner.join(h for member in self._members for h in member.header if not member.silent and member.history is not None ), file=f)
            for i in range(length):
                try:
                    line = []
                    for j, h in enumerate(history):
                        line.append(str(h[i]))
                except Exception as e:
                    log.error("Members are: {}".format(", ".join(map(str, self._members))))
                    log.exception("Could not dump complete history to file 'monitor.txt': {} occurred for history column {} line {}".format(type(e).__name__, j, i))
                    return
                else:
                    print( self._joiner.join(line), file=f)

        
                
            


class ROGStatistics(StatisticsCollector):
    """
    Store and print the Radius of Gyration.
    """        
    header=["ROG"]
    def update(self, sm, step):
        rog=ftur.radius_of_gyration(sm.bg.get_ordered_stem_poss())
        self.history[0].append(rog)
        return "{:6.2f} A".format(rog)
    @staticmethod
    def parse_value(stri):
        stri = stri.split()
        assert stri[1]=="A", stri
        return float(stri[0])

class AsphericityStatistics(StatisticsCollector):
    """
    Store and print the Asphericity.
    """        
    header=["Asphericity"]
    def update(self, sm, step):
        rog=ftur.asphericity(sm.bg.get_ordered_stem_poss())
        self.history[0].append(rog)
        return "{:6.2f}".format(rog)
    @staticmethod
    def parse_value(stri):
        return float(stri)
        
class AnisotropyStatistics(StatisticsCollector):
    """
    Store and print the Anisotropy.
    """        
    header=["Anisotropy"]
    def update(self, sm, step):
        rog=ftur.anisotropy(sm.bg.get_ordered_stem_poss())
        self.history[0].append(rog)
        return "{:6.2f}".format(rog)
    @staticmethod
    def parse_value(stri):
        return float(stri)
        
        
class ACCStatistics(StatisticsCollector):
    """
    Store and print the Adjacency Correlation Coefficient
    """            
    header=[ "ACC" ]
    def __init__(self, reference_sm):
        super(ACCStatistics, self).__init__()
        try:
            self._cm_calc = ftme.AdjacencyCorrelation(reference_sm.bg)
        except Exception as e:
            log.exception(e)
            warnings.warn("Cannot report ACC. {} in reference SM: {}".format(type(e), str(e)))
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

class RMSDStatistics(StatisticsCollector):
    """
    Store and print the Root Mean Square Deviation from the reference SM.
    """
    header = ["RMSD", "dRMSD", "maxRMSD", "maxdRMSD" ] #All allowed headers. __init__ sets self.header, which takes priority.
    def __init__(self, reference_sm, show_min_max=True, save_n_best=0, mode="RMSD"):
        """
        :param reference_sm: The reference spatial model, against which to collect statistics.
        :param show_min_max: Bool. Whether to show the minimal and maximal RMSD so far.
        :param mode: "RMSD" or "dRMSD"
        """
        super(RMSDStatistics, self).__init__()
        self.best_cgs=SortedCollection(key=lambda x: x[1], maxlen=save_n_best)
        self.mode=mode
        try:
            self._reference = reference_sm.bg.get_ordered_virtual_residue_poss()
        except:
            self.silent=True
            self.history=None
        else:
            self._showMinMax=show_min_max
            if show_min_max=="max":
                self.header = [ mode, "max"+mode ]
                self.history= [ [],[] ]
                self._maxRMSD = float("-inf")
            elif show_min_max:
                self.header = [ mode, "min"+mode, "max"+mode ]
                self.history= [ [],[],[] ]
                self._minRMSD = float("inf")
                self._maxRMSD = float("-inf")
            else:
                self.header = [mode]

    def update(self, sm, step):
        if not self.silent:
            curr_vress=sm.bg.get_ordered_virtual_residue_poss()
            if self.mode=="RMSD":
              rmsd=ftme.rmsd(self._reference, curr_vress)
            elif self.mode=="dRMSD":
              rmsd=ftme.drmsd(self._reference, curr_vress)
            if self.best_cgs.can_insert_right((None, rmsd)):
                self.best_cgs.insert_right((sm.bg.to_cg_string(),rmsd))
            if step % 10 == 0:
                for i, bg in enumerate(self.best_cgs):
                    with open(os.path.join(conf.Configuration.sampling_output_dir, 
                              'best_rmsd{:d}.coord'.format(i)), 'w') as f:
                        f.write(bg[0])
            if self._showMinMax:
                if rmsd>self._maxRMSD:
                    self._maxRMSD=rmsd
                if self._showMinMax==True:
                    if rmsd<self._minRMSD:
                        self._minRMSD=rmsd                
                    self.history[1].append( self._minRMSD )
                self.history[0].append( rmsd )
                self.history[-1].append( self._maxRMSD )
                if self._showMinMax=="max":
                      return "{:6.3f} A\t{:6.3f} A".format(rmsd, self._maxRMSD)
                else:
                    return "{:6.3f} A\t{:6.3f} A\t{:6.3f} A".format(rmsd, self._minRMSD, self._maxRMSD)
            else:
                self.history[0].append( rmsd )
                return "{:6.3f} A".format(rmsd)
        else: 
            return
    @staticmethod
    def parse_value(stri):
        stri = stri.split()
        assert stri[1]=="A"
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
        self.header = [ "Energy-Name", "Tracked-Energy" ]
        self.history = [ [], [] ]
    def update(self, sm, step):
        if self._background:
            energy=self._energy_function.eval_energy(sm.bg, background=True)
            self._energy_function.accept_last_measure()
        else:
            energy=self._energy_function.eval_energy(sm.bg, background=False)
        self.history[0].append(self._energy_function.shortname())
        self.history[1].append(energy)
        if isinstance(energy, np.ndarray) and len(energy)==1:
            energy="{:10.3f}".format(energy[0])
        elif isinstance(energy, float):
            energy="{:10.3f}".format(energy)
        sn = self._energy_function.shortname
        if sn=="":
            sn = self._energy_function.get_name()
        return "{}: {}".format(sn, energy)
    @property
    def header_str(self):
        return "Tracked Energy"
    
    @staticmethod
    def parse_value(stri):
        stri = stri.split()
        return stri[0][:-1],float(stri[1])

class EnergyMeasure(StatisticsCollector):
    """
    After every step, log the last accepted measure of the energy. Does not call eval_energy!
    """
    header = [ "measure_of_" ]
    def __init__(self, energy_function):
        """
        :param energy_function: The energy function from which the last measure will be logged.
        """
        super(EnergyMeasure, self).__init__()
        self._energy_function = energy_function
        self.header = [ "measure_of_"+self._energy_function.shortname ]
        self.history = [ [] ]
    def update(self, sm, step):
        measure=self._energy_function.accepted_measures[-1]
        self.history[0].append(measure)
        return "{:10.3f}".format(measure)
    @staticmethod
    def parse_value(stri):
        return float(stri)

class ShowTime(StatisticsCollector):
    """
    After every step, show the elapsed time (since start_time)
    """
    header = [ "time" ]
    def __init__(self, start_time):
        """
        :param start_time: A numeric value or the string "now". 
                           Elapsed time since start_time will be printed, collected.
        """
        super(ShowTime, self).__init__()
        if start_time == "now":
            self._start_time=time.time()
        else:
            self._start_time=start_time
        self.history= None #This does not save its history.
    def update(self, sm, step):
        elapsed=time.time()-self._start_time
        if elapsed<1000:
            return "{:5.1f} sec".format(elapsed)
        elif elapsed/60.0<1000:
            return "{:5.1f} min".format(elapsed/60.0)
        else:
            return "{:5.1f} h".format(elapsed/3600.0)
    @staticmethod
    def parse_value(stri):
        time, unit = stri.split()
        time=float(time)
        if unit == "h":
            time*=3600
        elif unit == "min":
            time*=60
        elif unit == "sec":
            pass
        else:
            assert False
        return time

class Delimitor(StatisticsCollector):
    """
    A dummy collector that always prints a delimitor and does not save any Statistics
    """        
    header = [] #Empty header, because no need to parse the dilimiter in file-parsing mode.
    def __init__(self, delimitor="|"):
        super(Delimitor, self).__init__()
        self._delimitor=delimitor
        self.history = None
        self.header= ["|"]
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
        self.header=["Distance_{}-{}".format(self._nuc1, self._nuc2)]
    def update(self, sm, step):
        dist = ftuv.vec_distance(sm.bg.get_virtual_residue(self._nuc1, True), 
                                 sm.bg.get_virtual_residue(self._nuc2, True))
        self.history[0].append(dist)
        return "{:6.2f} A".format(dist)
    @staticmethod
    def parse_value(stri):
        return float(stri.split()[0])


class LocalFragmentCoverage(StatisticsCollector):
    """
    What fraction of different stats (fragments) were accepted.
    
    Related: 10.1002/prot.24987
    """
    header_str = "Local-Coverage"
    header = [ header_str ]
    def __init__(self, stat_source, cg):
        super(LocalFragmentCoverage, self).__init__()
        self.history = [ [] for x in range(len(cg.defines)) ]
        if stat_source is None:
            self.silent=True
            self.history=None
            return
        #: A fess.builder.stat_container.StatStorage instance
        self._stat_source = stat_source
        #: Keep track of the coverage for all elements, so it only has to be updated if elements changed.
        self._coverage = {elem : 0. for elem in cg.defines}
        #: For each element, store the pdb_names of stats that have been accepted at least once
        self._known_fragments = defaultdict(set)
        #: The elem_defs as {elem : pdb_name} from the previouse iteration, to avoid recalculation.  
        self._prev_elem_defs = {elem : None for elem in cg.defines}
            
    def update(self, sm, step):
        if not self.silent:
            elem_defs = {elem: stat.pdb_name for elem, stat in sm.elem_defs.items()}
            for elem, stat_name in elem_defs.items():
                if stat_name != self._prev_elem_defs[elem]:
                    # We have a new stat for this elem. Update coverage
                    self._known_fragments[elem].add(stat_name)
                    new_cov = self._stat_source.coverage_for(self._known_fragments[elem], sm.bg, elem)
                    assert new_cov >= self._coverage[elem]
                    self._coverage[elem] = new_cov

            coverage = sum(self._coverage.values())/len(self._coverage)
            for i, elem in enumerate(sorted(self._coverage)):
                self.history[i].append(self._coverage[elem])
            return "{:.1%}".format(coverage)
    @staticmethod
    def parseValue(stri):
        assert stri[-1]=="%"
        return float(stri[:-1])/100

        
        
###################################################################################################
##########
###################################################################################################
_statisticsDefaultOptions={
    "showtime": "now",
    "rmsd":True,
    "drmsd":False,
    "acc":True,
    "rog":True,
    "name": True,
    "length": True,
    "step":True,
    "local_coverage":True,
    "extreme_energies":"lh",
    "reference_energy": None,
    "extreme_rmsd": True,
    "silent":False,
    "constituing_energies": "no_clash",
    "history": "all",
    "step_save": 0,
    "save_n_best": 0,
    "save_min_rmsd": 0,
    "measure" : [],
    "distance": [],
    "asphericity": True,
    "anisotropy": True
}


class SamplingStatistics:
    def __init__(self, sm_orig, energy_functions=[], stat_source=None, output_file=sys.stdout,
                 options="all", output_directory = None):
        """
        :param sm_orig: The Spatial Model against which to collect statistics. 
                        .. warning:: This should not be modified. 
                                     (Its best to create a new SpatialModel only for this 
                                     constructor with a deepcopy of the BulgeGraph used elsewhere.)
        :param enery_functions: A list of fbe.CombinedEnergy or fbe.EnergyFunction instances.
                      These energies will be evaluated with "background=False" every time
                      `update_statistics` is called.
                      If the energy used for sampling uses a background KDE, it makes sense to
                      pass a reference to this energy here as well.
        :param output_file: A opened file handle or stdout. 
                      If this is not stdout, and options["silent"] == False: Print to both 
                      stdout and the file.
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
        self.collector=None
        self.step = 0
        self.output_file = output_file

        collectors = []
        if options=="all":
            self.options=_statisticsDefaultOptions
        elif isinstance(options, dict):
            self.options=copy.copy(_statisticsDefaultOptions)
            self.options.update(options)
        else:
            raise ValueError("Options '{}' not recognized. Please use a  "
                             "dict or the string 'all'".format(options))


        if self.options["rog"]: collectors.append(ROGStatistics())
        if self.options["acc"]: collectors.append(ACCStatistics(sm_orig))
        if self.options["rmsd"]:
            r_col=RMSDStatistics(sm_orig, show_min_max = self.options["extreme_rmsd"],
                                 save_n_best=self.options["save_min_rmsd"])
            if not r_col.silent:
                collectors.append(Delimitor())
                collectors.append(r_col)

        if self.options["drmsd"]:
            dr_col=RMSDStatistics(sm_orig, show_min_max = False,
                                  mode="dRMSD")
            if not dr_col.silent:
                collectors.append(Delimitor())
                collectors.append(dr_col)
        
        if self.options["asphericity"]: 
            collectors.append(Delimitor())
            collectors.append(AsphericityStatistics())
        if self.options["anisotropy"]: 
            collectors.append(AnisotropyStatistics())
            
        if self.options["distance"]:
            collectors.append(Delimitor())
            for n1,n2 in self.options["distance"]:
                collectors.append(Distance(n1,n2))

        collectors.append(Delimitor())
        
        if self.options["local_coverage"]:
            collectors.append(LocalFragmentCoverage(stat_source, sm_orig.bg))
            collectors.append(Delimitor())

        for ef in energy_functions:
            collectors.append(EnergyTracking(ef))
        if energy_functions:
            collectors.append(Delimitor())

        for m in self.options["measure"]: #This has to be AFTER tracking energies!
            collectors.append(EnergyMeasure(m))
            
        if self.options["showtime"] is not False:
            collectors.append(ShowTime(self.options["showtime"]))
        
        #Note: SortedCollection(maxlen=0) is valid!
        # should contain tuples: `(bg, energy)`
        self.best_cgs = SortedCollection(key=lambda x: x[1], maxlen=self.options["save_n_best"]) 

        
        self.collector = CombinedStatistics(collectors)

        if output_directory is None:
            self.out_dir = conf.Configuration.sampling_output_dir
        else:
            self.out_dir = output_directory
            if not os.path.exists(self.out_dir):
                os.makedirs(self.out_dir)
                log.info ("Directory {} created.".format(self.out_dir))

    def printline(self, line):
        """Print to both STDOUT and the log file."""
        if self.output_file != sys.stdout and not self.options["silent"]:
            print (line)

        if self.output_file != None:
            print(line, file=self.output_file)
            self.output_file.flush()

    def print_header(self):
        """
        Print the header line
        """
        contribs=["Step\tSampling_Energy"]
        if self.options["constituing_energies"]:
            contribs.append("( Constituing_Energies )")
        contribs.append(self.collector.header_str)
        contribs.append("Sampling Move")
        self.printline("\t".join(contribs))

    def update_statistics(self, sm, energy, member_energies=[], change=""):
        """
        Add a new structure to the statistics.
      
        :param sm: The spatial model.
        :param energy: The energy used for sampling of the structure, or None. 
        :param member_energies: A list of tuples `(energy_shortname, value)`
        """
        self.step+=1
        line=["{:6d}\t{:10.3f}".format(self.step, energy)]
        if self.options["constituing_energies"]=="no_clash":
            ignore_names=[fbe.RoughJunctionClosureEnergy().shortname,
                         fbe.StemVirtualResClashEnergy().shortname]
        else:
            ignore_names=[]
        line.append("( "+" ".join("{} {:10.3f}".format(*x) for x in member_energies 
                                      if x[0] not in ignore_names)+" )")
        line.append(self.collector.update(sm, self.step))
        line.append(change)
        self.printline("\t".join(line))
        
        if self.best_cgs.can_insert_right((None, energy)):
            self.best_cgs.insert_right((sm.bg.to_cg_string(), energy))
        
        if self.step % 10 == 0:
            for i, cg_stri in enumerate(self.best_cgs):
                with open(os.path.join(self.out_dir, 
                                      'best{:d}.coord'.format(i)), 'w') as f:
                    f.write(cg_stri[0])
                              
        if self.options["step_save"]>0 and self.step % self.options["step_save"] ==0:
            cg_stri = sm.bg.to_cg_string()
            with open(os.path.join(self.out_dir, 
                              'step{:06d}.coord'.format(self.step)), "w") as f:
                f.write(cg_stri)


class OutfileParser(object):
    def __init__(self):
        self._collectors = None
    def _get_lookup_table(self):        
        """
        Return a dictionary with headers as keys and 
        StatisticCollector classes (not instances) as values
        """
        return {h : cls for cls in StatisticsCollector.__subclasses__() if not cls == CombinedStatistics for h in cls.header }
    
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
        with open(filepath) as  file:
            headers = None
            data = None
            for line in file:
                line=line.strip()
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
                    data=[]
                    for i in range(len(headers)):
                        data.append([])
                else:
                    fields = line.split('\t')
                    for i, field in enumerate(fields):
                        if i == 0: #Step
                            data[i].append(int(field))
                        elif i==1: #Sampling_Energy
                            data[i].append(float(field))
                        cls = self._collectors[i]
                        if cls is not None:
                            data[i].append(cls.parse_value(field))
            data_dic = {}
            for i, header in enumerate(headers):
                if data[i]:
                    if isinstance(data[i][0], tuple):
                        data_dic["{}_{}".format(header, data[i][0][0])] = [ x[1] for x in data[i]]
                    else:
                        data_dic[header] = data[i]
        return data_dic
