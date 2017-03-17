#!/usr/bin/python
from __future__ import absolute_import, division, print_function, unicode_literals
from builtins import (ascii, bytes, chr, dict, filter, hex, input,
                      int, map, next, oct, open, pow, range, round,
                      str, super, zip) #future package
from future.builtins.disabled import (apply, cmp, coerce, execfile,
                             file, long, raw_input, reduce, reload,
                             unicode, xrange, StandardError)
import forgi.threedee.model.stats as ftmstats
import sys
import random
from collections import defaultdict
import logging
log = logging.getLogger(__name__)
try:
    from functools import lru_cache #python 3.2
except ImportError:
    try:
        from backports.functools_lru_cache import lru_cache #pip install backports.functools_lru_cache
    except ImportError:
        lru_cache = lambda *args, **kwargs: lambda x: x #No-op decorator taking arguments


def parse_stats_file(file_handle):
    stats = {"stem": defaultdict(list),
             "angle": defaultdict(list),
             "loop": defaultdict(list), "3prime": defaultdict(list), "5prime": defaultdict(list)}
    for line in file_handle:
        line=line.strip()
        if "#" in line:
            line = line.split('#')[0]
        if not line: 
            continue
        if line.startswith("stem"):
            stem_stat = ftmstats.StemStat(line)
            stats["stem"][stem_stat.bp_length].append(stem_stat)
        elif line.startswith("angle"):
            angle_stat = ftmstats.AngleStat()
            try:
                angle_stat.parse_line(line)
            except:
                print("Could not parse line '{}'".format(line), file=sys.stderr)
                raise
            if len(angle_stat.define) > 0 and angle_stat.define[0] == 1: #An angle at the beginning of a structure
                #I guess this should never happen, if the stats do not stem from faulty bulge graphs.
                log.error("Ignoring angle stat {} because it is at the beginning of a structure."
                          " Does the stat come from a faulty BulgeGraph?".format(angle_stat.pdb_name))
                continue
            stats["angle"][(angle_stat.dim1, angle_stat.dim2, angle_stat.ang_type)].append(angle_stat)
            # Adding the reverse does not work as intended and produces a lot of structures 
            # that do not fulfill the constraint energy.
            # stats["angle"][(angle_stat.dim1, angle_stat.dim2, -angle_stat.ang_type)].append(angle_stat)
        else:
            key = line.split()[0]
            if key not in ["3prime", "5prime", "loop"]:
                raise ValueError("Illegal line in stats file: '{}'".format(line))
            stat = ftmstats.LoopStat(line)
            stats[key][stat.bp_length].append(stat)
    return stats

def read_stats_file(filename):
    with open (filename) as f:
        return parse_stats_file(f)
    
letter_to_stat_type = {
    "s": "stem",
    "h": "loop",
    "t": "3prime",
    "f": "5prime",
    "m": "angle",
    "i": "angle"
}
    
def _key_from_bg_and_elem(bg, elem):        
    dims = bg.get_node_dimensions(elem)
    if elem[0] in "i, m":
        ang_type = bg.get_angle_type(elem)
        return tuple([dims[0], dims[1], ang_type])
    else:
        return dims[0]
    
class StatStorage(object):
    def __init__(self, filename, fallback_filenames = None):
        self.filename = filename
        if fallback_filenames is None:
            fallback_filenames = []
        self.fallbacks = fallback_filenames
        self._sources = None
        self._has_warned = set() #Only emit warnings about insufficient stats once.
    def _iter_stat_sources(self):                
        if self._sources is None:
            self._sources = [read_stats_file(self.filename)]
        for i in range(len(self.fallbacks)+1):
            if i>=len(self._sources):
                self._sources.append(read_stats_file(self.fallbacks[i-1]))
            yield self._sources[i]
        raise StopIteration
    
    @lru_cache(maxsize = 128)
    def _possible_stats(self, stat_type, key, min_entries = 100):
        choose_from = []
        weights = []
        statfiles = self._iter_stat_sources()
        while len(choose_from)<min_entries:
            try:
                source = next(statfiles)[stat_type]
            except StopIteration: #All stat_files exhausted
                if (stat_type, key, min_entries) not in self._has_warned:
                    log.warning("Only {} stats found for {} with key {}".format(len(choose_from), stat_type, key))
                    self._has_warned.add((stat_type, key, min_entries))
                break
            if key in source:
                stats=source[key]
                log.debug("Appending {} stats from source.".format(len(stats)))
                choose_from.append(stats)
                if weights:
                    weights.append(min(min_entries - sum(weights), len(stats)))
                else:
                    weights.append(len(stats))
        if not choose_from:
            raise LookupError("No stats found for {} with key {}".format(stat_type, key))
        
        return weights, choose_from

    """
    def get_possible_stats(self, bg, elem, min_entries = 100):
        log.debug("Getting stats for {} of bg {}".format(elem, bg.name))
        return self._possible_stats(letter_to_stat_type[elem[0]], key, min_entries)
    """
    
    def sample_for(self, bg, elem, min_entries = 100):
        """
        Sample a stat for the given coarse-grained element elem of bulge graph bg.
        
        During sampling, fallback-files (if used) have a lower weight
        depending on the number of possible stats in the main file.
        
        .. note::
        
            A stat is a class holding an angle and distance based 
            description of a coarse-grained fragment.
            
        
        :param bg: A CoarseGrainRNA or BulgeGraph object
        :param elem: The element name, e.g. "s0"
        :param min_entries: If less than min-entries stats are found, try to use 
                    the fallback-files to gather enough stats to sample from.
        :param return: A singe Stat object, sampled from to possible stats.
        """
        key = _key_from_bg_and_elem(bg, elem)
        weights, stats = self._possible_stats(letter_to_stat_type[elem[0]], key, min_entries)
        r = random.randrange(sum(weights))
        for i, w in enumerate(weights):
            if r<w:
                try:
                    return random.choice(stats[i])
                except:
                    print("Stat file {}. weight = {}, r={}".format(i, w, r))
                    raise
                
            r-=w
        assert False
    def iterate_stats_for(self, bg, elem, min_entries = 100, cycle = False):
        """
        Iterate over all stats for the given element.
        
        If fallback-stats are present, a sample of the fallback_stats is 
        generated every time the iterator is created and after iterating the 
        main stats iteration continues over the sample of the fallback stats.
        """
        key = _key_from_bg_and_elem(bg, elem)
        weights, stats = self._possible_stats(letter_to_stat_type[elem[0]], key, min_entries)
        stat_samples = []
        for i, w in enumerate(weights):
            stat_samples.append(random.sample(stats[i], w))
        while True:
            for stats in stat_samples:
                for s in stats:
                    yield s
            if not cycle:
                break #Exhaust the generator
    
    @lru_cache(maxsize = 128)
    def _get_pdbname_statsets_for(self, elem, key, min_entries = 100):
        """
        Used in self.coverage_for, optimized by the use of lru_cache
        """
        weights, stats = self._possible_stats(letter_to_stat_type[elem[0]], key, min_entries)
        stat_sets = []
        for i, w in enumerate(weights):
            stat_set = set(stat.pdb_name for stat in stats[i])
            if len(stat_set) != len(stats[i]):
                raise ValueError("coverage_for only works if the pdb_names of all "
                                 "stats are unique.")
            stat_sets.append((w, stat_set))
        return sum(weights), stat_sets
        
    def coverage_for(self, sampled_stat_names, bg, elem, min_entries = 100):
        """
        Count, what percentage of the stats have been used during sampling.
        
        :param sampled_stat_names: A set of pdb_names of stats.
        
        For the other parameters, see `self.sample_for`
        """
        
        key = _key_from_bg_and_elem(bg, elem)
        total_weight, stat_sets = self._get_pdbname_statsets_for(elem, key, min_entries)
        coverage = 0.
        for weight, stat_set in stat_sets:
            found = sampled_stat_names & stat_set
            coverage += len(found)/len(stat_set) * weight/total_weight
        return coverage
