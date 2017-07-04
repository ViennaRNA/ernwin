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
import string
from logging_exceptions import log_to_exception


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
        elif line.startswith("angle") or line.startswith("open") or line.startswith("pseudo"):
            angle_stat = ftmstats.AngleStat()
            try:
                angle_stat.parse_line(line)
            except Exception as e:
                with log_to_exception(log, e):
                    log.error("Could not parse file due to error parsing line '{}'".format(line))
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
            # Note that CoarseGrainRNA.get_stats extracts two angle stats er angle.
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


class StatStorage(object):
    def __init__(self, filename, fallback_filenames = None):
        self.filename = filename
        if fallback_filenames is None:
            fallback_filenames = []
        self.fallbacks = fallback_filenames
        self._sources = None
        self._has_warned = set() #Only emit warnings about insufficient stats once.

    @staticmethod
    def key_from_bg_and_elem(bg, elem):
        dims = bg.get_node_dimensions(elem)
        if elem[0] in "i, m":
            ang_type = bg.get_angle_type(elem, allow_broken = True)
            return tuple([dims[0], dims[1], ang_type])
        elif elem[0]=="h" and dims[0]<3:
            return 3 #Hairpins<3 probably means missing residues. Just return the smalles possible dimension
        else:
            return dims[0]

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
        """
        :returns: Two lists, `weights` and `choose_from` of the same length.
                  weights is a list of floats, choose_from is a list of stats.
        """
        choose_from = []
        weights = []
        statfiles = self._iter_stat_sources()
        while sum(weights)<min_entries:
            try:
                source = next(statfiles)[stat_type]
            except StopIteration: #All stat_files exhausted
                if (stat_type, key, min_entries) not in self._has_warned:
                    log.warning("Only {} stats found for {} with key {}".format(len(choose_from), stat_type, key))
                    self._has_warned.add((stat_type, key, min_entries))
                break
            if key in source:
                stats=source[key]
                num_stats = len(stats)
                log.debug("Appending {} stats from source.".format(len(stats)))
                choose_from.extend(stats)
                if not weights:
                    weight = 1 #All stats from the first stat_source always has weight 1, even if there are more than min_entries stats.
                else:
                    remaining_total_weight = min_entries - sum(weights)
                    weight = min(1, remaining_total_weight/len(stats))
                weights += [weight]*num_stats
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
        key = self.key_from_bg_and_elem(bg, elem)
        log.debug("Calling _possible_stats with %r, %r", letter_to_stat_type[elem[0]], key)
        weights, stats = self._possible_stats(letter_to_stat_type[elem[0]], key, min_entries)
        r = random.uniform(0, sum(weights))
        for i, w in enumerate(weights):
            if r<=w:
                return stats[i]
            r-=w
        assert False
        return stats[0] #Fallback if asserts are disabled. Should be unreachable.

    def iterate_stats(self, elem, key, min_entries = 100, cycle = False):
        weights, stats = self._possible_stats(letter_to_stat_type[elem[0]], key, min_entries)
        stat_samples = []
        for i, w in enumerate(weights):
            r = random.random()
            if r<=w:
                stat_samples.append(stats[i])
        while True:
            for stat in stat_samples:
                yield stat
            if not cycle:
                break #Exhaust the generator

    def iterate_stats_for(self, bg, elem, min_entries = 100, cycle = False):
        """
        Iterate over all stats for the given element.

        If fallback-stats are present, a sample of the fallback_stats is
        generated every time the iterator is created and after iterating the
        main stats iteration continues over the sample of the fallback stats.
        """
        key = self.key_from_bg_and_elem(bg, elem)
        return self.iterate_stats(letter_to_stat_type[elem[0]], key, min_entries)


    def coverage_for(self, sampled_stat_names, bg, elem, min_entries = 100):
        """
        Count, what percentage of the stats have been used during sampling.

        :param sampled_stat_names: A set of pdb_names of stats.

        For the other parameters, see `self.sample_for`
        """

        key = self.key_from_bg_and_elem(bg, elem)
        weights, stats = self._possible_stats(letter_to_stat_type[elem[0]], key, min_entries)
        total_weight = sum(weights)
        coverage = 0.
        for i, weight in enumerate(weights):
            if stats[i].pdb_name in sampled_stat_names:
                coverage += weight/total_weight
        return coverage

def identitical_bases(seq1, seq2):
    if len(seq1) != len(seq2):
        # This is an interim solution while the stats still contain adjacent nucleotides for stems.
        # If stem stats have adjacent=False and all others adjacent=True, we can raise the ValueError again
        if len(seq2)==len(seq1)+2:
            seq2 = seq2[1:-1]
        else:
            log.debug("Sequences do not have the same length: %s %s.", seq1, seq2)
            return 0
            # raise ValueError("Pairwise identity is only defined for strings of same length. "
            #                  "Found {} and {}. Seq1={}".format(len(seq1), len(seq2), seq1[:10]))
    s = 0
    for i in range(len(seq1)):
        s+= (seq1[i]==seq2[i])
    return s

def seq_and_pyrpur_similarity(sequences, stat):
    translation = {ord(c): ord(t) for c, t in zip(u"AGCU", u"RRYY")}
    ib4 = 0
    ib2 = 0
    for i in range(len(sequences)): #single or double stranded
        #identical bases 4-letter and 2-letter alphabeth
        ib4 += identitical_bases(sequences[i], stat.seqs[i])
        ib2 += identitical_bases(sequences[i].translate(translation), stat.seqs[i].translate(translation))
    score = (ib2 + ib4 +1) / (sum(len(x) for x in sequences) * 2 +1) #^2 for ib4 and ib2
    #log.debug("Score is %f", score)
    return score

class SequenceDependentStatStorage(StatStorage):
    def __init__(self, filename, fallback_filenames = None, sequence_score = seq_and_pyrpur_similarity):
        self.sequence_score = sequence_score
        super(SequenceDependentStatStorage, self).__init__(filename, fallback_filenames)

    @staticmethod
    def key_from_bg_and_elem(bg, elem):
        dims = bg.get_node_dimensions(elem)
        if elem[0] in "i, m":
            ang_type = bg.get_angle_type(elem, allow_broken = True)
            return tuple([dims[0], dims[1], ang_type]), tuple(bg.get_define_seq_str(elem, adjacent = True))
        else:
            return dims[0], tuple(bg.get_define_seq_str(elem, adjacent = elem[0]!="s"))


    @lru_cache(maxsize = 256)
    def _possible_stats(self, stat_type, key, min_entries = 100):
        """
        :returns: Two lists, `weights` and `choose_from` of the same length.
                  weights is a list of floats, choose_from is a list of stats.
        """
        key, sequence = key
        choose_from = []
        weights = []
        statfiles = self._iter_stat_sources()
        while sum(weights)<min_entries:
            try:
                source = next(statfiles)[stat_type]
            except StopIteration: #All stat_files exhausted
                if (stat_type, key, min_entries) not in self._has_warned:
                    log.warning("Only {} stats found for {} with key {}".format(len(choose_from), stat_type, key))
                    self._has_warned.add((stat_type, key, min_entries))
                break
            if key in source:
                stats=source[key]
                num_stats = len(stats)
                choose_from.extend(stats)
                if not weights:
                    weight = 1
                else:
                    remaining_total_weight = min_entries - sum(weights)
                    weight = min(1, remaining_total_weight/len(stats))
                for stat in stats:
                    #log.debug("Evaluating score for %s %s %s", key, sequence, stat)
                    weights.append(weight*self.sequence_score(sequence, stat))
        if not choose_from:
            raise LookupError("No stats found for {} with key {}".format(stat_type, key))

        log.info("For key %s: %s stats with weight in range %s-%s", key, len(choose_from), min(weights), max(weights))
        if len(choose_from)<20:
            for i, w in enumerate(weights):
                log.debug("     Weight %f: %s", w, choose_from[i].seqs)
        return weights, choose_from
