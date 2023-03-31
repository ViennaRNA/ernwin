#!/usr/bin/python
from __future__ import absolute_import, division, print_function, unicode_literals
from builtins import (ascii, bytes, chr, dict, filter, hex, input,
                      int, map, next, oct, open, pow, range, round,
                      str, super, zip)  # future package
import sys
import random
import math
from collections import defaultdict
import logging
import string
import os.path as op

import numpy as np

from logging_exceptions import log_to_exception

import forgi.threedee.model.stats as ftmstats

from fess import data_file
from fess.builder import config
import fess.motif.annotate as fma
from six.moves import range
from six.moves import zip

log = logging.getLogger(__name__)
try:
    from functools import lru_cache #python 3.2
except ImportError:
    try:
        from backports.functools_lru_cache import lru_cache #pip install backports.functools_lru_cache
    except ImportError:
        lru_cache = lambda *args, **kwargs: lambda x: x #No-op decorator taking arguments


def key_to_human_readable(key):
    try:
        if key[2]==1:
            return "interior loop with {} and {} unpaired nucleotides".format(key[0], key[1])
        elif key[2]==6:
            return "multiloop segment with {} nucleotides".format(key[0])
    except Exception:
        pass
    return "key {}".format(key)


def patch_angtype(ang_type):
    """
    Instead of having 3 angle_types 2,3 and 4 for junjtion segments, use one generic angle type '6'
    """
    if 2<=abs(ang_type)<=3:
        return math.copysign(6, ang_type)
    if abs(ang_type)in [4,5]:
        return -math.copysign(6, ang_type)
    return ang_type

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
            angle_stat.ang_type = patch_angtype(angle_stat.ang_type)
            log.debug("Reading angle_stat with dimensions %s and %s, and type %s. With define %s", angle_stat.dim1, angle_stat.dim2, angle_stat.ang_type, angle_stat.define)
            stats["angle"][(angle_stat.dim1, angle_stat.dim2, angle_stat.ang_type)].append(angle_stat)
            # Adding the reverse does not work as intended and produces a lot of structures
            # that do not fulfill the constraint energy.
            # stats["angle"][(angle_stat.dim1, angle_stat.dim2, -angle_stat.ang_type)].append(angle_stat)
            # Note that CoarseGrainRNA.get_stats extracts two angle stats per angle.
        else:
            key = line.split()[0]
            if key not in ["3prime", "5prime", "loop"]:
                raise ValueError("Illegal line in stats file: '{}'".format(line))
            stat = ftmstats.LoopStat(line)
            stats[key][stat.bp_length].append(stat)
    return stats

def read_stats_file(filename):
    log.info("Reading stats-file %s", filename)
    with open (filename) as f:
        try:
            return parse_stats_file(f)
        except Exception as e:
            with log_to_exception(log, e):
                log.error("Failed to parse file %s", filename)
            raise

letter_to_stat_type = {
    "s": "stem",
    "h": "loop",
    "t": "3prime",
    "f": "5prime",
    "m": "angle",
    "i": "angle"
}

def make_continuous( discrete_angle_stats):
    '''
    Create a kernel density estimation of the statistics
    for each bulge dimension represented in the collection
    of discrete angle statistics.

    Each KDE will have six dimensions corresponding to the six
    dimensions necessary for describing the orientation of one
    helix with respect to another.

    :param discrete_angle_statistics: A dictionary of dictionaries,
        each one containing and AngleStats structure.
    '''
    import scipy.stats as ss
    data=[]
    for d in discrete_angle_stats:
        data += [[d.u, d.v, d.t, d.r1, d.u1, d.v1]]

    if len(data) < 3:
        raise ValueError("Insufficient stats to make continuouse")

    try:
        return ss.gaussian_kde(np.array(data).T)
    except np.linalg.LinAlgError as lae:
        print("Singular matrix, dimensions:", dims, file=sys.stderr)

class ContinuouseStatSampler:
    def __init__(self, all_stats, key):
        kde = make_continuous(all_stats)
        self.kde = kde
        self.key = key
        data=[]
        for d in all_stats:
            data += [[d.u, d.v, d.t, d.r1, d.u1, d.v1]]
        self.data = np.array(data)
    def sample(self):
        import scipy.stats as ss
        log.debug("DATA %s with shape %s", self.data, self.data.shape)
        r1 = ss.gaussian_kde(self.data[:, 3]).resample(1)[0][0]
        log.debug("r1 is %s", r1)
        rnd = np.random.rand(5) # 5 random values from 0 to 1
        u = rnd[0]*np.pi
        v = rnd[1]*2*np.pi-np.pi
        t = rnd[2]*2*np.pi-np.pi
        u1 = rnd[3]*np.pi
        v1 = rnd[4]*2*np.pi-np.pi

        #sample = self.kde.resample(1).T
        log.debug("continuouse stat sample %s", (u,v,t,r1,u1,v1))
        #u,v,t,r1,u1,v1 = sample[0]
        stat = ftmstats.AngleStat(stat_type="angle", pdb_name='cont-{:.1f}_{:.1f}_{:.1f}_{:.1f}_{:.1f}_{:.1f}'.format(u,v,t,r1,u1,v1),
                                dim1=self.key[0], dim2=self.key[1], u=u, v=v, t=t, r1=r1, u1=u1, v1=v1,
                                ang_type=self.key[2], define=[], seq="", vres={})
        return stat


class StatStorage(object):
    def __init__(self, filename, fallback_filenames = None, continuouse=None, blacklist=[]):
        self.filename = filename
        if fallback_filenames is None:
            fallback_filenames = []
        self.fallbacks = fallback_filenames
        self._sources = None
        self._has_reported = set() #Only emit warnings about insufficient stats once.
        if continuouse:
            if any(elem[0] not in "mi" for elem in continuouse):
                raise ValueError("Continuouse stats are only allowed for interior and multiloops.")
            self.continuouse = continuouse[:]
        else:
            self.continuouse = []
        self.blacklist=blacklist

    def in_blacklist(self, stat):
        statname = stat.pdb_name
        for pattern in self.blacklist:
            if pattern in statname:
                return True
        return False


    @staticmethod
    def key_from_bg_and_elem(bg, elem):
        dims = bg.get_node_dimensions(elem, with_missing=(elem[0]!="s"))
        if elem[0] in "i, m":
            ang_type = bg.get_angle_type(elem, allow_broken = True)
            ang_type = patch_angtype(ang_type)
            key = tuple([dims[0], dims[1], ang_type])
        elif elem[0]=="h" and dims[0]<3:
            key = 3 #Hairpins<3 probably means missing residues. Just return the smalles possible dimension
        else:
            key = dims[0]
        log.debug("Key for bg %s and element %s is %s (node dimensions without_missing "
                  "%s, with missing %s)", bg.name, elem, key,
                  bg.get_node_dimensions(elem, with_missing=False),
                  bg.get_node_dimensions(elem, with_missing=True))
        return key

    def _iter_stat_sources(self):
        if self._sources is None:
            self._sources = [read_stats_file(self.filename)]
        for i in range(len(self.fallbacks)+1):
            if i>=len(self._sources):
                self._sources.append(read_stats_file(self.fallbacks[i-1]))
            yield self._sources[i]

    @lru_cache(maxsize = 128)
    def _possible_stats(self, stat_type, key, min_entries=100, strict=False, enable_logging=True):
        try:
            return self._possible_stats_inner(stat_type, key, min_entries, strict, enable_logging=enable_logging)
        except RuntimeError as e:
            if "recursion" not in str(e):
                raise
            # Maximum recursion depth exceeded.
            raise LookupError("No stats found for {} with key {}, even after reducing length.".format(stat_type, key))

    def _possible_stats_inner(self, stat_type, key, min_entries = 100, strict=False, enable_logging=True):
        """
        :returns: Two lists, `weights` and `choose_from` of the same length.
                  weights is a list of floats, choose_from is a list of stats.
        """
        has_reported=False
        choose_from = []
        weights = []
        statfiles = self._iter_stat_sources()
        while sum(weights)<min_entries:
            try:
                sf = next(statfiles)
                source = sf[stat_type]
            except StopIteration: #All stat_files exhausted
                if enable_logging and len(choose_from)<min_entries and (stat_type, key, min_entries) not in self._has_reported:
                    log.warning("Only {} {}-stats found for {}".format(len(choose_from), stat_type, key_to_human_readable(key)))
                    has_reported = True
                break
            if key in source:
                stats=[ stat for stat in source[key] if not self.in_blacklist(stat) ]
                num_stats = len(stats)
                if enable_logging:
                    log.info("Appending {} stats from source {} for {}.".format(len(stats),id( sf), key))
                choose_from.extend(stats)
                if not weights:
                    weight = 1 #All stats from the first stat_source always has weight 1, even if there are more than min_entries stats.
                else:
                    remaining_total_weight = min_entries - sum(weights)
                    weight = min(1, remaining_total_weight/len(stats))
                weights += [weight]*num_stats
            else:
                if enable_logging:
                     log.info("Appending NO stats from source {} for {}.".format(id( sf), key))
        if enable_logging and (stat_type, key, min_entries) not in self._has_reported:
            log.info("Found {} stats for {} with key {}".format(len(choose_from), stat_type, key))
            has_reported = True

        if not choose_from:
            if not strict: # Fallback to the next smaller stat, until all options are exhausted.
                if stat_type == "loop" and key>3:
                    if enable_logging and (stat_type, key, min_entries) not in self._has_reported:
                        log.error("Trying %s instead of %s for %s-stat", key_to_human_readable(key-1), key_to_human_readable(key), stat_type)
                        self._has_reported.add((stat_type, key, min_entries))
                    return self._possible_stats_inner(stat_type, key-1, min_entries, strict, enable_logging)
                elif stat_type == "angle" and (key[0]>0 or 1000>key[1]>0):
                    if 1000>key[1]>key[0]:
                        new_key = (key[0], key[1]-1, key[2])
                    else:
                        new_key = (key[0]-1, key[1], key[2])
                    if enable_logging and (stat_type, key, min_entries) not in self._has_reported:
                        log.error("Trying key %s instead of %s for %s-stat", key_to_human_readable(new_key), key_to_human_readable(key), stat_type)
                        self._has_reported.add((stat_type, key, min_entries))
                    return self._possible_stats_inner(stat_type, new_key, min_entries, strict, enable_logging)
            # If everything else fails, raise an error even if strict was disabled.
            raise LookupError("No stats found for {} with key {}".format(stat_type, key))
        if has_reported:
            self._has_reported.add((stat_type, key, min_entries))
        return weights, choose_from

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
        :param return: A single Stat object, sampled from the possible stats.
        """
        if elem in self.continuouse:
            return self._sample_continuouse_stat(bg, elem, min_entries)
        else:
            key = self.key_from_bg_and_elem(bg, elem)
            log.debug("Calling _possible_stats with %r, %r", letter_to_stat_type[elem[0]], key)
            weights, stats = self._possible_stats(letter_to_stat_type[elem[0]], key, min_entries)
            r = random.uniform(0, sum(weights))
            for i, w in enumerate(weights):
                if r<=w:
                    # TODO: Penalize stats found with JARED, but for another loop
                    return stats[i]
                r-=w
            assert False
            return stats[0] #Fallback if asserts are disabled. Should be unreachable.

    def _sample_continuouse_stat(self, bg, elem, min_entries=100):
        key = self.key_from_bg_and_elem(bg, elem)
        kde = self._get_statsampler(letter_to_stat_type[elem[0]], key, min_entries)
        return kde.sample()

    def _iterate_continuouse_stat(self, bg, elem, min_entries=100):
        key = self.key_from_bg_and_elem(bg, elem)
        kde = self._get_statsampler(letter_to_stat_type[elem[0]], key, min_entries)
        log.debug("Starting to iterate continuouse stats")
        while True:
            yield kde.sample()

    def _get_statsampler(self, stat_type, key, min_entries):
        log.debug("Getting continuouse kde-based sampler for %s", key)
        all_stats = list(self.iterate_stats(stat_type, key, min_entries, False))
        return ContinuouseStatSampler(all_stats, key)

    def iterate_stats(self, stat_type, key, min_entries=100, cycle=False):
        weights, stats = self._possible_stats(stat_type, key, min_entries)
        stat_samples = []
        for i, w in enumerate(weights):
            r = random.random()
            if r <= w:
                stat_samples.append(stats[i])
        while True:
            for stat in stat_samples:
                yield stat
            if not cycle:
                break  # Exhaust the generator

    def load_stat_by_name(self, bg, elem, name):
        key = self.key_from_bg_and_elem(bg, elem)

        _, stats = self._possible_stats(letter_to_stat_type[elem[0]], key, min_entries=float('inf'),
                                        enable_logging=False)
        found = None
        for stat in stats:
            if stat.pdb_name == name:
                assert found is None
                found = stat
        if not found:
            raise RuntimeError("Cannot load stat {} for elem {}. Maybe a different stat file "
                               "or different cg was used?".format(name, elem))
        return found

    def iterate_stats_for(self, bg, elem, min_entries=100, cycle=False):
        """
        Iterate over all stats for the given element.

        If fallback-stats are present, a sample of the fallback_stats is
        generated every time the iterator is created and after iterating the
        main stats iteration continues over the sample of the fallback stats.
        """
        if elem in self.continuouse:
            for stat in self._iterate_continuouse_stat(bg, elem, min_entries):
                yield stat
        else:
            key = self.key_from_bg_and_elem(bg, elem)
            log.debug("Key is %s, elem is %s, elem[0] is %s", key, elem, elem[0])
            log.debug("letter_to_stat_type[elem[0]] is %s", letter_to_stat_type[elem[0]])
            for stat in self.iterate_stats(letter_to_stat_type[elem[0]], key, min_entries):
                yield stat



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
            ang_type = patch_angtype(ang_type)
            return tuple([dims[0], dims[1], ang_type]), tuple(bg.get_define_seq_str(elem, adjacent = True))
        else:
            return dims[0], tuple(bg.get_define_seq_str(elem, adjacent = elem[0]!="s"))


    @lru_cache(maxsize = 256)
    def _possible_stats(self, stat_type, key, min_entries = 100, enable_loging=True):
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
                sf = next(statfiles)
                source = sf[stat_type]
            except StopIteration: #All stat_files exhausted
                if enable_logging and (stat_type, key, min_entries) not in self._has_reported:
                    log.warning("Only {} stats found for {} with key {}".format(len(choose_from), stat_type, key))
                    self._has_reported.add((stat_type, key, min_entries))
                break
            if key in source:
                stats=[ stat for stat in source[key] if not self.in_blacklist(stat) ]
                num_stats = len(stats)
                if enable_logging:
                    log.info("Added %s stats ", num_stats)
                choose_from.extend(stats)
                if not weights:
                    weight = 1
                else:
                    remaining_total_weight = min_entries - sum(weights)
                    weight = min(1, remaining_total_weight/len(stats))
                for stat in stats:
                    #log.debug("Evaluating score for %s %s %s", key, sequence, stat)
                    weights.append(weight*self.sequence_score(sequence, stat))
            else:
                if enable_logging:
                    log.info("Nothing added from stat_source %s",id(sf))

        if not choose_from:
            raise LookupError("No stats found for {} with key {}".format(stat_type, key))
        if enable_logging:
            log.info("For key %s: %s stats with weight in range %s-%s", key, len(choose_from), min(weights), max(weights))
            if len(choose_from)<20:
                for i, w in enumerate(weights):
                    log.debug("     Weight %f: %s", w, choose_from[i].seqs)
        return weights, choose_from



def update_parser(parser):
    stat_options = parser.add_argument_group(
                            "Choosing stats",
                            description="These options control what stats ERNWIN\n"
                                        " uses for sampling")
    stat_options.add_argument('--stats-file', type=str,
                              default=data_file(config.Configuration.default_stats_file),
                              help= "A filename.\n"
                                    "A file containing all the stats to sample from\n"
                                    " for all coarse grained elements")
    stat_options.add_argument(
                        '--fallback-stats-files', nargs = '+', type=str,
                        help= "A list of fallback stats file that can be uses if insufficient stats "
                              "are found in the normal stats file for a coarse-grained element.\n"
                              "If more than one file is given, the files are used in the order specified.\n")
    stat_options.add_argument('--sequence-based', action="store_true",
                              help= "Take the sequence into account when choosing stats.")
    stat_options.add_argument('--jar3d', action="store_true",
                              help="Use JAR3D to restrict the stats \n"
                                   "for interior loops to matching motifs.\n"
                                   "Requires the correct paths to jar3d to be set in\n "
                                   "fess.builder.config.py"   )
    stat_options.add_argument('--continuouse-stats', type=str, help='A list of names of interior'
                                'or multiloops. The stats of these elements will be sampled from'
                                ' a continuouse distribution. EXPERIMENTAL, DONT USE THIS.')
    stat_options.add_argument('--blacklist-stats', type=str, help="A comma seperate list of pdb-ids. Disallow stats from these pdb ids.")

def from_args(args, cg):
    if args.sequence_based:
        StatSourceClass = SequenceDependentStatStorage
    else:
        StatSourceClass = StatStorage
    kwargs={}
    if args.continuouse_stats:
        kwargs["continuouse"] = args.continuouse_stats.split(",")
    if args.blacklist_stats:
        kwargs["blacklist"] = args.blacklist_stats.split(",")
    if args.jar3d:
        jared_out    = op.join(config.Configuration.sampling_output_dir, "jar3d.stats")
        jared_tmp    = op.join(config.Configuration.sampling_output_dir, "jar3d")
        motifs = fma.annotate_structure(cg, jared_tmp, cg.name.split('_')[0])
        fma.print_stats_for_motifs(motifs, filename = jared_out, temp_dir = config.Configuration.sampling_output_dir, args=args )
        #Here, all stats are fallback stats for the JAR3D hits.
        new_fallbacks = [args.stats_file]
        if args.fallback_stats_files is not None:
            new_fallbacks += args.fallback_stats_files
        stat_source = StatSourceClass(jared_out, new_fallbacks, **kwargs)
    else:
        stat_source = StatSourceClass(args.stats_file, args.fallback_stats_files, **kwargs)
    return stat_source
