#!/usr/bin/python
from __future__ import absolute_import, division, print_function, unicode_literals
from builtins import (ascii, bytes, chr, dict, filter, hex, input, #pip install future
                      int, map, next, oct, open, pow, range, round,
                      str, super, zip)
from future.builtins.disabled import (apply, cmp, coerce, execfile,
                             file, long, raw_input, reduce, reload,
                             unicode, xrange, StandardError)
__metaclass__=type
import itertools as it
import random
import logging
import inspect
from ..aux.utils import get_all_subclasses

log=logging.getLogger(__name__)

class Mover:
    HELPTEXT = ("{:25} Randomly replace a single fragment\n"
                "{:25} of the RNA.".format("Mover", ""))
    def __init__(self, stat_source):
        """
        A Mover class that randomly pickes an element and changes the stats for it.
        """
        self.stat_source = stat_source
        #: A list of tuples  (elemenmt_name, stat)
        self._prev_stats = None
        
    def _get_elem_and_stat(self, sm):
        possible_elements = sm.bg.get_mst()
        elem = random.choice(list(possible_elements))
        new_stat = self.stat_source.sample_for(sm.bg, elem)
        return elem, new_stat
    
    def move(self, sm):
        """
        Propose a single Monte Carlo move.
        
        :param sm: The spatial model
        """
        elem, new_stat = self._get_elem_and_stat(sm)
        self._prev_stats = {}
        movestring = self._move(sm, elem, new_stat)
        sm.new_traverse_and_build(start = elem, include_start = True)
        return movestring
        
    def _move(self, sm, elem, new_stat):
        prev_stat = sm.elem_defs[elem]
        self._prev_stats[elem] = prev_stat
        sm.elem_defs[elem]=new_stat        
        return "{}:{}->{};".format(elem, prev_stat.pdb_name, new_stat.pdb_name)

    def revert(self, sm):
        """
        Revert the last Monte Carlo move performed by this mover.
        """
        if self._prev_stats is None:
            raise RuntimeError("No (more) step(s) to revert.")
        for elem, stat in self._prev_stats.items():
            sm.elem_defs[elem] = stat
        self._prev_stats = {}
        sm.new_traverse_and_build(start='start')

class ExhaustiveMover(Mover):
    HELPTEXT = ("{:25} Try all fragment combinations for \n"
                "{:25} the given coarse-grained element(s).".format("ExhaustiveMover(elem[,elem,..])", ""))
    def __init__(self, stat_source, elems_of_interest, sm):
        """   
        Explore the conformation space exhaustively by iterating 
        through the stats for the given elements.
        
        :param stat_source: A StatContainer instance
        :param sm: The spatial model that will be explored
        :param elems_of_interest: A list of cg-element names.
                        The list gives the order, in which the structure space 
                        is explored in a bredth first search.
                        Example: "m0,i1,i2" Try all values for m0, then change 
                        i1 once and again try all values of m0, etc
                        Alternatively, a comma-seperated string of element names is supported.
                        
        ..warning ::
        
            If fallback-stats are used (and there are more fallback-stats than needed),
            the search is not really exhaustive but samples from the fallback stats.
        """
        super(ExhaustiveMover, self).__init__(stat_source)
        if isinstance(elems_of_interest, str):
            elems_of_interest = elems_of_interest.split(",")
        print(elems_of_interest)
        self._element_names = elems_of_interest
        self._move_iterator = self._iter_moves(sm)
    def _iter_moves(self, sm):
        for elem in self._element_names:
            for stat in self.stat_source.iterate_stats_for(sm.bg, elem):
                yield elem, stat
    def _get_elem_and_stat(self, sm):
        return next(self._move_iterator)

class NElementMover(Mover):
    HELPTEXT = ("{:25} In every move, randomly replace \n"
                "{:25} n different fragments.".format("NElementMover(N)", ""))

    def __init__(self, stat_source, n=2):
        """
        Change more than one element per iteration step.
        
        :param n: INT, How many elements should be moved in each step.
        """
        super(NElementMover, self).__init__(stat_source)
        self._n_moves = int(n)
    def move(self, sm):
        if self._n_moves>len(sm.bg.defines):
            raise ValueError("CG RNA has fewer elements than should be change per move.")
        self._prev_stats = {}
        movestring = []
        i=0
        while len(self._prev_stats)<self._n_moves:
            i+=1
            elem, new_stat = self._get_elem_and_stat(sm)
            if elem not in self._prev_stats:
                movestring.append(self._move(sm, elem, new_stat))
            if i>100000: 
                raise RuntimeError("Caught in an endless loop (?)")
        sm.new_traverse_and_build(start = "start")
        return "".join(movestring)

class OneOrMoreElementMover(NElementMover):
    HELPTEXT = ("{:25} In every move, randomly replace \n"
                "{:25} 1 to n fragments.".format("OneOrMoreElementMover[N]", ""))
    def __init__(self, stat_source, max_n=2):
        """
        Change more than one element per iteration step.
        
        :param n: INT, How many elements should be moved in each step.
        """
        super(OneOrMoreElementMover, self).__init__(stat_source)
        self._possible_n = list(range(1, int(max_n)+1))
    @property
    def _n_moves(self):
        return random.choice(self._possible_n)
    @_n_moves.setter
    def _n_moves(self, val):
        # The setter is called in the init function of the superclass.
        # Just ignore it
        pass
        
class ConnectedElementMover(NElementMover):
    """
    Change more than 1 connected elements per move step
    """
    HELPTEXT = ("{:25} In every move, randomly replace \n"
                "{:25} n adjacent fragments.".format("ConnectedElementMover(N)", ""))

    def _get_elem_and_stat(self, sm):
        if not self._prev_stats:
            return super(ConnectedElementMover, self)._get_elem_and_stat(sm)
        #Get all cg-elements connected to the changed element that have not yet been changed
        neighbors = { d for changed_elem in self._prev_stats for d in sm.bg.edges[changed_elem] if d not in self._prev_stats and d in sm.bg.mst}
        elem = random.choice(list(neighbors))
        new_stat = self.stat_source.sample_for(sm.bg, elem)
        return elem, new_stat

class WholeMLMover(Mover):
    """
    Pick a whole multiloop or a random other element and change the stats for it.
    
    Changing a whole multiloops means changing more than one stat. 
    The stats for all ml-segments of this multiloop are sampled independently.
    """
    HELPTEXT = ("{:25} In every move, either replace \n"
                "{:25} a single not-multiloop fragment,\n"
                "{:25} or all multiloop segments \n"
                "{:25} for a single multiloop.".format("WholeMLMover", "", "", ""))

    def move(self, sm):
        movestring = []
        self._prev_stats = {}
        while self._has_incomplete_ml(sm):
            elem, new_stat = self._get_elem_and_stat(sm)
            movestring.append(self._move(sm, elem, new_stat))
        sm.new_traverse_and_build(start = "start")
        return "".join(movestring)
    def _has_incomplete_ml(self, sm):
        if not self._prev_stats:
            return True
        if any(key[0]=="m" for key in self._prev_stats):
            log.info("_has_incomplete_ml: m in prev_stats {}".format(self._prev_stats.keys()))
            return self._get_missing_ml_nodes(sm)
        else:
            return False
    def _get_missing_ml_nodes(self, sm):
        sampled = set(k for k in self._prev_stats)
        whole_ml = set(sm.bg.shortest_mlonly_multiloop(list(sampled)[0]))
        whole_ml = whole_ml & sm.bg.mst
        missing = whole_ml - sampled
        log.info("_get_missing_ml_nodes: sampled {}, mst {}, missing {}".format(sampled, sm.bg.mst, missing))
        return missing

    def _get_elem_and_stat(self, sm):
        if not self._prev_stats:
            return super(WholeMLMover,self)._get_elem_and_stat(sm)
        else:
            elem = random.choice(list(self._get_missing_ml_nodes(sm)))
            new_stat = self.stat_source.sample_for(sm.bg, elem)
            return elem, new_stat

class MlRelaxationMover(Mover):
    """
    Change a random element and if it ewas a multiloop segment, change the other multiloop 
    segments of this multiloop until the loop closure energy is zero.
    """
    HELPTEXT = ("{:25} Not yet implemented\n".format("MlRelaxationMover"))
    def __init__(self, stat_source):
        raise NotImplementedError()
        



####################################################################################################
### Command line parsing
####################################################################################################

def mover_from_string(stri, stat_source, sm=None):
    """
    Return a single Mover instance from a string describing the mover.
    
    The string needs to contain the class name of the mover and optionally
    may contain one argument in square brackets: MOVERNAME[ARGUMENT]
    
    :param stat_source: The stat_container that shall be used for all moves.
    :param sm: The Spatiel model. Needed, for ExhaustiveMover
    
    """
    if '[' in stri:
        movername, option = stri.split('[')
        if option[-1]!="]":
            raise ValueError("Missing closing bracket at end of '[{}'".format(option))
        args = [option[:-1]]
    else:
        args = []
        movername = stri
    kwargs = {}
    mover_classes = { cls.__name__: cls for cls in get_all_subclasses(Mover, include_base = True) }
    try:
        cls = mover_classes[movername]
        required = inspect.getargspec(cls.__init__).args 
        if "sm" in required:
            kwargs["sm"] = sm
        mover = cls(stat_source, *args, **kwargs)
    except TypeError as e: #Give more helpful error messages
        if "arguments" not in e.message or "__init__" not in e.message:
                raise
        # Something like TypeError: __init__() takes exactly 4 arguments (3 given)
        signature = inspect.getargspec(cls.__init__)
        called = "[self], stat_source,"
        called+=",".join(map(str,args))
        if kwargs:
            called+=","
            for k,v in kwargs.items():
                called+="{}={},".format(k,v)
        raise TypeError("__init__ with signature ({}) cannot be called with ({})".format(signature, called))

    return mover

def get_argparse_help(): 
    """
    Return a pre-formatted string that can be passed to "help" in argparse.add_argument,
    if the resulting argument is parsed with `mover_from_string`
    
    .. warning::
    
        Please specify the same default in `add_argument` as described in this help-text
    """
    help= ("Which types of Moves to use during sampling.\n"
           "Default: Mover\n"
           "One of the following:\n")
    for mover_class in get_all_subclasses(Mover, include_base = True):
        help+=mover_class.HELPTEXT+"\n"
    return help
