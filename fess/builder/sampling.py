#!/usr/bin/python
from __future__ import absolute_import, division, print_function, unicode_literals
from builtins import (ascii, bytes, chr, dict, filter, hex, input,
                      int, map, next, oct, open, pow, range, round,
                      str, super, zip)

__metaclass__ = type
import logging
import random
import math

log=logging.getLogger(__name__)
class NoopRevertWarning(UserWarning):
    pass

class MCMCSampler:
    '''
    Sample using tradition accept/reject sampling.
    '''
    def __init__(self, sm, energy_function, mover, stats_collector, rerun_prev_energy=False):
        """
        :param sm: A fess.builder.models.SpatialModel instance. The RNA that will be sampled.
        :param energy_function: A fess.builder.energy.CombinedEnergy instance.
                                Used to evaluate the structure during the accept/reject step
        :param mover: A fess.builder.move.Mover instance.
                      It generated the next SpatialModel from the previous.
        """
        self.sm = sm
        self.mover = mover
        self.energy_function = energy_function
        self.stats_collector =  stats_collector
        self.rerun_prev_energy = rerun_prev_energy
        self.last_clashes=[]
        self.last_bad_mls=[]

        #: Store the previous energy.
        log.debug("MCMCSampler __init__ calling eval_energy")
        self.prev_energy = energy_function.eval_energy(sm.bg, sampled_stats=sm.elem_defs)
        log.info("Initial energy of the SpatialModel is {}".format(self.prev_energy))
        log.info("Junction energy is %s", {k:v.shortname for k,v in sm.junction_constraint_energy.items()})
        #: Store the previouse constituing energies (for StatisticsCollector)
        self.prev_constituing = self.energy_function.constituing_energies

        self.energy_function.accept_last_measure()

        #: BT: I do not know, if this is still needed!
        self.sm.get_sampled_bulges() #Store in sm which bulges are sampled (vs broken ml-segments)

        #: Keep track of the number of performed sampling steps.
        self.step_counter = 0

    def eval_energy(self):
        if self.sm.fulfills_constraint_energy():
            if self.sm.constraint_energy is None:
                self.last_clashes="no_energy"
            else:
                self.last_clashes=[]
            if self.sm.junction_constraint_energy is None:
                self.last_bad_mls="no_energy"
            else:
                self.last_bad_mls=[]
            self.last_bad_mls=[]
            return self.energy_function.eval_energy(self.sm.bg, sampled_stats=self.sm.elem_defs)
        else:
            if self.sm.constraint_energy is None:
                self.last_clashes="no_energy"
            else:
                self.last_clashes=self.sm.constraint_energy.bad_bulges
            if self.sm.junction_constraint_energy is None:
                self.last_bad_mls="no_energy"
            else:
                self.last_bad_mls=set( ml for e in self.sm.junction_constraint_energy.values() for ml in e.bad_bulges)

            log.info("Constraint energy not fulfilled!")
            return float("inf")

    def step(self):
        """
        Make a single sampling move and accept or reject the resulting RNA conformation.
        """
        self.step_counter += 1
        if self.rerun_prev_energy:
            # The energy of staying may get worse with every reject step
            self.prev_energy = self.energy_function.eval_energy(self.sm.bg , sampled_stats=self.sm.elem_defs)
        #Make a sinle move (i.e. change the Spatial Model)
        movestring = self.mover.move(self.sm)
        # Accept or reject the new spatial model based on the energy.
        # This stores the new energy as self.prev_energy
        ms, accepted = self.accept_reject()
        movestring += ms
        self.stats_collector.update_statistics( self.sm, self.prev_energy,
                                                self.prev_constituing, movestring,
                                                self.last_clashes,
                                                self.last_bad_mls )
        return accepted

    def accept_reject(self):
        """
        Evaluate the energy of self.sm and either accept or reject the new conformation.
        """
        log.debug("MCMCSampler accept_reject calling eval_energy")
        energy = self.eval_energy()

        movestring=[]
        movestring.append("{:.3f}".format(self.prev_energy))
        movestring.append("->")
        movestring.append("{:.3f};".format(energy))

        if energy <= self.prev_energy:
            movestring.append("A")
            # lower energy means automatic acceptance accordint to the
            # metropolis hastings criterion
            self.accept(energy)
            accepted = True
        else:
            # calculate a probability
            r = random.random()
            if r > math.exp(self.prev_energy - energy):
                movestring.append("R")
                # reject the sampled statistic and replace it with the old one
                self.reject()
                accepted = False
            else:
                movestring.append("A")
                self.accept(energy)
                accepted = True
        return "".join(movestring), accepted

    def accept(self, energy):
        """
        :param energy: The energy of the accepted state.
                       This is to avoid expensive recalculation of the energy
        """
        # accept the new statistic
        self.prev_energy = energy
        self.prev_constituing =  self.energy_function.constituing_energies
        self.energy_function.accept_last_measure()
        for e in self.energy_function.iterate_energies():
            if hasattr(e, "accepted_projDir"):
                self.sm.bg.project_from=e.accepted_projDir
        for e in self.energy_function.iterate_energies():
            if hasattr(e, "used_stat"):
                self.sm.elem_defs[e.element] = e.used_stat
        self.sm.bg.infos["Energy"]=["{} {}".format(energy, self.energy_function.shortname)]

    def reject(self):
        self.energy_function.reject_last_measure()
        try:
            self.mover.revert(self.sm)
        except RuntimeError as e:
            raise
            #This warning will be ignored in ReplicaExchangeSimulations
            warnings.warn(e.message, NoopRevertWarning)
        # We need to recaluculate the prev_energy, because Energy might have been recalibrated.
        log.debug("MCMCSampler After rejecting: reject calling eval_energy again")
        self.prev_energy = self.energy_function.eval_energy(self.sm.bg, sampled_stats=self.sm.elem_defs)
