from __future__ import absolute_import, division, print_function, unicode_literals
from builtins import (ascii, bytes, chr, dict, filter, hex, input,
                      int, map, next, oct, open, pow, range, round,
                      str, super, zip)
from future.builtins.disabled import (apply, cmp, coerce, execfile,
                             file, long, raw_input, reduce, reload,
                             unicode, xrange, StandardError)
                      
import itertools
from multiprocessing import Process, Pipe
import fess.builder.models as fbm
import forgi.threedee.model.coarse_grain as ftmc

import logging
log = logging.getLogger(__name__)
import os
import numpy as np
import random
import cProfile

def try_replica_exchange(sampler1, sampler2):
    sm1 = sampler1.sm
    sm2 = sampler2.sm
    sm1e2 = sampler2.energy_function.eval_energy(sm1)
    sm2e1 = sampler1.energy_function.eval_energy(sm2)
    sm1e1 = sampler1.prev_energy
    sm2e2 = sampler2.prev_energy
    total_prev_e = sm1e1 + sm2e2
    total_exchanged_e = sm1e2 + sm2e1
    p = np.exp(total_prev_e - total_exchanged_e) #np.exp can return inf instead of raising an Overflow error
    p = min(1, p)
    r =random.random()
    if r<=p:
        sampler1.sm = sm2
        sampler2.sm = sm1
        sampler1.accept(sm2e1)
        sampler2.accept(sm1e2)
        movestring = "RE {}<->{};A".format(hex(id(sampler1)), hex(id(sampler2)))
    else:
        sampler1.reject()
        sampler2.reject()
        movestring = "RE {}xxx{};R".format(hex(id(sampler1)), hex(id(sampler2)))

    for sampler in [sampler1, sampler2]:
        sampler.stats.update_statistics(sampler.sm, sampler.prev_energy, sampler.prev_constituing, movestring)
    
    
    
class ReplicaExchange(object):
    def __init__(self, sampler_list):
        self.sampler_list = sampler_list
    
    def run(self, steps):
        for i in range(steps):
            for j, sampler in enumerate(self.sampler_list):
                log.info("Sampler {}: changing element".format(j))
                sampler.step()
            for j, sampler1 in enumerate(self.sampler_list):
                if j+1<len(self.sampler_list):
                    log.info("Trying Replica exchange")
                    try_replica_exchange(sampler1, self.sampler_list[j+1])

class MultiPipeConnection(object):
    def __init__(self, names, connections):
        self._connections = { n:c for n,c in zip(names, connections)}
            
    def send(self, name, object):
        log.warning("PID {}: Sending {}".format(os.getpid(),name))
        self._connections[name].send(object)
    def recv(self, name):
        log.warning("PID {}: Waiting to receive {}".format(os.getpid(),name))
        r =  self._connections[name].recv()        
        log.warning("PID {}: Received {}.".format(os.getpid(), name))
        return r

        
def MultiPipe(names):
    conns1, conns2 = [],[]
    for i in range(len(names)):
        c1, c2 = Pipe()
        conns1.append(c1)
        conns2.append(c2)
    return MultiPipeConnection(names, conns1), MultiPipeConnection(names, conns2)
    
class MultiprocessingReProcess(Process):
    def __init__(self, pipe_lower, pipe_higher, sampler, prev_sm, next_sm, steps, idnr):
        """
        :param pipe_lower:  The end of the pipe connecting to the neighboring process with the lower 
                            temperature
        :param pipe_higher: The end of the pipe connecting to the neighboring sampler process with 
                            higher temperature.
        :param prev_sm: We keep a copy of the neighboring process's sm to avoid the overhead of 
                        always sending a the sm between processes via the pipe.
                        This can save overhead because we expect changes to the sm to be 
                        rather infrequent.
        :paramnext_sm: See prev_sm, only for the other neighbor process.
        """
        Process.__init__(self)
        self.pipe_lower=pipe_lower
        self.pipe_higher=pipe_higher
        self.sm_lower = prev_sm
        self.sm_higher = next_sm
        self.sampler = sampler
        self.steps = steps
        self.id = idnr

    def run_exchange(self):
        log.warning("Sampler {} running with pid {}. Is lowest? {}".format(self.id, os.getpid(), self.pipe_lower is None))
        for step in range(self.steps):
            sm_changed = self.sampler.step() #Return a boolean indicating if the step was accepted.
            replica_exchanged = False
            movestring ="RE:"
            if self.pipe_lower is not None:
                energy_curr_with_lower = self.recv_energy_with_different_sampler(sm_changed) 
                energy_lower_with_lower = self.recv_step_result() 
                energy_lower_with_curr = self.sampler.energy_function.eval_energy(self.sm_lower) 
                
                old_energy = self.sampler.prev_energy + energy_lower_with_lower
                total_exchanged_e = energy_curr_with_lower + energy_lower_with_curr
                p = np.exp(old_energy - total_exchanged_e) #np.exp can return inf instead of raising an Overflow error
                r=random.random()
                if r<=p:
                    #EXCHANGE
                    self.sampler.sm, self.sm_lower = self.sm_lower, self.sampler.sm
                    self.sampler.accept(energy_lower_with_curr)
                    self.notify_lower_about_exchange(True)
                    sm_changed = True
                    replica_exchanged = True
                    movestring += "lower<>"
                else:
                    #NO EXCHANGE
                    self.sampler.reject()
                    self.notify_lower_about_exchange(False)
                    movestring += "lowerXX"
            movestring+="this"
            if self.pipe_higher is not None:
                self.send_step_result(sm_changed)
                energy_higher_with_curr = self.send_energy_of_different_sm() #Stores sm_higher
                if self.recv_exchange_notification():
                    self.sampler.sm, self.sm_higher = self.sm_higher, self.sampler.sm
                    self.sampler.accept(energy_higher_with_curr)
                    replica_exchanged = True
                    movestring+="<>higher"
                else:
                    self.sampler.reject()
                    movestring+="XXhigher"
            if replica_exchanged:
                movestring+=";A"
            else:
                movestring+=";R"
            log.error("PID {}: {}, {}".format(os.getpid(), self.sampler.prev_energy, type(self.sampler.prev_energy)))
            self.sampler.stats.update_statistics(self.sampler.sm, self.sampler.prev_energy, self.sampler.prev_constituing, movestring)

    def run(self):
        cProfile.runctx('self.run_exchange()', globals(), locals(), 'prof%d.prof' %self.id)
                
    def notify_lower_about_exchange(self, exchange_happened=True):
        self.pipe_lower.send("exchange_result", exchange_happened)
    
    def recv_exchange_notification(self):
        return self.pipe_higher.recv("exchange_result")
    
    @staticmethod
    def sm_from_cg_string(cg_string):
        bg = ftmc.CoarseGrainRNA()
        bg.from_cg_string(cg_string)
        sm = fbm.SpatialModel(bg)
        sm.load_sampled_elems()
        sm.new_traverse_and_build()
        return sm

    def recv_step_result(self):
        cg_string, energy = self.pipe_lower.recv("step_result")
        if cg_string is not None:
            self.sm_lower = self.sm_from_cg_string(cg_string)
        log.warning("PID {} recv_step_result: Energy {} received".format(os.getpid(), energy))
        return energy
    
    def send_step_result(self, changed):
        if not changed:
            self.pipe_higher.send("step_result", (None, self.sampler.prev_energy))
        else:
            self.pipe_higher.send("step_result", (self.sampler.sm.bg.to_cg_string(), self.sampler.prev_energy))
        
    
    def recv_energy_with_different_sampler(self, changed = True):
        #Inform lower neighbor process about the changes in the sm
        if changed:
            self.pipe_lower.send("exchanged_energy", self.sampler.sm.bg.to_cg_string())
        else:
            self.pipe_lower.send("exchanged_energy", None)            
        #Receive energy of our sm with other sampler Temperature
        energy = self.pipe_lower.recv("exchanged_energy") 
        log.warning("PID {} recv_energy_with_different_sampler: Energy {} received".format(os.getpid(), energy))
        return energy

    def send_energy_of_different_sm(self):
        #Receive changed sm from higher temperature neighbor
        cg_string = self.pipe_higher.recv("exchanged_energy") # RECV <-A
        if cg_string is not None:
            self.sm_higher = self.sm_from_cg_string(cg_string)
        #Calculate energy and send it back
        e = self.sampler.energy_function.eval_energy(self.sm_higher)
        self.pipe_higher.send("exchanged_energy", e) 
        return e
    
def start_parallel_replica_exchange(sampler_list, num_steps):
    prev_conn = None
    processes = []
    for i, sampler in enumerate(sampler_list):
        if i<len(sampler_list)-1:
            conn1, conn2 = MultiPipe(["step_result", "exchanged_energy", "exchange_result"])
        else:
            conn1=None
        if i>0:
            prev_sm=sampler_list[i-1].sm
        else:
            prev_sm = None
        if i+1<len(sampler_list):
            next_sm = sampler_list[i+1].sm
        else:
            next_sm = None
        processes.append(MultiprocessingReProcess(prev_conn, conn1, sampler, prev_sm, next_sm, num_steps, i))
        prev_conn = conn2
    for p in processes:
        p.start()
    for p in processes:
        p.join() #Block until all subprocesses are done.