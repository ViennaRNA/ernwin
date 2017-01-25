from builtins import zip
from threading import Thread
from Queue import Queue, Empty
import random
import numpy as np
import logging
log = logging.getLogger(__name__)
# Multiprocess Queue tutorial: http://toastdriven.com/blog/2008/nov/11/brief-introduction-multiprocessing/

class NoThreadsRunning(RuntimeError):
    """An exception that will be called if all sub-threads have died."""
    pass

class ReplicaExchange(object):
    def __init__(self, sampler, num_threads=1):
        self._replica_count = len(sampler)
        self._sample_queue = Queue()
        self._exchange_queue = Queue()
        for i, sampler in enumerate(sampler):
            self._sample_queue.put((i,sampler))
        #Spawn processes that listen at the sample queue and push samplers into the exchange queue!
        self._processes = []
        for i in range(num_threads):
            p = Thread(target=self._worker, args=(self._sample_queue, self._exchange_queue))
            self._processes.append(p)
        
    @staticmethod
    def _worker(sample_queue, exchange_queue):
        while True:
            i, sampler = sample_queue.get()
            if sampler is None:
                return
            sampler.step() #Change self.spatial model
            exchange_queue.put((i, sampler))
            
    def run(self, steps=1000):
        try:
            for p in self._processes:
                p.start()
            for i in range(steps):
                self._finish_iteration()
        finally:
            # Do the cleanup of the multiple processes!
            for i in range(len(self._processes)):
                self._sample_queue.put((-1, None)) #Poison pill to end all workers.
            for p in self._processes:
                p.join()
        
    def _finish_iteration(self):
        done = set()
        to_process = {}
        while len(done) != self._replica_count:
            while True:
                try:
                    i, sampler = self._exchange_queue.get(timeout = 5)
                    break
                except Empty:
                    log.debug("ExchangeQueue is empty")
                    for t in self._processes:
                        if t.isAlive():
                            log.debug(t," is alive")
                            break
                    else:
                        log.debug("No Thread is alive!")
                        raise NoThreadsRunning("All worker threads have died! Exiting!")
            if i in done: # Result is for the next step. 
                self._exchange_queue.put((i,sampler))
            else:
                assert i not in to_process
                to_process[i] = sampler
                self._exchange_where_possible(to_process, done)
        assert not to_process
        
    def _exchange_where_possible(self, to_process, done):
        for i, sampler1 in list(to_process.items()):
            if (i-1 in to_process or i-1 in done or i==0) and (i+1 in to_process or i+1 in done or i==self._replica_count-1):
                if i-1 in to_process:
                    self._try_exchange(sampler1, to_process[i-1])
                if i+1 in to_process:
                    self._try_exchange(sampler1, to_process[i+1])
                del to_process[i]
                done.add(i)
                self._sample_queue.put((i, sampler1))

    @staticmethod
    def _try_exchange(sampler1, sampler2):
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
        else:
            sampler1.reject()
            sampler2.reject()
