from builtins import zip
import multiprocessing

def sample_single(sample_queue, exchange_queue)

class ReplicaExchange(object):
    def __init__(self, sampler):
        self.sampler = sampler
        self.sample_queue = multiprocessing.Queue()
        self.exchange_queue = multiprocessing.Queue()
        for i, sampler in enumerate(self.sampler):
            self.sample_queue.append((i,sampler))
        #Spawn processes that listen at the sample queue and push samplers into the exchange queue!
    def run(self, steps=1000):
        for i in range(steps):
            self.run_single()
        # Do the cleanup of the multiple processes!
    def run_single(self):
        for pair_of_samplers in self.exchange_queue: # Problem: I cannot guarantee the order samplers are put into the exchange queue!!! I need to keep track of the sampler ids.
            try_exchange
            self.sample_queue.extend(pair_of_samplers)
