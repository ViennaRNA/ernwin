#!/usr/bin/python

from corgy.utilities.vector import vec_distance

class LongRangeDistanceEnergy:
    def __init__(self):
        pass

    def eval_energy(self, bg):
        keys = bg.defines.keys()
        energy = 0.

        for i in range(len(keys)):
            for j in range(i+1, len(keys)):
    
                d1 = keys[i]
                d2 = keys[j]

                point1 = bg.get_point(d1)
                point2 = bg.get_point(d2)

                dist = vec_distance(point1, point2)

                if dist > 6.0 and dist < 25.0:
                    energy += 1

        return -energy

    def naive_energy(self, bg):
        keys = bg.defines.keys()
        energy = 0.

        for i in range(len(keys)):
            for j in range(i+1, len(keys)):
                energy += 1

        return -energy
