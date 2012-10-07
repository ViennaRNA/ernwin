import math

def clock_angle(a1, a2):
    '''
    The amount one needs to rotate from angle a1
    to angle a2 in a clockwise direction. I.e.
    increasing a1 until it reaches a2.

    @param a1: The first angle.
    @param a2: The second angle.
    '''
    if a2 >= a1:
        return a2 - a1

    else:
        a2 += 2. * math.pi

    return a2 - a1

