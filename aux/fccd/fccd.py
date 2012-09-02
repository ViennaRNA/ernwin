# Python stdlib
from math import sqrt, pi
from random import randint

# Numpy stuff
#from Numeric import *
#from MLab import eye
from numpy.matlib import eye, zeros
from numpy.random import random
from numpy import array, dot
from numpy.linalg import det, svd

from corgy.utilities.vector import vec_angle
#from LinearAlgebra import singular_value_decomposition, determinant
#from RandomArray import random, normal

# Biopython stuff
from Bio.PDB import *

import sys

S=eye(3,3)
S[2,2]=-1


def fccd(moving, fixed, threshold=0.1, maxit=10000):
    if len(moving)<6:
        raise "Moving too short"
    if len(fixed)!=3:
        raise "Fixed should have length 3"
    lng=len(moving)
    # Coordinates along COLUMNS
    moving_coords=zeros((3, 3))
    fixed_coords=zeros((3, 3))
    it=0
    while 1:
        if it==maxit:
            # Stop - max iterations reached
            return "MAXIT", rmsd, it
        for i in range(2, lng-2):
            it=it+1
            center=moving[i]
            # move to pivot origin
            for j in range(0, 3):
                index=-(3-j)
                v=moving[index]-center
                moving_coords[:,j]=array([v.get_array()]).transpose()
            for j in range(0, 3):
                v=fixed[j]-center
                fixed_coords[:,j]=array([v.get_array()]).transpose()
            # Do SVD
            #a=matrixmultiply(fixed_coords, transpose(moving_coords))
            a = dot(fixed_coords, moving_coords.transpose())
            u, d, vt=svd(a)
            # Check reflection
            if (det(u)*det(vt))<0:
                u=dot(u, S)
            # Calculate rotation
            rot_matrix=dot(u, vt)
            assert(det(rot_matrix)>0)
            # Apply rotation
            for j in range(i+1, len(moving)):
                v=moving[j]-center
                v = dot(rot_matrix, v.get_array())
                v = Vector(array(v)[0])
                v=v+center
                moving[j]=v
            # Calculate RMSD
            rmsd=0.0
            for j in range(1, 4):
                m=moving[-j]
                f=fixed[-j]
                n=(m-f).norm()
                rmsd+=n*n
            rmsd=sqrt(rmsd/3.0)
            if rmsd<threshold:
                # stop - RMSD threshold reached
                return "SUCCESS", rmsd, it


if __name__=="__main__":

    def make_random_chain(n=12):
        """
        Return a list of random vectors, each with distance
        3.8 from the previous vector.
        """
        v=Vector(0,0,0)
        l=[v]
        for i in range(0, n-1):
            nv=Vector(random(3))
            nv.normalize()
            nv=v+nv**3.8
            l.append(nv)
            v=nv
        return l

    def rotate_last_three(chain):
        """
        Take the last three vectors of chain, copy them,
        and apply a random rotation and translation.
        """
        l=[]
        rotation_axis=Vector(random(3))
        angle=0.1+0.1*pi*random()
        m=rotaxis(angle, rotation_axis)
        t=Vector(-0.5*random(3))
        for i in range(-3, 0):
            v=chain[i]-t
            v=v.left_multiply(m)
            l.append(v)
        return l

    # Moving segment
    moving=make_random_chain(20)
    # Fixed segment 
    # Last three residues of the moving segment
    # after applying a random rotation/translation
    fixed=rotate_last_three(moving)

    # Close the gap
    print "moving before:", moving[:4]
    print "fixed before:", fixed

    angles = []
    for i in range(2, len(moving)):
        angles += [vec_angle(moving[i-1] - moving[i-2], moving[i] - moving[i-1])]

    status, rmsd, it=fccd(moving, fixed, 0.01)

    angles1 = []
    for i in range(2, len(moving)):
        angles1 += [vec_angle(moving[i-1] - moving[i-2], moving[i] - moving[i-1])]

    print "moving after:", moving[:4]
    print "fixed after:", fixed

    for i in range(len(angles)):
        print "angles[i] - angles1[i]:", angles[i] - angles1[i]

    #print "fixed:", fixed

    # Print result
    print "Final RMSD ", rmsd
    print "Number of iterations ", it
    print "Status ", status


            




    

    



