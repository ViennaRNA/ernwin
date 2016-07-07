from __future__ import print_function, absolute_import, division
from builtins import (ascii, bytes, chr, dict, filter, hex, input,
                      map, next, oct, open, pow, range, round,
                      str, super, zip) 


import argparse, random, math
import forgi.utilities.stuff as fus
import forgi.threedee.model.coarse_grain as ftmc
import forgi.threedee.model.stats as ftms
import forgi.threedee.utilities.rmsd as ftur
import forgi.threedee.utilities.vector as ftuv
import forgi.threedee.utilities.graph_pdb as ftug
import numpy as np
import scipy.cluster.hierarchy as hcluster
import scipy.stats
import matplotlib.pyplot as plt
import copy, sys
import collections as col
import itertools
from collections import defaultdict

def generateParser():
    parser=argparse.ArgumentParser( description="Report coaxial stacking.")
    parser.add_argument("files", type=str, nargs="+", help="One or more cg files that all have the same bulge graph!")
    parser.add_argument("-v", "--verbose", action="store_true", help="Verbose output")
    parser.add_argument("-p", "--plot", action="store_true", help="Plot occurrence of coaxial stacking over time.")
    parser.add_argument("-m", "--per-multiloop", action="store_true", help="Print statistics per multiloop.")
    parser.add_argument("-r", "--per-residue", action="store_true", help="Report stacking residues.")
    return parser

def is_close(coords1, coords2, cutoff=14, verbose=False):
    """
    See  doi:  10.1261/rna.305307, PMCID: PMC1894924 for cuttoffs for stacking
    """
    returnvalue=0
    min_dist=float("inf")
    if ftuv.magnitude(coords1[0]-coords2[0])<min_dist:
        returnvalue=1
        min_dist=ftuv.magnitude(coords1[0]-coords2[0])
    if ftuv.magnitude(coords1[0]-coords2[1])<min_dist:
        returnvalue=2
        min_dist=ftuv.magnitude(coords1[0]-coords2[1])
    if ftuv.magnitude(coords1[1]-coords2[0])<min_dist:
        returnvalue=3
        min_dist=ftuv.magnitude(coords1[1]-coords2[0])
    if ftuv.magnitude(coords1[1]-coords2[1])<min_dist:
        returnvalue=4
        min_dist=ftuv.magnitude(coords1[1]-coords2[1])
    if verbose: print("Distance {}".format(min_dist), file=sys.stderr)
    if min_dist<cutoff:
        lsd=ftuv.line_segment_distance(coords1[0], coords1[1], coords2[0], coords2[1])
        lsd=ftuv.magnitude(lsd[1]-lsd[0])
        if lsd<min_dist:
            if verbose: print("line segment distance {}<{}. Stems are overlapping".format(lsd, min_dist), file=sys.stderr)
            return False, min_dist
        return returnvalue, min_dist
    return False, min_dist

def is_coaxial(coords1, coords2, stack_type, cutoff=math.acos(0.75), verbose=False):
    c1=coords1[(stack_type-1)//2]
    c2=coords2[(stack_type-1)%2]
    f1=coords1[int(not (stack_type-1)//2)]
    f2=coords2[int(not (stack_type-1)%2)]
    v1=c1-f1
    v2=f2-c2
    #v1=coords1[1]-coords1[0]
    #v2=coords2[1]-coords2[0]
    angle = ftuv.vec_angle(v1, v2)
    if verbose: print("Angle is {} ({}), cutoff: {} ".format(angle, math.degrees(angle), cutoff))
    if angle>cutoff:
        return False, angle
    return True, angle

def shear_ok(coords1, coords2, stack_type, cutoff_a=1.05, cutoff_o=10, verbose=False):
  v1=coords1[1]-coords1[0]
  v2=coords2[1]-coords2[0]

  c1=coords1[(stack_type-1)//2]
  c2=coords2[(stack_type-1)%2]
  diff=c2-c1
  angle=ftuv.vec_angle(v1, diff)
  if angle>math.pi/2:
        angle=math.pi-angle
  if angle>cutoff_a:
      if verbose: print("Shear angle is {} > {} ".format(abs(angle), cutoff_a))
      return False 
  s_dist = point_line_distance(c2, c1, v1)
  if s_dist>cutoff_o:
      if verbose: print("Shear distance is {} > {} ".format(s_dist, cutoff_o))
      return False
  if verbose: print("Shear ok: {}, {} ".format(abs(angle), s_dist))
  return True

def point_line_distance(p, m, v):
    """
    :param p: the point
    :param m: a point on the line
    :param v: the direction of the line
    Formula from http://onlinemschool.com/math/library/analytic_geometry/p_line/
    """
    return ftuv.magnitude(np.cross((m-p), v))/ftuv.magnitude(v)
parser=generateParser()

def is_stacking(stem1, stem2, cg, verbose=False):
    stack_type, distance = is_close(cg.coords[stem1], cg.coords[stem2], verbose=verbose)
    is_coax, ang= is_coaxial(cg.coords[stem1], cg.coords[stem2], stack_type, verbose=verbose)
    if stack_type and is_coax:
        if shear_ok(cg.coords[stem1], cg.coords[stem2], stack_type, verbose=verbose):
            return True
    return False

def report_all_stacks(cg, verbose=False):
    not_close=0
    not_coaxial=0
    sheared=0
    stacks=col.defaultdict(list)
    for stem1, stem2 in itertools.combinations(cg.stem_iterator(), 2):
        if verbose: print("===== {}, {} =====".format(stem1, stem2), file=sys.stderr)
        stack_type, _ = is_close(cg.coords[stem1], cg.coords[stem2], verbose=verbose)
        if stack_type:
            is_coax, ang = is_coaxial(cg.coords[stem1], cg.coords[stem2], stack_type, verbose=verbose)
            if is_coax:
                if shear_ok(cg.coords[stem1], cg.coords[stem2], stack_type, verbose=verbose):
                    stacks[stem1].append(stem2)
                    stacks[stem2].append(stem1)
                else: sheared+=1
            else: not_coaxial+=1
        else: not_close+=1
    if verbose:
        print( "{} stems: {} combinations not close, {} not coaxial, "
               "{} sheared.".format(len(x for x in cg.defines if x[0]=="s"), 
                                    not_close, not_coaxial, sheared))
    return stacks

def report_stacks_per_multiloop(cgs):
    stacks=defaultdict(int)
    no_stacks=defaultdict(int)
    ml_types=defaultdict(int)
    for cg in cgs:
        loops=[]
        for d in cg.defines:
            if d[0]=="m":
                loops.append(tuple(cg.find_bulge_loop(d, 200)))
        loops=set(loops)
        for loop in loops:        
            stack_type=[]
            defines = [l for l in loop if l[0]=="m" ]
            degree = len(defines)
            lengths = [ (cg.get_length(l)) for l in defines  ]
            connections = [ cg.connection_type(l, cg.connections(l)) for l in defines ]
            for i, l in enumerate(defines):
                stem1, stem2 = cg.connections(l)                    
                key=(degree, i, tuple(lengths), tuple(connections))
                if cg.is_stacking(l):
                    stack_type.append(connections[i])
                    stacks[key]+=1
                else:
                    no_stacks[key]+=1
            stack_type.sort()
            key=tuple([degree, tuple(lengths), tuple(stack_type)])
            ml_types[key]+=1
    return stacks, no_stacks, ml_types

def print_collapsed_keys(stack, no_stacks, collapse, name):
        stacks_collapsed = defaultdict(int)
        no_stacks_collapsed = defaultdict(int)
        for key in stacks:
            newkey = collapse(key)
            stacks_collapsed[newkey]+=stacks[key]
        for key in no_stacks:
            newkey = collapse(key)
            no_stacks_collapsed[newkey]+=no_stacks[key]
        print(name)
        allkeys = set(stacks_collapsed) | set(no_stacks_collapsed)
        for key in sorted(allkeys):
            print("  {} {}: {} occurrences, {} stacks, {} not stacking, {}%".format( 
                          name, key, 
                          stacks_collapsed[key]+no_stacks_collapsed[key],
                          stacks_collapsed[key],
                          no_stacks_collapsed[key],
                          stacks_collapsed[key]/(stacks_collapsed[key]+no_stacks_collapsed[key])*100))

def print_signature(stacks):
    if not stacks:
        return "no stacks"
    sig=[]
    for stack in stacks:
        sig.append("{}+{}".format(stack[0], stack[1]))
    return "; ".join(sig)

class color(object):
  def __getitem__(self, i):
      return ((31*i)%50/100.+0.5,(23*i)%50/100.+0.5, (29*i)%50/100.+0.5)
  def __len__(self):
      return 10**5
if __name__=="__main__":
    args = parser.parse_args()
    if args.per_multiloop:
        cg_files=[]
        for filename in args.files:
            cg_files.append( ftmc.CoarseGrainRNA(filename) )
        stacks, no_stacks, ml_types = report_stacks_per_multiloop(cg_files)
        name="Connection Type"
        collapse=lambda key: key[3][key[1]]
        print_collapsed_keys(stacks, no_stacks, collapse, name)
        name="Loop Degree"
        collapse=lambda key: key[0]
        print_collapsed_keys(stacks, no_stacks, collapse, name)
        name="Loop Length"
        collapse=lambda key: key[2][key[1]]
        print_collapsed_keys(stacks, no_stacks, collapse, name)
        name="Loop Degree, Connection Type"
        collapse=lambda key: (">5", key[3][key[1]]) if key[0]>5 else (key[0], key[3][key[1]])
        print_collapsed_keys(stacks, no_stacks, collapse, name)
        name="Loop Length, Connection Type"
        collapse=lambda key: (">2", key[3][key[1]]) if key[2][key[1]]>2 else (key[2][key[1]], key[3][key[1]])
        print_collapsed_keys(stacks, no_stacks, collapse, name)
        print("3-way junction classification")
        for mlt in sorted(ml_types):
            if mlt[0]!=3: continue
            print("   linker lengths: {}, stacking {}: {} ({}%)".format(mlt[1], mlt[2],ml_types[mlt],(100* ml_types[mlt]/sum(ml_types[k] for k in ml_types if k[0]==mlt[0] and k[1]==mlt[1]))))
        m3 = col.defaultdict(int)
        for mlt in ml_types:
            if mlt[0]!=3: continue
            key=tuple([abs(x) for x in mlt[2]])
            m3[key]+=ml_types[mlt]
        print("Without linker-length classification:")
        for mlt in sorted(m3):
            print("   stacking {}: {} ({}%)".format(mlt ,m3[mlt],(100* m3[mlt]/sum(m3.values()))))
    elif args.per_residue:
        for filename in args.files:
            cg = ftmc.CoarseGrainRNA(filename)
            print(cg.name)
            for d in cg.defines:
                if d[0] in "mi":
                    if cg.is_stacking(d):
                        stem1, stem2 = cg.connections(d)
                        side_nts = cg.get_connected_residues(stem1, stem2, d)
                        for nts in side_nts:
                            print ("{} and {} are stacking".format(cg.seq_ids[nts[0]-1], cg.seq_ids[nts[1]-1]))
    else:
        stacks=[]
        not_close=0
        not_coaxial=0
        sheared=0
        angles=col.defaultdict(list)
        distances=col.defaultdict(list)
        structures=col.Counter()
        for filename in args.files:
            stacks.append(col.defaultdict(int))
            if args.verbose: print("File {}".format(filename), file=sys.stderr)
            cg = ftmc.CoarseGrainRNA(filename)
            for stem1, stem2 in itertools.combinations(cg.stem_iterator(), 2):
                if args.verbose: print("===== {}, {} =====".format(stem1, stem2), file=sys.stderr)
                stack_type, distance = is_close(cg.coords[stem1], cg.coords[stem2], verbose=args.verbose)
                if stem1<stem2:
                    distances[(stem1,stem2)].append(distance)
                else:
                    distances[(stem2,stem1)].append(distance)
                is_coax, ang= is_coaxial(cg.coords[stem1], cg.coords[stem2], stack_type, verbose=args.verbose)
                ang=math.degrees(ang)
                if stem1<stem2:
                    angles[(stem1,stem2)].append(ang)
                else:
                    angles[(stem2,stem1)].append(ang)
                if stack_type:
                    if is_coax:
                        if shear_ok(cg.coords[stem1], cg.coords[stem2], stack_type, verbose=args.verbose):
                            if stem1<stem2:
                                stacks[-1][(stem1,stem2)]+=1
                            else:
                                stacks[-1][(stem2,stem1)]+=1
                        else: sheared+=1
                    else: not_coaxial+=1
                else: not_close+=1
            stack_signature = tuple(sorted(stacks[-1].keys()))
            structures[stack_signature]+=1
        newstacks=col.defaultdict(int)
        for ddict in stacks:
            for k,v in ddict.items():
                newstacks[k]+=v
        if args.plot:
            to_plot=col.defaultdict(list)
            for i, ddict in enumerate(stacks):
                for key in  [("s0","s1"), ("s1","s2"),("s2","s3"), ("s0","s3"), ("s0","s2"), ("s1","s3")]:
                    if key in ddict:
                        to_plot[key].append(ddict[key])
                    else:
                        to_plot[key].append(0)            
            fig, ax = plt.subplots(3)
            a=-0.2
            all_angles=[x for k in angles for x in angles[k]]
            all_distances=[x for k in distances for x in distances[k]]
            xs_ang=np.linspace(min(all_angles), max(all_angles), 500)
            xs_dist=np.linspace(min(all_distances), max(all_distances), 500)
            for key in [("s0","s1"), ("s1","s2"),("s2","s3"), ("s0","s3"), ("s0","s2"), ("s1","s3")]:
                if key in [("s0","s1"), ("s1","s2"),("s2","s3"), ("s0","s3")]:
                    width=2
                else:
                    width=1
                if key==("s0","s3"):
                    style="dashed"
                else:
                    style="solid"
                a+=0.05
                ax[0].plot(list(range(len(to_plot[key]))), np.array(to_plot[key])+a, "o", markeredgewidth=0.0, label="{}".format(key))
                d=scipy.stats.gaussian_kde(angles[key])
                ax[1].plot(xs_ang, d(xs_ang), linewidth=width, linestyle=style)
                d=scipy.stats.gaussian_kde(distances[key])
                ax[2].plot(xs_dist, d(xs_dist),linewidth=width, linestyle=style)
            ax[0].set_ylim([-0.5,1.5])
            ax[0].legend(loc="lower right")
            ax[0].set_xlabel("simulation step")
            ax[0].set_ylabel("~1: present; ~0: absent")
            ax[1].set_xlabel("angle")
            ax[1].set_ylabel("frequency")
            ax[2].set_xlabel("distance")
            ax[2].set_ylabel("frequency")
            plt.show()
            
        print ("{} not close, {} not coaxial, {} sheared".format(not_close, not_coaxial, sheared))
        for k, v in sorted(newstacks.items()):
            print("{}\t\t{}".format(k,v))
        print("The following stacking pattetrs were observed:")
        for k, v in sorted(structures.items()):
            print("{:20}\t{}".format(print_signature(k),v))
        if args.plot:
            colors=color()
            fig,ax = plt.subplots(2,2)
            stru_values = sorted(structures.values(), key=lambda x: x)
            stru_keys   = sorted(structures.keys(), key=lambda x: structures[x])
            ax[0,0].pie(newstacks.values(), labels=list(map(lambda x: "{}+{}".format(*x), newstacks.keys())), autopct="%1.1f", colors=colors)
            ax[0,1].pie(stru_values, labels=list(map(print_signature, stru_keys)), autopct="%1.1f", colors=colors)
            ax[1,1].pie(stru_values[:-1], labels=list(map(print_signature, stru_keys[:-1])), autopct="%1.1f", colors=colors)
            plt.show()
