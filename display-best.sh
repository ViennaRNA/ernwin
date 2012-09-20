#!/bin/bash

while getopts "xp" OPTION
do
    case $OPTION in
        x)
            set -x
            ;;
        p)
            options='out.pdb'
    esac
done

n=0; ./fess/scripts/coordinates_to_pymol.py -x bobbins/best/best${n}.coord > ss.pym; pymol $options temp.pml
