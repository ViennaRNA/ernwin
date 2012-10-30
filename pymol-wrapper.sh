#!/bin/bash

if [ $# -ne 1 ]; then
    echo "Usage.sh: ./pymol-wrapper.sh struct_id"
    exit
fi

./fess/scripts/create_pymol.sh fess/structures/$1.pdb
pymol fess/output/$1/prepare/temp.pdb fess/output/$1/pymol/cartoon.pml
