#!/bin/bash

if [ $# -ne 1 ]; then
    echo "Usage.sh: ./pymol-wrapper.sh struct_id"
    exit
fi


./scripts/create_pymol.sh structures/$1.pdb
pymol output/$1/prepare/temp.pdb output/$1/pymol/cartoon.pml
