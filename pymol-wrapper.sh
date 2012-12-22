#!/bin/bash

PDBNAME=
DEBUG=
TEXT=
BASE_NORMALS=
MAX_STEM_DISTANCE=0

while getopts "xtbm:p:" OPTION
do
    case $OPTION in
        x)
            set -x
            DEBUG=1
            ;;
        t)
            TEXT=-t
            ;;
        b)
            BASE_NORMALS=-b
            ;;
        m)
            MAX_STEM_DISTANCE=$OPTARG
            ;;
        p)
            PDBNAME=$OPTARG
            echo $PDBNAME
            ;;
    esac
done

if [[ -z $PDBNAME ]]
then
    echo "Usage.sh: ./pymol-wrapper.sh -p struct_id"
    exit 1
fi

OUTPUT_DIR=~/data/ernwin/processed/${PDBNAME}
OUTPUT_PYMOL_DIR=${OUTPUT_DIR}/pymol
OUTPUT_PREPARE_DIR=${OUTPUT_DIR}/prepare

if [[ -z $DEBUG ]]
then
    ./fess/scripts/create_pymol.sh $TEXT $BASE_NORMALS -m $MAX_STEM_DISTANCE -p fess/structures/$PDBNAME.pdb
else
    ./fess/scripts/create_pymol.sh $TEXT $BASE_NORMALS -m $MAX_STEM_DISTANCE -xp fess/structures/$PDBNAME.pdb
fi

pymol $OUTPUT_PREPARE_DIR/temp.pdb $OUTPUT_PYMOL_DIR/cartoon.pml 
