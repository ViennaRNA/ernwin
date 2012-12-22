#!/bin/bash

PDBNAME=
TEXT=
BASE_NORMALS=
MAX_STEM_DISTANCE=0

while getopts "xtbm:p:" OPTION
do
    case $OPTION in
        x)
            set -x
            ;;
        t)
            TEXT=-x
            ;;
        b)
            BASENORMALS=1
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

base=${PDBNAME##*/}
pdb=${base%\.*}

OUTPUT_DIR=~/data/ernwin/processed/${pdb}

echo output_dir $OUTPUT_DIR
echo $PDBNAME

OUTPUT_GRAPH_DIR=${OUTPUT_DIR}/graph
OUTPUT_PYMOL_DIR=${OUTPUT_DIR}/pymol
OUTPUT_PREPARE_DIR=${OUTPUT_DIR}/prepare

LOCAL_SCRIPT_DIR=fess/scripts

./$LOCAL_SCRIPT_DIR/coordinates_to_pymol.py -m $MAX_STEM_DISTANCE $TEXT $OUTPUT_GRAPH_DIR/temp.comp > $OUTPUT_PYMOL_DIR/coarse_grain.pym 
./$LOCAL_SCRIPT_DIR/graph_to_pymol.py $OUTPUT_GRAPH_DIR/temp.comp $OUTPUT_PYMOL_DIR/coarse_grain.pym > $OUTPUT_PYMOL_DIR/cartoon.pml

if [[ ! -z $BASENORMALS ]]
then
    ./$LOCAL_SCRIPT_DIR/base_normals_to_pymol.py $OUTPUT_PREPARE_DIR/temp.pdb >> $OUTPUT_PYMOL_DIR/coarse_grain.pym

fi
