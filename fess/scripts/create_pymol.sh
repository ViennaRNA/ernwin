#!/bin/bash

base=${1##*/}
pdb=${base%\.*}

OUTPUT_DIR=fess/output/${pdb}

OUTPUT_GRAPH_DIR=${OUTPUT_DIR}/graph
OUTPUT_PYMOL_DIR=${OUTPUT_DIR}/pymol

LOCAL_SCRIPT_DIR=fess/scripts

./$LOCAL_SCRIPT_DIR/coordinates_to_pymol.py -x $OUTPUT_GRAPH_DIR/temp.comp > $OUTPUT_PYMOL_DIR/coarse_grain.pym 
./$LOCAL_SCRIPT_DIR/graph_to_pymol.py $OUTPUT_GRAPH_DIR/temp.comp $OUTPUT_PYMOL_DIR/coarse_grain.pym > $OUTPUT_PYMOL_DIR/cartoon.pml

