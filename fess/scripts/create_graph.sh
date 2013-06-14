#!/bin/bash

set -x 

base=${1##*/}
pdb=${base%\.*}

OUTPUT_DIR=~/data/ernwin/processed/${pdb}

OUTPUT_PREPARE_DIR=${OUTPUT_DIR}/prepare
OUTPUT_GRAPH_DIR=${OUTPUT_DIR}/graph
OUTPUT_REPORT_DIR=${OUTPUT_DIR}/report
OUTPUT_NEATO_DIR=${OUTPUT_DIR}/neato

LOCAL_SCRIPT_DIR=fess/scripts

./$LOCAL_SCRIPT_DIR/create_bulge_graph.py $OUTPUT_PREPARE_DIR/temp.dotplot $pdb >  $OUTPUT_GRAPH_DIR/temp.bulge
./$LOCAL_SCRIPT_DIR/check_bulge_integrity.py $OUTPUT_GRAPH_DIR/temp.bulge

./$LOCAL_SCRIPT_DIR/graph_to_neato.py $OUTPUT_GRAPH_DIR/temp.bulge > $OUTPUT_NEATO_DIR/temp_bulge.neato
#neato -Tpng $OUTPUT_NEATO_DIR/temp_bulge.neato -o $OUTPUT_REPORT_DIR/temp_bulge.png

./$LOCAL_SCRIPT_DIR/collapse_graph.py $OUTPUT_GRAPH_DIR/temp.bulge > $OUTPUT_GRAPH_DIR/temp.graph
./$LOCAL_SCRIPT_DIR/check_bulge_integrity.py $OUTPUT_GRAPH_DIR/temp.graph

./$LOCAL_SCRIPT_DIR/graph_to_neato.py $OUTPUT_GRAPH_DIR/temp.graph > $OUTPUT_NEATO_DIR/temp_graph.neato
./$LOCAL_SCRIPT_DIR/add_sequence.py $OUTPUT_GRAPH_DIR/temp.graph $OUTPUT_PREPARE_DIR/temp.seq > $OUTPUT_GRAPH_DIR/temp.graph.seq
./$LOCAL_SCRIPT_DIR/add_3d_information.py $OUTPUT_GRAPH_DIR/temp.graph.seq $OUTPUT_PREPARE_DIR/temp.pdb > $OUTPUT_GRAPH_DIR/temp.graph.coord
./$LOCAL_SCRIPT_DIR/get_long_range_interactions.py $OUTPUT_GRAPH_DIR/temp.graph.coord $OUTPUT_PREPARE_DIR/temp.mcannotate > $OUTPUT_GRAPH_DIR/temp.comp

./$LOCAL_SCRIPT_DIR/graph_to_neato.py $OUTPUT_GRAPH_DIR/temp.comp > $OUTPUT_NEATO_DIR/temp_comp.neato
#neato -Tpng $OUTPUT_NEATO_DIR/temp_comp.neato -o $OUTPUT_REPORT_DIR/temp_comp.png
