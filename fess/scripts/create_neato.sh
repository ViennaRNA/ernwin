#!/bin/bash

base=${1##*/}
pdb=${base%\.*}

OUTPUT_DIR=~/data/ernwin/processed/${pdb}
SCRIPTS_DIR=fess/scripts

OUTPUT_PREPARE_DIR=${OUTPUT_DIR}/prepare
OUTPUT_GRAPH_DIR=${OUTPUT_DIR}/graph
OUTPUT_REPORT_DIR=${OUTPUT_DIR}/report
OUTPUT_NEATO_DIR=${OUTPUT_DIR}/neato



$SCRIPTS_DIR/graph_to_neato.py $OUTPUT_GRAPH_DIR/temp.bulge > $OUTPUT_NEATO_DIR/temp_bulge.neato
neato -Tpng $OUTPUT_NEATO_DIR/temp_bulge.neato -o $OUTPUT_REPORT_DIR/temp_bulge.png

$SCRIPTS_DIR/graph_to_neato.py $OUTPUT_GRAPH_DIR/temp.graph > $OUTPUT_NEATO_DIR/temp_graph.neato
neato -Tpng $OUTPUT_NEATO_DIR/temp_graph.neato -o $OUTPUT_REPORT_DIR/temp_graph.png

$SCRIPTS_DIR/graph_to_neato.py $OUTPUT_GRAPH_DIR/temp.comp > $OUTPUT_NEATO_DIR/temp_comp.neato
neato -Tpng $OUTPUT_NEATO_DIR/temp_comp.neato -o $OUTPUT_REPORT_DIR/temp_comp.png
