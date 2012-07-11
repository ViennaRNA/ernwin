#!/bin/bash

base=${1##*/}
pdb=${base%\.*}

OUTPUT_DIR=output/${pdb}

OUTPUT_PREPARE_DIR=${OUTPUT_DIR}/prepare
OUTPUT_GRAPH_DIR=${OUTPUT_DIR}/graph
OUTPUT_REPORT_DIR=${OUTPUT_DIR}/report
OUTPUT_NEATO_DIR=${OUTPUT_DIR}/neato
OUTPUT_STATS_DIR=${OUTPUT_DIR}/stats

#./scripts/graph_to_angles.py $OUTPUT_GRAPH_DIR/temp.comp > $OUTPUT_STATS_DIR/temp.angles
#./scripts/graph_to_distances.py $OUTPUT_GRAPH_DIR/temp.comp > $OUTPUT_STATS_DIR/temp.gdistances
#./scripts/all_long_range_distances.py $OUTPUT_GRAPH_DIR/temp.comp > $OUTPUT_STATS_DIR/temp.longrange.all
#./scripts/long_range_distances.py $OUTPUT_GRAPH_DIR/temp.comp > $OUTPUT_STATS_DIR/temp.longrange.contact
./scripts/long_range_interaction_eval.py $OUTPUT_GRAPH_DIR/temp.comp > $OUTPUT_STATS_DIR/temp.energy
#./scripts/bobbins_energy.py $OUTPUT_GRAPH_DIR/temp.comp > $OUTPUT_STATS_DIR/energy.long_range_distance

