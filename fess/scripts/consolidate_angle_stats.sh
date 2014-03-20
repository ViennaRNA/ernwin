#!/bin/bash
OUTPUT_DIR=~/data/ernwin/processed/

#find output/ -name "temp.angles" | xargs cat > fess/stats/angles.csv
find $OUTPUT_DIR -name "temp.angles" | xargs cat > fess/stats/temp.stats
find $OUTPUT_DIR -name "temp.energy" | xargs cat > fess/stats/temp.energy
find $OUTPUT_DIR -name "temp.longrange.contact" | xargs cat > fess/stats/temp.longrange.contact
find $OUTPUT_DIR -name "temp.longrange.all" | xargs cat > fess/stats/temp.longrange.all
