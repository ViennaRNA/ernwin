#!/bin/bash

#find output/ -name "temp.angles" | xargs cat > fess/stats/angles.csv
find fess/output/ -name "temp.angles" | xargs cat > fess/stats/temp.stats
find fess/output/ -name "temp.energy" | xargs cat > fess/stats/temp.energy
find fess/output/ -name "temp.longrange.contact" | xargs cat > fess/stats/temp.longrange.contact
find fess/output/ -name "temp.longrange.all" | xargs cat > fess/stats/temp.longrange.all
