#!/bin/bash

#find output/ -name "temp.angles" | xargs cat > stats/angles.csv
find output/ -name "temp.angles" | xargs cat > stats/temp.stats
find output/ -name "temp.energy" | xargs cat > stats/temp.energy
find output/ -name "temp.longrange.contact" | xargs cat > stats/temp.longrange.contact
find output/ -name "temp.longrange.all" | xargs cat > stats/temp.longrange.all
