#!/bin/bash

#find output/ -name "temp.angles" | xargs cat > stats/angles.csv
find output/ -name "temp.angles" | xargs cat > stats/temp.stats
