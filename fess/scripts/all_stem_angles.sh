#!/bin/bash

rm fess/stats/stem_angles.csv
find fess/output -name "temp.comp" | xargs -n 1 ./fess/scripts/stem_angles.py | tee -a fess/stats/stem_angles.csv
