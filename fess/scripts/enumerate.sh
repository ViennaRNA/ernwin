#/bin/bash

#export PYTHONPATH=$PYTHONPATH:/home/mescalin/pkerp/projects/

if [ $# -ne 1 ]; then
    echo "Usage.sh: ./enumerate.sh pdb_file"
    exit
fi

base=${1##*/}
pdb=${base%\.*}

echo pwd
echo $(pwd)

OUTPUT_DIR=~/data/ernwin/processed/${pdb}

if [ -d $OUTPUT_DIR ]; then
    rm -rf $OUTPUT_DIR
fi

mkdir -p ${OUTPUT_DIR}/{prepare,graph,neato,report,stats,pymol}

pdb_name=$(echo $1 | awk -F '_' '{print $4}')

./fess/scripts/prepare.sh $1
#cp /scr/plastilin/pkerp/data/ernwin/processed/${pdb_name}/prepare/temp.dotplot ${OUTPUT_DIR}/prepare
./fess/scripts/create_graph.sh $1
./fess/scripts/create_pymol.sh -p $1
./fess/scripts/create_stats.sh $1

