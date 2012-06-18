#/bin/bash

export PYTHONPATH=$PYTHONPATH:/home/mescalin/pkerp/projects/

if [ $# -ne 1 ]; then
    echo "Usage.sh: ./enumerate.sh pdb_file"
    exit
fi

base=${1##*/}
pdb=${base%\.*}

if [ -d output/$pdb ]; then
    rm -rf output/$pdb
fi

mkdir -p output/${pdb}/{prepare,graph,neato,report,stats,pymol}

./scripts/prepare.sh $1
./scripts/create_graph.sh $1
./scripts/create_pymol.sh $1
./scripts/create_stats.sh $1

