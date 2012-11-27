#/bin/bash

#export PYTHONPATH=$PYTHONPATH:/home/mescalin/pkerp/projects/

if [ $# -ne 1 ]; then
    echo "Usage.sh: ./enumerate.sh pdb_file"
    exit
fi

base=${1##*/}
pdb=${base%\.*}

if [ -d fess/output/$pdb ]; then
    rm -rf fess/output/$pdb
fi

mkdir -p fess/output/${pdb}/{prepare,graph,neato,report,stats,pymol}

./fess/scripts/prepare.sh $1
./fess/scripts/create_graph.sh $1
#./fess/scripts/create_pymol.sh $1
./fess/scripts/create_stats.sh $1

