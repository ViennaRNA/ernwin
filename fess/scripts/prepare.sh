#!/bin/bash

#set -x

ROSETTA_TOOLS_DIR=/scr/plastilin/pkerp/apps/rosetta_TRUNK/rosetta_tools/
MCANNOTATE_BIN=/scr/plastilin/pkerp/apps/mcannotate/MC-Annotate
MCANNOTATE_TO_DOTPLOT=../assembly-tests/scripts/mcannotate-to-dot-plot.py
SCRIPT_DIR=../assembly-tests/scripts
LOCAL_SCRIPT_DIR=scripts
K2N_PATH=/scr/plastilin/pkerp/apps/k2n_standalone/knotted2nested.py

base=${1##*/}
pdb=${base%\.*}

OUTPUT_DIR=output/${pdb}/prepare

#echo PYTHONPATH, $PYTHONPATH
export PYTHONPATH=$PYTHONPATH:/home/mescalin/pkerp/projects/fess


function convert_original_pdb_to_rosetta_ready {
    base=${1##*/}
    pdb=${base%\.*}

    # convert the regular pdb structure to rosetta's RNA pdb format
    $SCRIPT_DIR/make_rna_rosetta_ready.py $1 2> $OUTPUT_DIR/out.err > $OUTPUT_DIR/temp.pdb

    # extract the sequence
    #${SCRIPT_DIR}/pdb_to_seq.py ${1}_RNA.pdb > tests/${1}/rosetta_inputs/add.fasta

    #We're only interested in the first chain of the molecule
    ./$LOCAL_SCRIPT_DIR/get_biggest_chain.py $OUTPUT_DIR/temp.pdb $OUTPUT_DIR/temp.pdb.1
    $SCRIPT_DIR/make_rna_rosetta_ready.py $OUTPUT_DIR/temp.pdb.1 2> $OUTPUT_DIR/out.err > $OUTPUT_DIR/temp.pdb
}


function annotate_secondary_structure {
    base=${1##*/}
    pdb=${base%\.*}

    # Replace the rosetta-style nucleotide names with the original ones
    # so that we can annotate the secondary structure

    sed -i 's/rG/ G/g' $OUTPUT_DIR/temp.pdb 
    sed -i 's/rC/ C/g' $OUTPUT_DIR/temp.pdb 
    sed -i 's/rU/ U/g' $OUTPUT_DIR/temp.pdb 
    sed -i 's/rA/ A/g' $OUTPUT_DIR/temp.pdb 
    
    $MCANNOTATE_BIN $OUTPUT_DIR/temp.pdb > $OUTPUT_DIR/temp.mcannotate
    #create a dot plot without pseudoknots
    #$SCRIPT_DIR/mcannotate-to-dot-plot.py temp.pdb temp.mcannotate > temp.dotplot
    #$SCRIPT_DIR/mcannotate-to-bpseq.py temp.pdb temp.mcannotate > temp.bpseq
    ./$LOCAL_SCRIPT_DIR/mcannotate_to_bpseq.py $OUTPUT_DIR/temp.mcannotate > $OUTPUT_DIR/temp.bpseq

    python $K2N_PATH -f bpseq -F vienna $OUTPUT_DIR/temp.bpseq > $OUTPUT_DIR/temp.ss
    cat $OUTPUT_DIR/temp.ss | tail -n 2 | tr '.' 'x' > $OUTPUT_DIR/temp.tofold
    cat $OUTPUT_DIR/temp.ss | tail -n 2 | head -n 1 > $OUTPUT_DIR/temp.seq
    cat $OUTPUT_DIR/temp.ss | tail -n 1 | tr '.' 'x' > $OUTPUT_DIR/temp.const
    cat $OUTPUT_DIR/temp.ss | tail -n 1 > $OUTPUT_DIR/temp.dotplot
    #./get-helical-regions.py temp.pdb temp.mcannotate > ss.py

}

if [ $# -ne 1 ]; then
    echo "Usage.sh: ./enumerate.sh pdb_file"
    exit
fi

if [ -f $OUTPUT_DIR/temp.pdb ]; then
    rm $OUTPUT_DIR/temp.pdb
fi

echo 'processing', $1

convert_original_pdb_to_rosetta_ready $1
annotate_secondary_structure $1

res=`grep '[()]' $OUTPUT_DIR/temp.dotplot`

#if [ -z $res ]; then
    #echo removing, $1
    #mv $1 useless
#fi

#exit
echo $1, $(cat $OUTPUT_DIR/temp.dotplot) >&2

mv $OUTPUT_DIR/temp.dotplot $OUTPUT_DIR/temp.dotplot.2
./scripts/remove_length_one_stems.py $OUTPUT_DIR/temp.dotplot.2 > $OUTPUT_DIR/temp.dotplot

rm $OUTPUT_DIR/temp.dotplot.2

