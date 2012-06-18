#!/bin/bash

#set -x

ROSETTA_TOOLS_DIR=/scr/plastilin/pkerp/apps/rosetta_TRUNK/rosetta_tools/
MCANNOTATE_BIN=/scr/plastilin/pkerp/apps/mcannotate/MC-Annotate
MCANNOTATE_TO_DOTPLOT=../assembly-tests/scripts/mcannotate-to-dot-plot.py
SCRIPT_DIR=../assembly-tests/scripts
K2N_PATH=/scr/plastilin/pkerp/apps/k2n_standalone/knotted2nested.py
OUTPUT_DIR=output


function convert_original_pdb_to_rosetta_ready {
    base=${1##*/}
    pdb=${base%\.*}

    # convert the regular pdb structure to rosetta's RNA pdb format
    $SCRIPT_DIR/make_rna_rosetta_ready.py $1 2> $OUTPUT_DIR/out.err > $OUTPUT_DIR/temp.pdb

    # extract the sequence
    #${SCRIPT_DIR}/pdb_to_seq.py ${1}_RNA.pdb > tests/${1}/rosetta_inputs/add.fasta

    #We're only interested in the first chain of the molecule
    ./get_biggest_chain.py $OUTPUT_DIR/temp.pdb $OUTPUT_DIR/temp.pdb.1
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
    ./mcannotate_to_bpseq.py $OUTPUT_DIR/temp.mcannotate > $OUTPUT_DIR/temp.bpseq

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

res=`grep '[()]' temp.dotplot`

#if [ -z $res ]; then
    #echo removing, $1
    #mv $1 useless
#fi

#exit
echo $1, $(cat temp.dotplot) >&2

base=${1##*/}
pdb=${base%\.*}

mv $OUTPUT_DIR/temp.dotplot $OUTPUT_DIR/temp.dotplot.2
./remove_length_one_stems.py $OUTPUT_DIR/temp.dotplot.2 > $OUTPUT_DIR/temp.dotplot
./create_bulge_graph.py $OUTPUT_DIR/temp.dotplot $pdb >  $OUTPUT_DIR/temp.bulge
./check_bulge_integrity.py $OUTPUT_DIR/temp.bulge
./graph_to_neato.py $OUTPUT_DIR/temp.bulge > $OUTPUT_DIR/temp_bulge.neato
neato -Tpng $OUTPUT_DIR/temp_bulge.neato -o $OUTPUT_DIR/temp_bulge.png
./bulge_to_fork_graph.py $OUTPUT_DIR/temp.bulge > $OUTPUT_DIR/temp.graph
./graph_to_neato.py $OUTPUT_DIR/temp.graph > $OUTPUT_DIR/temp_graph.neato
./check_bulge_integrity.py $OUTPUT_DIR/temp.graph
neato -Tpng $OUTPUT_DIR/temp_graph.neato -o $OUTPUT_DIR/temp_graph.png
./graph_to_angles.py $OUTPUT_DIR/temp.graph $OUTPUT_DIR/temp.pdb > $OUTPUT_DIR/temp.angles 
./graph_to_pymol.py $OUTPUT_DIR/temp.graph > $OUTPUT_DIR/cartoon.pml
#./graph_to_com_pymol.py temp.graph temp.pdb temp.mcannotate 1 > ss.py
./graph_to_coordinates.py $OUTPUT_DIR/temp.graph $OUTPUT_DIR/temp.pdb > $OUTPUT_DIR/temp.coords
./coordinates_to_pymol.py $OUTPUT_DIR/temp.coords > $OUTPUT_DIR/ss.py

#cat temp.bulge | awk -f graph_to_neato.awk > temp_bulge.neato
#cat temp.graph | awk -f graph_to_neato.awk > temp_graph.neato

#./graph_to_simplified_pymol.py temp.graph temp.pdb temp.mcannotate > ss.py
#./graph_to_com_pymol.py temp.graph temp.pdb temp.mcannotate > ss.py

./get_long_range_interactions.py $OUTPUT_DIR/temp.graph $OUTPUT_DIR/temp.mcannotate > $OUTPUT_DIR/temp_long_range.neato 2>> $OUTPUT_DIR/longrange.stats
cat temp_graph.neato $OUTPUT_DIR/temp_long_range.neato | awk -f $OUTPUT_DIR/combine-neato.awk > $OUTPUT_DIR/temp_combined.neato
/scr/plastilin/pkerp/apps/ViennaRNA/bin/RNAfold -C < $OUTPUT_DIR/temp.tofold
./output_protein_chains.py $1 $OUTPUT_DIR/temp_protein.pdb 2> /dev/null

neato -Tpng $OUTPUT_DIR/temp_combined.neato -o $OUTPUT_DIR/temp_combined.png
neato -Teps $OUTPUT_DIR/temp_bulge.neato -o $OUTPUT_DIR/temp_bulge.eps
neato -Teps $OUTPUT_DIR/temp_graph.neato -o $OUTPUT_DIR/temp_graph.eps
neato -Teps $OUTPUT_DIR/temp_combined.neato -o $OUTPUT_DIR/temp_combined.eps


./graph_to_angles.py $OUTPUT_DIR/temp.graph $OUTPUT_DIR/temp.pdb >> $OUTPUT_DIR/angles.csv
./bulge_to_distances.py $OUTPUT_DIR/temp.bulge $OUTPUT_DIR/temp.pdb > $OUTPUT_DIR/temp.distances
./graph_to_distances.py $OUTPUT_DIR/temp.graph $OUTPUT_DIR/temp.pdb >> $OUTPUT_DIR/temp.gdistances




#./graph_to_junction_orientations.py temp.graph temp.pdb >> temp.horientations

#./graph_to_constraints.py temp.bulge adj_lengths.csv stem_lengths.csv bulge_lengths.csv trans_lengths.csv > dg.data
#dgsol -s20 > dgsol.err
#./get_best_sol.py dg.sol > out.csv
#./solution_to_pymol.py temp.bulge out.csv > distance.py
#echo 'run distance.py' > distance.pml

num_angles=$(cat temp.angles | wc -l)
#echo "hide cartoon, temp" >> cartoon.pml

#echo "num_angles: $num_angles"

#if [ $num_angles -gt 0 ]; then
    #cat temp.angles | xargs -n 8 ./helix_angle.py temp.pdb
#fi

#print $params
#cmd="./helix_angle.py temp.pdb $params"
#$cmd


