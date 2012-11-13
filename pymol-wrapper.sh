#!/bin/bash

PDBNAME=
DEBUG=
TEXT=
BASE_NORMALS=

while getopts "xtbp:" OPTION
do
    case $OPTION in
        x)
            set -x
            DEBUG=1
            ;;
        t)
            TEXT=-t
            ;;
        b)
            BASE_NORMALS=-b
            ;;
        p)
            PDBNAME=$OPTARG
            echo $PDBNAME
            ;;
    esac
done

if [[ -z $PDBNAME ]]
then
    echo "Usage.sh: ./pymol-wrapper.sh -p struct_id"
    exit 1
fi

if [[ -z $DEBUG ]]
then
    ./fess/scripts/create_pymol.sh $TEXT $BASE_NORMALS -p fess/structures/$PDBNAME.pdb
else
    ./fess/scripts/create_pymol.sh $TEXT $BASE_NORMALS -xp fess/structures/$PDBNAME.pdb
fi

pymol fess/output/$PDBNAME/prepare/temp.pdb fess/output/$PDBNAME/pymol/cartoon.pml > /dev/null
