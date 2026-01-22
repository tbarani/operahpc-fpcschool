#!/bin/bash

export MMM_LIB_DIR=$HOME/mm-opera-hpc/build/bubble/src
compil_dir=$MMM_LIB_DIR

#######################
phydro=0.15
porosity=0.096
internalPressure=1e8
###############

mpirun -n 6 ./BubbleImposedMacroDisplacements -f ./bubbles.txt -m ./mesh_5um_10_pct.msh --dmin 0.6 --pref $internalPressure -v 1 -l $compil_dir/libBehaviour.so -sf 0.2 -mg 1 -hp $phdro -por $porosity

