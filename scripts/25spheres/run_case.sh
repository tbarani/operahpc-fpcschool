#!/bin/bash

export MMM_LIB_DIR=$HOME/mm-opera-hpc/build/bubble/src
compil_dir=$MMM_LIB_DIR
exec_path=$HOME/mm-opera-hpc/build/bubble/test-bubble

mpirun -n 6 $exec_path --dmin 0.6 --pref 1e8 -v 1 -l $compil_dir/libBehaviour.so -f ./bubbles.txt -m ./mesh_5um_10_pct.msh -sf 0.2
#./test-bubble --dmin 0.6 --pref 1e8 -v 1 -l $compil_dir/libBehaviour.so -f ./bubbles.txt -m ./mesh_5um_10_pct.msh -sf 0.2
