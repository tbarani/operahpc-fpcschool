#!/bin/bash
compil_dir=$MMM_LIB_DIR

mpirun -n 6 ./test-bubble --dmin 0.6 --pref 1e8 -v 1 -l $compil_dir/libBehaviour.so -f ./bubbles.txt -m ./mesh_5um_10_pct.msh -sf 0.2
#./test-bubble --dmin 0.6 --pref 1e8 -v 1 -l $compil_dir/libBehaviour.so -f ./bubbles.txt -m ./mesh_5um_10_pct.msh -sf 0.2
