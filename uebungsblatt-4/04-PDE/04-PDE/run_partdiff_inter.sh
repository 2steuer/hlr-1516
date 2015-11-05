#!/bin/bash

HOST=$(hostname)

export OMP_NUM_THREADS=12

./partdiff-openmp-zeile 12 2 $1 2 2 100 > out_inter/t_${2}_${1}.out

echo "[$2] done with 12 cores and $1 interlines on ${HOST}"
