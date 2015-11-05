#!/bin/bash

HOST=$(hostname)

export OMP_NUM_THREADS=$1

./partdiff-openmp-zeile $1 2 512 2 2 1024 > out/t_${2}_${1}.out

echo "[$2] done with $1 cores on ${HOST}"
