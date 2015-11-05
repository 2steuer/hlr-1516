#!/bin/sh
./run_partdiff_inter.sh 1 $1
./run_partdiff_inter.sh 2 $1
./run_partdiff_inter.sh 4 $1
./run_partdiff_inter.sh 8 $1
./run_partdiff_inter.sh 16 $1
./run_partdiff_inter.sh 32 $1
./run_partdiff_inter.sh 64 $1
./run_partdiff_inter.sh 128 $1
./run_partdiff_inter.sh 256 $1
./run_partdiff_inter.sh 512 $1
./run_partdiff_inter.sh 1024 $1

