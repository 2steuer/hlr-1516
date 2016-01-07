#!/bin/bash
shopt -n nullglob
FILES=slurm2/*.job
for f in $FILES
do
        sbatch $f
done
