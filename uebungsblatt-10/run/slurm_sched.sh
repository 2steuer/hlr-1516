#!/bin/bash
shopt -n nullglob
for f in *.job
do
	sbatch f
done

