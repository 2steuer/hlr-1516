#!/bin/bash

if [ ! -f "ergebnisse" ]
then 
	touch ergebnisse
fi

for ((i=7;i<11;i++)); do
	for ((j=1;j<13;j++)); do
		echo "run: $i threads: $j" >> ergebnisse;
		cat t_${i}_${j}.out | grep zeit >> ergebnisse;
	done
done
