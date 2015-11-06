#!/bin/bash

if [ ! -f "ergebnisse" ]
then 
	touch ergebnisse
fi

#i nummer des laufs
#j threadanzahl

for ((i=1;i<10;i++)); do
	for ((j=1;j<13;j++)); do
		if [ ! -f "t_${i}_${j}.out" ]
			then 
			echo "run: $i threads: $j" >> ergebnisse;
			cat t_${i}_${j}.out | grep zeit >> ergebnisse;
		fi
	done
done
