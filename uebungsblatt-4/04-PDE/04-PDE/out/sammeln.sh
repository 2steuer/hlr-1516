#!/bin/bash

if [ ! -f "ergebnisse" ]
then 
	touch ergebnisse
fi

for ((i=1;i<11;i++)); do
	for ((j=1;j<13;j++)); do
		file=t_${i}_${j}.out
                if [ -f "$file" ]
                        then 
                        echo "run: $i threads: $j" >> ergebnisse;
                        cat ${file} | grep zeit >> ergebnisse;
                fi
	done
done
