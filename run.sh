#!/bin/bash

java -cp build TestGetGenomes -d data -g 8400,9050,10997,19515 -p 3,3,3,1 \
        -s subGenomeRegions.txt -o work

java -cp build TestGetContigInput -g 8400,9050,10997,19515 -w 2,2,2,3 -wa 4 \
        -i work/genomesInString_8400_9050_10997_19515.txt -o work

python mwmatching-io.py -i work/contigInput_8400_9050_10997_19515.txt -o work/contigOutput.txt

java -cp build TestGetContigOutputAndScaffoldInput \
        -mml 3 -p 3,3,3,1 -w 2,2,2,3 -g 8400,9050,10997,19515 \
        -co work/contigOutput.txt \
        -s work/subgenomeRangesInGeneOrder_8400_9050_10997_19515.txt \
        -gf work/genomesInString_8400_9050_10997_19515.txt -o work

python mwmatching-io.py -i work/scaffolds -o work/scaffolds_output
