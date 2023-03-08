#!/bin/bash

FILE="/lustre/gamma/gbartram/DESPEC_S450_NEARLINE/Cluster_Submission/Include_Runs.txt"

while IFS= read -r line
 do
 RUN="$line"
 cd Run_$RUN
 hadd -f -j 20 run${RUN}_merged.root run$RUN*
 cd ..
done < "$FILE"
