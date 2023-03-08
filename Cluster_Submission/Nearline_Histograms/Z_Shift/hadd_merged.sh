#!/bin/bash

FILE="/lustre/gamma/gbartram/DESPEC_S450_NEARLINE/Cluster_Submission/Include_Runs.txt"

while IFS= read -r line
 do
 RUN="$line"
 cp Run_$RUN/run${RUN}_merged.root merge
done < "$FILE"

cd merge
hadd -j 20 Z_shift_all_merged.root run*_merged.root
