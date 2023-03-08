#!/bin/bash

FILE="/lustre/gamma/gbartram/DESPEC_S450_NEARLINE/Cluster_Submission/Include_Runs_Os.txt"

IFS=","

while read -r line z1shift z2shift; do
 RUN="$line"
 cp Run_$RUN/run${RUN}_merged.root merge_os
done < "$FILE"

cd merge_os
hadd -j 20 merge_os.root run*_merged.root
