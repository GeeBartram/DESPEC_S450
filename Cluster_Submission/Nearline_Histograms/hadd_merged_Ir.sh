#!/bin/bash

FILE="/lustre/gamma/gbartram/DESPEC_S450_NEARLINE/Cluster_Submission/Include_Runs_Ir.txt"
IFS=","

while read -r line z1shift z2shift; do 
 RUN="$line"
 cp Run_$RUN/run${RUN}_merged.root merge_ir
done < "$FILE"

cd merge_ir
hadd -j 20 merge_ir.root run*_merged.root
