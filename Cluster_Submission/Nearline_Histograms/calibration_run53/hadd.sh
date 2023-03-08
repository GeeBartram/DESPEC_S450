#!/bin/bash

FILE="/lustre/gamma/gbartram/DESPEC_S450_NEARLINE/Cluster_Submission/Nearline_Histograms/calibration_run53/Include_Runs.txt"

while IFS= read -r line
 do
 RUN="$line"
 cd Run_$RUN
 hadd -j 20 run${RUN}_merged.root run$RUN*
 cd ..
done < "$FILE"
