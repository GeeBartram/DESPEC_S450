#!/bin/bash

FILE="/lustre/gamma/gbartram/DESPEC_S450_NEARLINE/Cluster_Submission/Include_Runs.txt"
DESTINATION="/u/gbartram/S450/Runs/calibration_run85"

while IFS= read -r line
 do
 RUN="$line"
 cp Run_$RUN/run${RUN}_merged.root $DESTINATION
done < "$FILE"

