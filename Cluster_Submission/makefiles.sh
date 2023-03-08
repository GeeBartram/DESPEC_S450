#!/bin/bash

#CREATES LIST OF FILES
FILE="/lustre/gamma/gbartram/DESPEC_S450_NEARLINE/Cluster_Submission/Include_Runs.txt"
DATA_PATH=/lustre/gamma/DESPEC_S450_FILES/ts

while IFS= read -r line
 do
 RUN=($line)
 if [ ! -f /lustre/gamma/gbartram/DESPEC_S450_NEARLINE/Cluster_Submission/List_Files/files_$RUN.txt ]; then
  ls $DATA_PATH/s450f00$RUN* > /lustre/gamma/gbartram/DESPEC_S450_NEARLINE/Cluster_Submission/List_Files/files_$RUN.txt
 fi
 LISTFILE="/lustre/gamma/gbartram/DESPEC_S450_NEARLINE/Cluster_Submission/List_Files/files_$RUN.txt"


done < "$FILE"


