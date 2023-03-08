#!/bin/bash
 #SETUP ENVIORNMENT
 source /cvmfs/eel.gsi.de/bin/go4login
 
NAME=14
 
LISTFILE="/lustre/gamma/gbartram/DESPEC_S450_NEARLINE/Cluster_Submission/List_Files/files_$NAME.txt"
 export ROOT_INCLUDE_PATH=/lustre/gamma/gbartram/DESPEC_S450_NEARLINE
 echo "DESPEC Kronos Started at date"

 #MAKES DIRECTORY FOR ROOT FILES
 OUTPUT_PATH=/lustre/gamma/gbartram/DESPEC_S450_NEARLINE/Cluster_Submission/Nearline_Histograms/Run_$NAME
 mkdir -p ${OUTPUT_PATH}

 NFILES=$(cat ${LISTFILE} | wc -l)
 echo "Analysing" $NFILES "Files"

##READ NAMES OF FILES IN ARRAY
 declare -a array
 while IFS= read -r line
 do
    array+=($line)
 done < "$LISTFILE"

 #FASTER WAY TO READ FILE NAMES
 #readarray -t array < $LISTFILE

 echo "Array is $SLURM_ARRAY_TASK_ID"
 part=(  "${array[@]:$SLURM_ARRAY_TASK_ID:2}" )

 echo "Running Go4!"

 #echo "$NAME" >> "$TEMP" 

 go4analysis -file ${part[*]} -enable-asf 1800 -asf ${OUTPUT_PATH}/run${NAME}_${SLURM_ARRAY_TASK_ID}.root

 echo "Run" $NAME "Completed :)"
