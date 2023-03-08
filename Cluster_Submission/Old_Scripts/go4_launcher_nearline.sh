#!/bin/bash

##SETUP ENVIORNMENT
source /cvmfs/eel.gsi.de/bin/go4login
export ROOT_INCLUDE_PATH=/lustre/gamma/gbartram/DESPEC_S450_NEARLINE
echo "DESPEC Kronos Started at `date`"

##SET DATA LOCATION
#dpath=~/lustre/gamma/d004/ts/aida/

##READ IN LIST OF FILES TO RUN
#LISTFILE="u/gbartram/DESPEC_S450_NEARLINE/Cluster_Submission/file_list_Ir.txt"
LISTFILE="/lustre/gamma/gbartram/DESPEC_S450_NEARLINE/Cluster_Submission/file_list_Os.txt"

##COUNT NUMBER OF FILES
NFILES=$(cat ${LISTFILE} | wc -l)
echo "Analysing" $NFILES "Files"

OUTPUT_PATH=/lustre/gamma/gbartram/DESPEC_S450_NEARLINE/Cluster_Submission/Nearline_Histograms/Run_Os

mkdir -p ${OUTPUT_PATH}

##READ NAMES OF FILES IN ARRAY
declare -a array
while IFS= read -r line
do
    array+=($line)
done < "$LISTFILE"

echo "Array is $SLURM_ARRAY_TASK_ID"
part=(  "${array[@]:$SLURM_ARRAY_TASK_ID:2}" ) # :5 number of files to put together -> Has to be the same in the 2 .sh scripts

echo "Running Go4!"

#go4analysis -file ${part[*]} -enable-asf 1800 -asf /lustre/gamma/DESPEC_S450_NEARLINE/Cluster_Submission/Nearline_Histograms/Run_f0043_f0072/Files_202Os_${SLURM_ARRAY_TASK_ID}.root
go4analysis -file ${part[*]} -enable-asf 1800 -asf ${OUTPUT_PATH}/Files_202Os_${SLURM_ARRAY_TASK_ID}.root



