#!/bin/bash

##Setup environment
source /cvmfs/eel.gsi.de/bin/go4login
export ROOT_INCLUDE_PATH=/lustre/gamma/DESPEC_S450_NEARLINE
echo "DESPEC Kronos Started at `date`"

##Set data location
#dpath=~/lustre/gamma/d004/ts/aida/

##Read in list of files to run. Format names seperated by space,tab,newline
#LISTFILE="/lustre/gamma/amistry/DESPEC_S452_NEARLINE_200421/Cluster_Submission/file_list_f187188_testfatgain.txt"

#LISTFILE="/lustre/gamma/amistry/DESPEC_S452_NEARLINE_200421/Cluster_Submission/file_list_full_f74_f75.txt"

##LISTFILE="/lustre/gamma/DESPEC_S450_NEARLINE/Cluster_Submission/file_list_Ir.txt"
LISTFILE="/lustre/gamma/DESPEC_S450_NEARLINE/Cluster_Submission/file_list_f0043-f0046.txt"

##Count number of files
NFILES=$(cat ${LISTFILE} | wc -l)
        echo "Analysing" $NFILES "Files"

##Read names from list file
        declare -a array
        while IFS= read -r line
        do
array+=($line)
        done < "$LISTFILE"

        echo "Array is $SLURM_ARRAY_TASK_ID"
        part=(  "${array[@]:$SLURM_ARRAY_TASK_ID:2}" ) # :5 number of files to put together -> Has to be the same in the 2 .sh scripts

########added by z.chen, 13 May,2022#########
        OUTPUTFILE="Files_202Os_$SLURM_ARRAY_TASK_ID.root"
        OUTPUTPATH="/lustre/gamma/DESPEC_S450_NEARLINE/Cluster_Submission/Nearline_Histograms/Run_f0043_f0046"

        if [ ! -d "$OUTPUTPATH" ]; then
        echo "path is not correct, please check"
        echo "$OUTPUTPATH"
        exit 1
        fi




        echo "Running Go4!"
        go4analysis -file ${part[*]} -enable-asf 1800 -asf ${OUTPUTPATH}/${OUTPUTFILE}
#go4analysis -file ${part[*]} -enable-asf 1800 -asf /lustre/gamma/DESPEC_S450_NEARLINE/Cluster_Submission/Nearline_Histograms/Run_f0043_f0046/Files_202Os_$SLURM_ARRAY_TASK_ID.root

#go4analysis -file ${part[*]} -enable-asf 1800 -asf /lustre/gamma/amistry/DESPEC_S452_NEARLINE_200421/Cluster_Submission/Nearline_Histograms/ANALYf74f75_$SLURM_ARRAY_TASK_ID.root

