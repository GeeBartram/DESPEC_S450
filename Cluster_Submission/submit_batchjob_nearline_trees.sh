#!/bin/bash

#CREATES LIST OF FILES
FILE="/lustre/gamma/gbartram/DESPEC_S450_NEARLINE/Cluster_Submission/Include_Runs.txt"
DATA_PATH=/lustre/gamma/DESPEC_S450_FILES/ts


while IFS= read -r line
 do
 RUN="$line"

 echo "Creating files for Run_$RUN..."

  ls $DATA_PATH/s450f00$RUN* > /lustre/gamma/gbartram/DESPEC_S450_NEARLINE/Cluster_Submission/List_Files/files_${RUN}.txt

 echo '#!/bin/bash
 #SETUP ENVIORNMENT
 source /cvmfs/eel.gsi.de/bin/go4login
 ' > go4_launchers/go4_launcher_nearline_tree_$RUN.sh

 echo "NAME=$RUN
 " >> go4_launchers/go4_launcher_nearline_tree_$RUN.sh

 echo 'LISTFILE="/lustre/gamma/gbartram/DESPEC_S450_NEARLINE/Cluster_Submission/List_Files/files_$NAME.txt"
 export ROOT_INCLUDE_PATH=/lustre/gamma/gbartram/DESPEC_S450_NEARLINE
 echo "DESPEC Kronos Started at `date`"

 #MAKES DIRECTORY FOR TREES FILES
OUTPUT_PATH=/lustre/gamma/gbartram/DESPEC_S450_NEARLINE/Cluster_Submission/Trees/Run_$NAME
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

 go4analysis -file ${part[*]} -step 1 -store ${OUTPUT_PATH}/tree${NAME}_${SLURM_ARRAY_TASK_ID}.root 99 -disable-asf -maxtreesize 100g -step 2 -disable-step





 echo "Run" $NAME "Completed :)"' >> go4_launchers/go4_launcher_nearline_tree_$RUN.sh
 

#LOOKS FOR NO OF FILES
 LISTFILE="/lustre/gamma/gbartram/DESPEC_S450_NEARLINE/Cluster_Submission/List_Files/files_$RUN.txt"

 readarray -t size < $LISTFILE


#SUBMIT JOB


  echo "Processing Run_$RUN..."
 sbatch -p long -J gbartram_go4_s450 -D /lustre/gamma/gbartram/DESPEC_S450_NEARLINE/ -o logs/go4_%A_%a_trees.out.log -e logs/go4_%A_%a_trees.err.log \
  --time=7-00:00:00 --mem-per-cpu=4G \
 --array=0-${#size[@]}:2 -- /lustre/gamma/gbartram/DESPEC_S450_NEARLINE/Cluster_Submission/go4_launchers/go4_launcher_nearline_tree_$RUN.sh	

 unset size


done < "$FILE"
