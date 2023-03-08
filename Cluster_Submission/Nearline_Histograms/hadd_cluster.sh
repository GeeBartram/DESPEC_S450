#!/bin/bash

#CREATES LIST OF FILES
FILE="/lustre/gamma/gbartram/DESPEC_S450_NEARLINE/Cluster_Submission/Include_Runs.txt"

IFS=","

while read -r line z1shift z2shift; do 

 RUN="$line"

 echo '#!/bin/bash
 ' > hadds/hadd_$RUN.sh

 echo "NAME=$RUN
 " >> hadds/hadd_$RUN.sh

 echo ' cd /lustre/gamma/gbartram/DESPEC_S450_NEARLINE/Cluster_Submission/Nearline_Histograms/Run_${NAME}/
 rm -rf run${NAME}_merged.root
 rm -rf slurm*
 hadd -j 20 run${NAME}_merged.root run${NAME}*
' >> hadds/hadd_$RUN.sh
 

#SUBMIT JOB


  echo "hadding Run_$RUN..."
  
  sbatch -J gbartram_hadd_s450 -D  /lustre/gamma/gbartram/DESPEC_S450_NEARLINE/ -o logs/hadd_%A_%a.out.log -e logs/hadd_%A_%a.err.log \
  --time=8:00:00 --mem-per-cpu=4G \
  -- /lustre/gamma/gbartram/DESPEC_S450_NEARLINE/Cluster_Submission/Nearline_Histograms/hadds/hadd_$RUN.sh

 unset size


done < "$FILE"
