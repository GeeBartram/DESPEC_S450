#!/bin/bash

##Read in list of files to run
LISTFILE="/lustre/gamma/gbartram/DESPEC_S450_NEARLINE/Cluster_Submission/test/files_84.txt"
declare -a size
while IFS= read -r line
do
    size+=($line)
done < "$LISTFILE"

##Submit job

sbatch -J gbartram_go4_s450 -D /lustre/gamma/gbartram/DESPEC_S450_NEARLINE/ -o logs/go4_%A_%a.out.log -e logs/go4_%A_%a.err.log \
  --time=8:00:00 --mem-per-cpu=4G \
  --array=0-${#size[@]}:2 -- /lustre/gamma/gbartram/DESPEC_S450_NEARLINE/Cluster_Submission/go4_launcher_nearline_imp.sh

  unset size
