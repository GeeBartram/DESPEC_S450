#!/bin/bash
 
NAME=52
 
 cd /lustre/gamma/gbartram/DESPEC_S450_NEARLINE/Cluster_Submission/Nearline_Histograms/Run_${NAME}/
 rm -rf run${NAME}_merged.root
 rm -rf slurm*
 hadd -j 20 run${NAME}_merged.root run${NAME}*

