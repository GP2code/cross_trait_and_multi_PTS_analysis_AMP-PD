#!/bin/bash
#$ -N mungesumstats_for_TRS_project
#$ -pe smp 1
#$ -l h_vmem=100G
#$ -l h_rt=72:0:0
#$ -cwd
#$ -j y
#$ -o /path/to/home/TRS_project/sumstats/logs/mungesumstats_for_TRS_project.log

module load R/4.4.1

Rscript /path/to/home/TRS_project/sumstats/scripts/MungeSumstats_for_TRS_project.R


## Submitted using: qsub /path/to/home/TRS_project/sumstats/scripts/submit_MungeSumstats_for_TRS_project.sh