#!/bin/bash
#$ -cwd
#$ -o "log/impute18.o"
#$ -e "log/impute18.e"
#$ -j y
#$ -m e
#$ -v "resnum=18"

source /etc/profile.d/modules.sh

## Load app
module load gcc/9.2.0
module load R/4.0.4

## Edit with your job command
Rscript impute_model.R

