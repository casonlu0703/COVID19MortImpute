#!/bin/bash
#$ -cwd
#$ -o "log/impute8.o"
#$ -e "log/impute8.e"
#$ -j y
#$ -m e
#$ -v "resnum=8"

source /etc/profile.d/modules.sh

## Load app
module load gcc/9.2.0
module load R/4.0.4

## Edit with your job command
Rscript impute_model.R

