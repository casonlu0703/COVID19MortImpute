#!/bin/bash
#$ -cwd
#$ -o "log/impute1.o"
#$ -e "log/impute1.e"
#$ -j y
#$ -m e
#$ -v "resnum=1"

source /etc/profile.d/modules.sh

## Load app
module load gcc/9.2.0
module load R/4.0.4

## Edit with your job command
Rscript impute_model.R

