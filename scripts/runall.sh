#!/bin/bash

# SLURM directives
#SBATCH --partition=cpu
#SBATCH --ntasks=1
#SBATCH --mem=64G
#SBATCH --cpus-per-task=32
#SBATCH --time 00:45:00
#SBATCH --job-name="A_without23andme_CONSERV"
#SBATCH --output="outputNewGwasSub.out.txt"
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=sennar.pierp@gmail.com

#module load r/4.2.2-gcc-10.3.0-withx-rmath-standalone-python3+-chk-version
module load r/4.3.0-gcc-13.2.0-withx-rmath-standalone-python-3.11.6
#INPUTS here, 1 2 3 4 5

#CorrectOLS_Gro_withoutGC_nocovar_variable

Rscript ~/proj/NewGwasSub/code.R "A_without23andme_CONSERV" 5
#Rscript ~/proj/NewGwasSub/matrix.R
