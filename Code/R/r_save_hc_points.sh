#!/bin/bash
#
#SBATCH --job-name=save_hc_points
#SBATCH -o save_hc_points.out
#SBATCH -e save_hc_points.err
#
#SBATCH --cpus-per-task=2 #Note: you are always allocated an even number of cpus
#SBATCH --mem=1G
#SBATCH --time=10:00

module purge
module load GCC/11.2.0 OpenMPI/4.1.1
module load R/4.2.0   

Rscript HPC_save_hc_points.R