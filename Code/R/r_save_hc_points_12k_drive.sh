#!/bin/bash
#
#SBATCH --job-name=save_hc_points_12k_drive
#SBATCH -o save_hc_points_12k_drive.out
#SBATCH -e save_hc_points_12k_drive.err
#
#SBATCH --partition=parallel
#SBATCH --cpus-per-task=32 #Note: you are always allocated an even number of cpus
#SBATCH --mem=64G
#SBATCH --time=10-00:00:0

#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=rasidilhani@gmail.com

module purge
module load GCC/11.2.0 OpenMPI/4.1.1
module load R/4.2.0   

Rscript HPC_save_hc_points_12k_drive.R