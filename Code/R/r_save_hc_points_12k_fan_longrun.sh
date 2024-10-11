#!/bin/bash
#
#SBATCH --job-name=save_hc_points_12k_fan_longrun
#SBATCH -o save_hc_points_12k_fan_longrun.out
#SBATCH -e save_hc_points_12k_fan_longrun.err
#
#SBATCH --partition=longrun
#SBATCH --cpus-per-task=16 #Note: you are always allocated an even number of cpus
#SBATCH --mem=32G
#SBATCH --time=10-00:00:0

module purge
module load GCC/11.2.0 OpenMPI/4.1.1
module load R/4.2.0   

Rscript HPC_save_hc_points_12k_fan_longrun.R