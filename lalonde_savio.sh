#!/bin/bash
# Job name: 
#SBATCH --job-name=lalonde
#
# Partition:
#SBATCH --partition=savio
#
# Account:
#SBATCH --account=ac_biostat
#
# QoS:
#SBATCH --qos=condo_biostat
#
# Processors:
#SBATCH --ntasks=20
#
# Memory requirement:
###SBATCH --mem-per-cpu=6G
#
# Wall clock limit:
#SBATCH --time=36:00:00
#
# Mail type:
#SBATCH --mail-type=all
#
# Mail user:
#SBATCH --mail-user=kecolson@berkeley.edu
#
## Run command
module unload intel
module load gcc openmpi r Rmpi
mpirun -n 1 R --vanilla < run_analysis_on_lalonde.R > run_analysis_on_lalonde.Rout
