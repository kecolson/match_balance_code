#!/bin/bash
#
#$ -cwd
#$ -V
#$ -j y
#$ -S /bin/bash
#$ -M kecolson@berkeley.edu
#$ -m beas
#
mpirun -n 1 R --vanilla < match_balance_and_mse2_rmpi.R > match_balance_and_mse2_rmpi.Rout