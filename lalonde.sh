#!/bin/bash
#
#$ -cwd
#$ -V
#$ -j y
#$ -S /bin/bash
#$ -M kecolson@berkeley.edu
#$ -m beas
#
mpirun -n 1 R --vanilla < run_analysis_on_lalonde.R > run_analysis_on_lalonde.Rout