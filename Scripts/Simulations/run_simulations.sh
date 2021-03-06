#!/bin/bash

# #PBS -l walltime=12:00:00,mem=16gb,nodes=1:ppn=8
#PBS -l walltime=24:00:00,mem=72gb,nodes=1:ppn=32
# #PBS -N gencor_recurrent_selection_simulation
# #PBS -N gencor_selection_simulation
#PBS -N gencor_prediction_simulation
# #PBS -N gencor_prediction_space_simulation
#PBS -M neyha001@umn.edu
#PBS -m abe
#PBS -r n

# Change the working directory
cd /panfs/roc/groups/6/smithkp/neyha001/Genomic_Selection/GenCorPrediction/Scripts/Simulations

module load R/3.5.0

# For genetic correlation
Rscript popvar_gencor_simulation.R

# For genetic correlation and recurrent selection
# Rscript popvar_gencor_recurrent_selection_simulation.R

