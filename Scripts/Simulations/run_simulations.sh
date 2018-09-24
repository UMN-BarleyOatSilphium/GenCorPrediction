#!/bin/bash

#PBS -l walltime=24:00:00,mem=22gb,nodes=1:ppn=8
# #PBS -l walltime=12:00:00,mem=62gb,nodes=1:ppn=16
# #PBS -N PopVar_suitability_simulation
# #PBS -N PopVar_gencor_selection_simulation
#PBS -N PopVar_gencor_prediction_simulation
#PBS -M neyha001@umn.edu
#PBS -m abe
#PBS -r n

# Change the working directory
cd /panfs/roc/groups/6/smithkp/neyha001/Genomic_Selection/PopVarVal/Scripts/Simulations

module load R/3.5.0

## Run the simulation
# Rscript popvar_suitability_simulation.R

# For genetic correlation
Rscript popvar_gencor_simulation.R

# For genetic correlation and selection
# Rscript popvar_gencor_selection_simulation.R
