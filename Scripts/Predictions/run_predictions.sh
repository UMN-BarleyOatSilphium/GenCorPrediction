#!/bin/bash

#PBS -l walltime=40:00:00,mem=62gb,nodes=1:ppn=24
#PBS -N popvar_predictions_full
#PBS -M neyha001@umn.edu
#PBS -m abe
#PBS -r n

# Change the working directory
cd /panfs/roc/groups/6/smithkp/neyha001/Genomic_Selection/PopVarVal/Scripts/Predictions

module load R/3.5.0

# Prediction of all BP families
Rscript popvar_predictions.R
