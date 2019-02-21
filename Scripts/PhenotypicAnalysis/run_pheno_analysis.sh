#!/bin/bash

#PBS -l walltime=24:00:00,mem=48gb,nodes=1:ppn=24
#PBS -N tp_gencor_permutation
#PBS -M neyha001@umn.edu
#PBS -m abe
#PBS -r n

# Change the working directory
cd /panfs/roc/groups/6/smithkp/neyha001/Genomic_Selection/PopVarVal/Scripts/PhenotypicAnalysis

module load R/3.5.0

# Permutation test to estimate SE of genetic correlation
Rscript phenotypic_correlation_permutation.R