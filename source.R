## GenCorPrediction
## 
## A script that automatically loads the data relevant for the GenCorPrediction project

library(tidyverse)
library(readxl)
library(rrBLUP)
library(neyhart)
library(boot)

## Directories
proj_dir <- repo_dir


# Geno, pheno, and enviro data
geno_dir <-  "path/to/folder/with/genotype/data" 

# Other directories
fig_dir <- file.path(proj_dir, "Figures")
data_dir <- pheno_dir <- file.path(proj_dir, "Data")
result_dir <- file.path(proj_dir, "Results")

# Create a vector of colors
all_colors <- colorRampPalette(umn_palette(2)[3:5])



######
# MSI Source starts here
######


# Source the project functions
source(file.path(proj_dir, "source_functions.R"))


# Load the genotypic data
load(file.path(geno_dir, "S2_genos_mat.RData")) # Available from the Triticeae Toolbox

# Relevant traits
traits <- c("FHBSeverity", "HeadingDate", "PlantHeight")
traits_replace <- setNames(c("FHB Severity", "Heading Date", "Plant Height"), traits)


# Load an entry file
entry_list <- read_csv(file.path(data_dir, "project_entries.csv"))
# Load the trial metadata
trial_info <- read_csv(file.path(data_dir, "trial_metadata.csv"))


# Grab the entry names that are not checks
tp <- subset(entry_list, Group == "S2TP", Line_name, drop = T)

pot_pars <- subset(entry_list, Group == "Potential_Parent", Line_name, drop = T)
pars <- subset(entry_list, Note == "Parent", Line_name, drop = T)

exper <- subset(entry_list, Group == "Experimental", Line_name, drop = T)

# Find the tp and vp that are genotypes
tp_geno <- intersect(tp, row.names(s2_imputed_mat))
pot_pars_geno <- intersect(pot_pars, row.names(s2_imputed_mat))
pars_geno <- intersect(pars, row.names(s2_imputed_mat))

# Define the checks
checks <- unique(subset(entry_list, str_detect(Group, "Check"), Line_name, drop = T))

# Extract the tp and vp from the G matrix
s2_imputed_mat_use <- s2_imputed_mat[c(tp_geno, pot_pars_geno),]


## Format a list of the crosses and their experiment name
cross_list <- entry_list %>% 
  filter(Group == "Experimental") %>% 
  separate(Pedigree, c("parent1", "parent2"), "/") %>% 
  rename_all(str_to_lower) %>% 
  distinct(family, parent1, parent2, note)

