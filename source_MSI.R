## GenCorPrediction
## 
## A script that automatically loads the data relevant for the GenCorPrediction project


# Load packages
packages <- c("dplyr", "tidyr", "tibble", "stringr", "readxl", "readr", "parallel",
              "purrr", "pbsim", "Matrix", "pbsimData", "rrBLUP", "EMMREML", "modelr")

## The packages pbsim and pbsimData are available from GitHub:
## pbsim: https://github.com/neyhartj/pbsim
## pbsimData: https://github.com/neyhartj/pbsimData


invisible(lapply(packages, library, character.only = TRUE))


## Directories
proj_dir <- "/path/to/directory/on/supercomputer"


geno_dir <-  "/path/to/genotype/data/on/supercomputer"
pheno_dir <- "/path/to/phenotype/data/on/supercomputer"

# Other directories
fig_dir <- file.path(proj_dir, "Figures")
data_dir <- file.path(proj_dir, "Data")
result_dir <- file.path(proj_dir, "Results")


## Source the 'source.R' script from a define starting point
source_lines <- readLines(file.path(repo_dir, "source.R"))
source_lines_discard <- seq(which(grepl(pattern = "^# MSI", x = source_lines)))
source_lines_run <- source(textConnection(source_lines[-source_lines_discard]))

