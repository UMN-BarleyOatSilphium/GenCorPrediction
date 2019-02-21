## Genetic correlation permutation test
## 
## This script will generate a null distribution of the genetic correlation between traits using a permutation test
## 
## Author: Jeff Neyhart
## Date: February 21, 2019
## 

# Run the source script
repo_dir <- "/panfs/roc/groups/6/smithkp/neyha001/Genomic_Selection/GenCorPrediction/"
source(file.path(repo_dir, "source_MSI.R"))


# # Run on a local machine
# repo_dir <- getwd()
# source(file.path(repo_dir, "source.R"))
# # Additional libraries
# library(modelr)
# library(EMMREML)



## Number of cores
n_cores <- 24
n_cores <- detectCores()

# Load the PVV BLUEs
load(file.path(data_dir, "PVV_BLUE.RData"))

## Create a matrix of trait pairs
trait_pairs <- combn(x = traits, m = 2, simplify = F)


## Calculate the genetic correlation in the TP
K <- A.mat(X = s2_imputed_mat_use[tp_geno,], min.MAF = 0, max.missing = 1)

tp_tomodel <- trait_pairs %>%
  map(~filter(tp_prediction_BLUE, trait %in% ., line_name %in% tp_geno) %>% spread(trait, value))

# Model
tp_corG <- tp_tomodel %>% 
  map_df(~{
    df <- .
    
    f <- as.formula(paste("cbind(", paste(tail(names(df), 2), collapse = ","), ") ~ line_name"))
    mf <- model.frame(f, df)
    
    Y <- model.response(mf)
    X <- model.matrix(~ 1, mf)
    Z <- model.matrix(~ -1 + line_name, mf)
    
    # Fit the model
    fit <- emmremlMultivariate(Y = t(Y), X = t(X), Z = t(Z), K = K)
    # fit <- sommer::mmer(Y = Y, X = X, Z = list(g = list(Z = Z, K = K)))
    
    
    ## Return variance components
    varcomp <- fit$Vg
    # varcomp <- fit$var.comp$g
    data.frame(trait1 = colnames(Y)[1], trait2 = colnames(Y)[2], correlation = varcomp[1,2] / prod(sqrt(diag(varcomp))),
               row.names = NULL, stringsAsFactors = FALSE)
    
  })



## Generate permutations

# Number of permutations
nPerm <- 1000

set.seed(1224)
permutations <- tp_tomodel %>% 
  map(~rename_at(., vars(-line_name), funs(c("trait1", "trait2")))) %>% 
  map2(.x = ., .y = trait_pairs, ~mutate(.x, traits = paste0(.y, collapse = "/"))) %>%
  map_df(~permute(data = ., n = nPerm, trait1, trait2) %>% mutate(traits = map_chr(perm, ~as.data.frame(.)$traits[1])))


## Assign cores and split
permutations_split <- permutations %>%
  assign_cores(n_cores) %>%
  split(.$core)

## Iterate over cores
permutation_out <- mclapply(X = permutations_split, mc.cores = n_cores, FUN = function(core_df) {
  
  results_out <- vector("list", nrow(core_df))
  
  ## Iterate over rows
  for (i in seq_along(results_out)) {
    
    perm_data <- core_df$perm[[i]]
    
    f <- cbind(trait1, trait2) ~ line_name
    mf <- model.frame(f, perm_data)
    
    Y <- model.response(mf)
    X <- model.matrix(~ 1, mf)
    Z <- model.matrix(~ -1 + line_name, mf)
    
    # Fit the model
    fit <- emmremlMultivariate(Y = t(Y), X = t(X), Z = t(Z), K = K)

    
    ## Return variance components
    varcompG <- fit$Vg
    varcompE <- fit$Ve

    results_out[[i]] <- data.frame(corG = varcompG[1,2] / prod(sqrt(diag(varcompG))), corE = varcompE[1,2] / prod(sqrt(diag(varcompE))))
    
  }
  
  
  # Append to df and return
  core_df %>% 
    mutate(out = results_out) %>% 
    select(traits, .id, out)
  
})


## Bind rows
corG_permutation_out <- bind_rows(permutation_out)

## Save
save_file <- file.path(result_dir, "tp_corG_permutation_out.RData")
save("corG_permutation_out", file = save_file)







