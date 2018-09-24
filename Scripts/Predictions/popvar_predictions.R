## PopVar Predictions
## 
## These are predictions of genetic variance using all available C1 lines (to be
## listed below).
## 
## This script can be run locally or on MSI
## 

# Run the source script
repo_dir <- "/panfs/roc/groups/6/smithkp/neyha001/Genomic_Selection/PopVarVal/" 
source(file.path(repo_dir, "source_MSI.R"))

# # Run on a local machine
# repo_dir <- getwd()
# source(file.path(repo_dir, "source.R"))
# library(PopVar)

# Load genotypic and phenotypic data
load(file.path(geno_dir, "S2_genos_mat.RData"))
load(file.path(pheno_dir, "PVV_BLUE.RData"))


# Subset the genos
genos_use <- s2_imputed_mat[c(tp_geno, pot_pars_geno),]

# Reformat the genotypic data
# Create data.frame for output
G_in <- as.data.frame(cbind( c("", row.names(genos_use)), rbind(colnames(genos_use), genos_use)) )

# Format the phenos
phenos_use <- tp_prediction_BLUE %>%
# phenos_use <- tp_relevant_BLUE %>%
  spread(trait, value) %>%
  as.data.frame()

# Format the map
map_use <- snp_info %>% 
  select(marker = `rs#`, chrom, cM_pos) %>%
  as.data.frame()
  
## Subset the crosses that were made
cross_list <- entry_list %>%
  filter(Group == "Experimental") %>%
  separate(Pedigree, c("parent1", "parent2"), sep = "/") %>%
  rename_all(str_to_lower)

# Format the crossing block
crossing_block <- cross_list %>%
  distinct(family, parent1, parent2) %>%
  as.data.frame() %>%
  column_to_rownames("family")


# Determine the number of cores
n_cores <- detectCores()


## Predict all possible bp families

# Create the crossing block of all pp_geno individuals
# Remove the crosses that were already predicted
all_par_crossing_block <- combn(x = pot_pars_geno, m = 2) %>%
  t() %>%
  as.data.frame() %>%
  rename(parent1 = V1, parent2 = V2)


# Add cores to the crossing block and split by core
all_par_crossing_block_split <- all_par_crossing_block %>% 
  assign_cores(df = ., n_core = n_cores) %>%
  split(.$core)



# Apply the function by core
all_par_pred_out <- mclapply(X = all_par_crossing_block_split, FUN = function(core_df) {
    
    # Return predictions
    out <- pop.predict(G.in = G_in, y.in = phenos_use, map.in = map_use, 
                crossing.table = core_df, tail.p = 0.1, nInd = 150,
                min.maf = 0, mkr.cutoff = 1, entry.cutoff = 1, remove.dups = FALSE,
                nSim = 25, nCV.iter = 1, models = "rrBLUP", impute = "pass")
    
    tidy.popvar(out)

  }, mc.cores = n_cores)

all_family_pred <- bind_rows(all_par_pred_out)


## Save
save_file <- file.path(result_dir, "all_family_prediction_results.RData")
save("all_family_pred", file = save_file)





## Separate predictions of FHB severity using TP data from each location
phenos_use_FHB <- tp_prediction_BLUE_FHB %>%
  rename(FHBSeverity = value) %>%
  as.data.frame() %>%
  split(.$location)
  

cross_pred_out_FHB <- phenos_use_FHB %>%
  map_df(~{
    
    df <- .
    df1 <- select(df, -location)
    
    # Return predictions
    out <- pop.predict(G.in = G_in, y.in = df1, map.in = map_use, 
                       crossing.table = crossing_block, tail.p = 0.1, nInd = 150,
                       min.maf = 0, mkr.cutoff = 1, entry.cutoff = 1, remove.dups = FALSE,
                       nSim = 25, nCV.iter = 1, models = "rrBLUP", impute = "pass")
  
    out$predictions %>% 
      mutate_all(unlist) %>% 
      as_data_frame() %>%
      mutate(trait = "FHBSeverity", location = unique(df$location)) %>%
      select(trait, location, names(.))
    
  })


pred_results_FHB <- cross_pred_out_FHB

## Save
save_file <- file.path(result_dir, "prediction_results_FHB.RData")
save("pred_results_FHB", file = save_file)











# ### Predictions using different population sizes
# # Cross data
# 
# 
# ## List of TP sizes
# tp_size_list <- seq(25, 150, by = 25)
# n_samples <- 50
# 
# ## Generate TP samples
# tp_samples <- map(tp_size_list, ~replicate(n = n_samples, sort(sample(tp_geno, size = .))))
# 
# ## Iterate over the population sizes
# tp_size_predictions <- tp_samples %>%
#   map(~{
#     
#     sample_mat <- .
#     
#     apply(X = sample_mat, MARGIN = 2, FUN = function(tp_geno_sample) {
#     
#       # Create data.frame for output
#       G_in_sample <- G_in[c(tp_geno_sample, pot_pars_geno),]
#       
#       # Format the phenos
#       # phenos_use <- tp_prediction_BLUE %>% 
#       phenos_use <- tp_relevant_BLUE %>%
#         spread(trait, value) %>%
#         as.data.frame() %>%
#         filter(line_name %in% tp_geno_sample)
#       
#       # Format the map
#       map_use <- snp_info %>% 
#         select(marker = `rs#`, chrom, cM_pos) %>%
#         as.data.frame()
#   
#   
#       ## Predict
#       pred_out <- pop.predict(G.in = G_in_sample, y.in = phenos_use, map.in = map_use,
#                               crossing.table = crossing_block, tail.p = 0.1, nInd = 150,
#                               min.maf = 0, mkr.cutoff = 1, entry.cutoff = 1, remove.dups = FALSE,
#                               nSim = 25, nCV.iter = 1, models = "rrBLUP", impute = "pass")
#   
#   
#       # Convert to DF
#       tidy.popvar(pred_out) %>%
#         left_join(rownames_to_column(crossing_block, "family"), ., by = c("parent1" = "Par1", "parent2" = "Par2"))
#       
#     })
#     
#   })
# 
















