## GenCorPrediction
## 
## Script to make predictions of the genetic correlation using PopVar.
## 
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


## Set seed
set.seed(1009)

# Apply the function by core
gencor_popvar_list <- mclapply(X = all_par_crossing_block_split, FUN = function(core_df) {
    
    # Return predictions
    out <- pop.predict(G.in = G_in, y.in = phenos_use, map.in = map_use, 
                crossing.table = core_df, tail.p = 0.1, nInd = 150,
                min.maf = 0, mkr.cutoff = 1, entry.cutoff = 1, remove.dups = FALSE,
                nSim = 25, nCV.iter = 1, models = "rrBLUP", impute = "pass")
    
    tidy.popvar(out)

  }, mc.cores = n_cores)

gencor_popvar_out <- bind_rows(gencor_popvar_list)


## Save
save_file <- file.path(result_dir, "gencor_popvar_prediction_results.RData")
save("gencor_popvar_out", file = save_file)
 
# 
# ## Is there a difference in predicitons when simulation small crosses many times or larger crosses few times?
# gencor_popvar_small_fam <- pop.predict(G.in = G_in, y.in = phenos_use, map.in = map_use, 
#                                        crossing.table = crossing_block, tail.p = 0.1, nInd = 150,
#                                        min.maf = 0, mkr.cutoff = 1, entry.cutoff = 1, remove.dups = FALSE,
#                                        nSim = 25, nCV.iter = 1, models = "rrBLUP", impute = "pass")
# 
# gencor_popvar_small_fam_tidy <- tidy.popvar(gencor_popvar_small_fam) %>%
#   select(Par1, Par2, trait, contains("cor")) %>% 
#   rename_at(vars(contains("cor")), funs(str_split(string = ., pattern = "_") %>% map(3))) %>% 
#   gather(trait2, correlation, -Par1:-trait) %>% 
#   filter(!is.na(correlation))
# 
# gencor_popvar_large_fam <- pop.predict(G.in = G_in, y.in = phenos_use, map.in = map_use, 
#                                        crossing.table = crossing_block, tail.p = 0.1, nInd = 1000,
#                                        min.maf = 0, mkr.cutoff = 1, entry.cutoff = 1, remove.dups = FALSE,
#                                        nSim = 1, nCV.iter = 1, models = "rrBLUP", impute = "pass")
# 
# gencor_popvar_large_fam_tidy <- tidy.popvar(gencor_popvar_large_fam) %>%
#   select(Par1, Par2, trait, contains("cor")) %>% 
#   rename_at(vars(contains("cor")), funs(str_split(string = ., pattern = "_") %>% map(3))) %>% 
#   gather(trait2, correlation, -Par1:-trait) %>% 
#   filter(!is.na(correlation))
# 
# ## Correlate
# full_join(gencor_popvar_small_fam_tidy, gencor_popvar_large_fam_tidy, by = c("Par1", "Par2", "trait", "trait2")) %>% 
#   group_by(trait, trait2) %>%
#   summarize(prec = cor(correlation.x, correlation.y))
# 
# ## Plot
# full_join(gencor_popvar_small_fam_tidy, gencor_popvar_large_fam_tidy, by = c("Par1", "Par2", "trait", "trait2")) %>% 
#   qplot(x = correlation.x, y = correlation.y, data = .) +
#   facet_grid(trait ~ trait2)







