## PopVarVal simulations of prediction accuracy for the genetic
## correlation between two traits
## 
## 
## Author: Jeff Neyhart
## Last modified: September 2, 2018
## 

# Run the source script
repo_dir <- "/panfs/roc/groups/6/smithkp/neyha001/Genomic_Selection/GenCorPrediction/"
source(file.path(repo_dir, "source_MSI.R"))

# Load the two-row simulation genotypes
load(file.path(geno_dir, "s2_cap_simulation_data.RData"))

## Check if the results are present - if so only simulate the missing combinations
check_results <- T


# # Run on a local machine
# repo_dir <- getwd()
# source(file.path(repo_dir, "source.R"))
# # Additional libraries
# library(pbsim)
# 
# # Load the two-row simulation genotypes
# load(file.path(gdrive_dir, "BarleyLab/Projects/SideProjects/Resources/s2_cap_simulation_data.RData"))


# Remove monomorphic SNPs
# Filter for breeding programs relevant to my data
s2_cap_genos <- s2_cap_genos[str_detect(string = row.names(s2_cap_genos), pattern = "AB|BA|WA|N2|MT"),]
s2_cap_genos <- s2_cap_genos[,!colMeans(s2_cap_genos) %in% c(0, 2)]
s2_snp_info <- subset(s2_snp_info, rs %in% colnames(s2_cap_genos))





# Number cores
n_cores <- 24 # Local machine for demo
n_cores <- detectCores()


#### Simulation to test different genetic architectures and training population parameters


## Fixed parameters
sim_pop_size <- 150
n_iter <- 50
n_env <- 3
n_rep <- 1
n_crosses <- 50
k_sp <- 1.76

## Outline the parameters to perturb
trait1_h2_list <- trait2_h2_list <- c(0.3, 0.6, 1)
nQTL_list <- c(30, 100)
tp_size_list <- seq(150, 600, by = 150)
gencor_list <- c(-0.5, 0, 0.5)
probcor_list <- data_frame(arch = c("pleio", "close_link", "loose_link"),
                           input = list(cbind(0, 1), cbind(5, 1), rbind(c(25, 0), c(35, 1)) ))
model_list <- c("RRBLUP", "BayesC")

# Create a data.frame of parameters
param_df <- crossing(trait1_h2 = trait1_h2_list, trait2_h2 = trait2_h2_list, nQTL = nQTL_list, tp_size = tp_size_list,
                     gencor = gencor_list, probcor = probcor_list, model = model_list, iter = seq(n_iter))

map_sim <- s2_snp_info %>%
  split(.$chrom) %>%
  map(~setNames(.$cM_pos, .$rs)) %>%
  map(~structure(., class = "A"))  %>%
  structure(., class = "map") %>%
  # Jitter
  qtl::jittermap(.) %>%
  `names<-`(., seq_along(.))

# Create the base genome and the pedigree - these are fixed for all simulations
genome <- sim_genome(map = map_sim)
# Pedigree to accompany the crosses
ped <- sim_pedigree(n.ind = sim_pop_size, n.selfgen = Inf)



## Check the results file
save_file <- file.path(result_dir, "popvar_gencor_simulation_prediction_results.RData")

# If it exists, load it and create the missing combinations
if (file.exists(save_file) & check_results) {
  load(save_file)

  missing <- popvar_prediction_simulation_out %>% 
    select(-input, -results) %>% 
    mutate_all(as.factor) %>% 
    anti_join(x = complete_(., names(.)), y = .) %>%
    mutate_all(as.character) %>%
    mutate_all(parse_guess)

  # Build a new parameter set
  param_df <- left_join(missing, param_df)

}




# Split the parameter df
param_df_split <- param_df %>%
  assign_cores(n_cores) %>%
  split(.$core)

# Parallelize
simulation_out <- mclapply(X = param_df_split, FUN = function(core_df) {

  # # For local machine
  # i <- 11
  # core_df <- param_df_split[[i]]
  # #


  # Create a results list
  results_out <- vector("list", nrow(core_df))

  # Iterate over the rows of the param_df
  for (i in seq_along(results_out)) {
  # for (i in seq(i2, length(results_out))) {

    trait1_h2 <- core_df$trait1_h2[i]
    trait2_h2 <- core_df$trait2_h2[i]
    L <- core_df$nQTL[i]
    maxL <- max(param_df$nQTL)
    tp_size <- core_df$tp_size[i]
    gencor <- core_df$gencor[i]
    probcor <- core_df$input[[i]]
    model <- core_df$model[i]
    # nMar <- core_df$marker_density[i]
    
    maxL <- ifelse(probcor[1] == 0, maxL + (L / 2), maxL)
    

    # Simulate QTL
    qtl_model <- replicate(n = 2, matrix(NA, ncol = 4, nrow = L), simplify = FALSE)
    genome1 <- sim_multi_gen_model(genome = genome, qtl.model = qtl_model, add.dist = "geometric", max.qtl = maxL,
                                   corr = gencor, prob.corr = probcor)

    ## Adjust the genetic architecture - only if pleiotropy is not present
    if (probcor[1] != 0) {
      genome1 <- adj_multi_gen_model(genome = genome1, geno = s2_cap_genos, gencor = gencor)
    }

    # ## Trim the markers to the desired density
    # marker_names <- markernames(genome1)
    # 
    # # Determine the number of markers to draw from each chromosome
    # # Weight each chromosome by the number of markers
    # nMarChr <- round((nmar(genome1, by.chr = T) / nmar(genome1)) * nMar)
    # # Add or subtract from the shortest or longest chromosomes to comply with nMar
    # if (sum(nMarChr) > nMar) {
    #   nsub <- sum(nMarChr) - nMar
    #   chr <- order(chrlen(genome1), decreasing = TRUE)[seq(nsub)]
    #   nMarChr[chr] <- nMarChr[chr] - 1
    #   
    # } else if (sum(nMarChr) < nMar) {
    #   nadd <- nMar - sum(nMarChr)
    #   chr <- order(chrlen(genome1), decreasing = FALSE)[seq(nadd)]
    #   nMarChr[chr] <- nMarChr[chr] + 1
    #   
    # }
    # 
    # # Get the marker map
    # marker_map <- find_markerpos(genome1, marker_names) %>% 
    #   table_to_map()
    # 
    # # Use K-means to determine the center of each cluster
    # marker_subset <- marker_map %>% 
    #   map2(.x = ., .y = nMarChr, ~kmeans(x = .x, centers = .y, nstart = 500, iter.max = 50)) %>%
    #   map("centers") %>% 
    #   map2(.x = ., .y = marker_map, function(.x, .y) {
    #     ctn <- .x
    #     cm <- .y
    #     map(ctn, ~which.min(abs(cm - .))) %>% map_chr(names) })
    # 
    # ## Add the markers to the QTL, then subset the population
    # markers_qtl <- pull_qtl(genome1, unique = T) %>% 
    #   split(.$chr) %>% 
    #   map("qtl_name") %>%
    #   map2(.x = ., .y = marker_subset, c)
    
    # Edit the genome
    genome2 <- genome1
    # genome2$map <- map2(.x = genome1$map, .y = markers_qtl, ~.x[.y]) %>% map(sort)
    
     
    # Create the TP by random selection
    tp1 <- create_pop(genome = genome2, geno = s2_cap_genos[sort(sample(nrow(s2_cap_genos), size = tp_size)), markernames(genome = genome2, include.qtl = TRUE)]) %>% 
      # Phenotype the base population
      sim_phenoval(pop = ., h2 = c(trait1_h2, trait2_h2), n.env = n_env, n.rep = n_rep) 

    
    # Measure the genetic variance, genetic covariance, and genetic correlation in the tp1
    tp_summ <- tp1$geno_val %>%
      mutate(cor = cor(trait1, trait2)) %>% 
      gather(trait, value, trait1, trait2) %>% 
      group_by(trait) %>% 
      summarize_at(vars(cor, value), funs(mean, var)) %>% 
      select(trait, mean = value_mean, var = value_var, cor = cor_mean)
      
    par_pop <- tp1

    # ##### Create a set of cycle 1 progeny from which to predict crosses
    # ## Select the best tp individuals for both traits
    # tp_select <- select_pop(pop = tp1, intensity = 0.1, index = c(1, 1), type = "phenotypic")
    # # Randomly create crosses from the TP individuals
    # crossing_block <- sim_crossing_block(parents = indnames(tp_select), n.crosses = 40)
    # # Pedigree to accompany the crosses
    # ped <- sim_pedigree(n.ind = 25, n.selfgen = Inf)
    # 
    # # Make theses crosses
    # par_pop <- sim_family_cb(genome = genome1, pedigree = ped, founder.pop = tp_select, crossing.block = crossing_block)
    # #####

    # Randomly create crosses from the cycle1 individuals
    crossing_block <- sim_crossing_block(parents = indnames(par_pop), n.crosses = n_crosses)

    
    ## Predict the marker effects
    tp1 <- pred_mar_eff(genome = genome2, training.pop = tp1, method = model, n.iter = 1500, burn.in = 500, 
                        save.at = file.path(proj_dir, "Scripts/Simulations/bglr_out/", paste0("out", core_df$core[i])))


    ## Predict genetic variance and correlation
    pred_out <- pred_genvar(genome = genome2, pedigree = ped, training.pop = tp1, founder.pop = par_pop,
                            crossing.block = crossing_block, method = model) %>%
      mutate(pred_musp = pred_mu + (k_sp * sqrt(pred_varG))) # %>%
      # # Add predictions of covariance
      # group_by(parent1, parent2) %>%
      # mutate(pred_covG = pred_corG[1] * (prod(sqrt(pred_varG)))) %>%
      # ungroup()



    # ## Check the predictions of correlation versus PopVar
    # # Convert items for PopVar
    # # Pheno
    # pheno_use <- tp1$pheno_val$pheno_mean
    # geno_use <- genotype(genome1, combine_pop(list(tp1, par_pop)))
    # map_use <- genome2$gen_model$trait1 %>%
    #   select(qtl_name, chr, pos)
    #
    #
    # # # Pass to PopVar
    # # Convert genotypes into something useable for PopVar
    # geno_use <- as.data.frame(cbind( c("", row.names(geno_use)), rbind(colnames(geno_use), geno_use)) )
    #
    # pred_out_pv <- pop.predict(G.in = geno_use, y.in = pheno_use, map.in = map_use,
    #                           crossing.table = crossing_block, tail.p = 0.1, nInd = sim_pop_size,
    #                           min.maf = 0, mkr.cutoff = 1, entry.cutoff = 1, remove.dups = FALSE,
    #                           nSim = n_pops, nCV.iter = 1, models = "rrBLUP", impute = "pass")
    #
    # # Tidy
    # tidy_pred_out <- pred_out_pv$predictions %>%
    #   map(as.data.frame) %>%
    #   map(~mutate_all(., unlist) %>% rename_at(vars(contains("cor")), ~"pred.corG")) %>%
    #   list(., names(.)) %>%
    #   pmap_df(~mutate(.x, trait = str_extract(.y, "trait[0-9]{1}"))) %>%
    #   select(parent1 = Par1, parent2 = Par2, trait, pred.varG, pred.corG)
    #
    # ## Correlate predictions of genetic variance and correlation
    # pred_out %>%
    #   left_join(., tidy_pred_out) %>%
    #   distinct(parent1, parent2, pred_corG, pred.corG) %>%
    #   summarize(acc = cor(pred_corG, pred_corG))
    #
    # pred_out %>%
    #   left_join(., tidy_pred_out) %>%
    #   group_by(trait) %>%
    #   summarize(acc = cor(pred_varG, pred_varG))
    #
    # ## Completely accurate!




    ## Calculate the expected genetic variance in these populations
    expected_var <- calc_exp_genvar(genome = genome2, pedigree = ped, founder.pop = par_pop, crossing.block = crossing_block) %>%
      mutate(exp_musp = exp_mu + (k_sp * sqrt(exp_varG))) # %>%
      # # Add expectations of covariance
      # group_by(parent1, parent2) %>%
      # mutate(exp_covG = exp_corG[1] * (prod(sqrt(exp_varG)))) %>%
      # ungroup()
    
    ## Combine the expected and predicted results
    expected_predicted <- full_join(
      x = pred_out %>% rename_all(~str_replace(., "pred_", "")) %>% gather(parameter, prediction, -parent1:-trait),
      y = expected_var %>% rename_all(~str_replace(., "exp_", "")) %>% gather(parameter, expectation, -parent1:-trait),
      by = c("parent1", "parent2", "trait", "parameter")
    ) %>% 
      # If looking at the correlation, change the trait to simply trait1
      filter(!(trait == "trait2" & parameter %in% c("corG", "covG")))

    # Summarize
    results_summ <- expected_predicted %>% 
      group_by(trait, parameter) %>% 
      summarize(accuracy = cor(prediction, expectation, use = "complete.obs"), 
                bias = mean(prediction - expectation) / mean(expectation))

    ## Add the accuracy results to the results list
    results_out[[i]] <- list(
      summary = results_summ,
      metadata = list(tp_summ = tp_summ, mar_eff_meta = tp1$mar_eff_meta, 
                      pred_exp = subset(expected_predicted, trait == "trait1" & parameter == "corG", c(prediction, expectation)))
    )

  }

  # Add the results to the core_df, remove core
  core_df %>%
    mutate(results = results_out) %>%
    select(-core)

}, mc.cores = n_cores)

# Bind and save
popvar_prediction_simulation_out <- bind_rows(simulation_out)

# Save
if (check_results) {
  save_file <- file.path(result_dir, "popvar_gencor_simulation_prediction_results_missing.RData")
  save("popvar_prediction_simulation_out", file = save_file)
  
  
} else {
  save_file <- file.path(result_dir, "popvar_gencor_simulation_prediction_results.RData")
  save("popvar_prediction_simulation_out", file = save_file)
  
}













