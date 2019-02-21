## PopVarVal simulations of correlated response to selection
## by taking advantage of the predicted genetic correlation
## 
## 
## Author: Jeff Neyhart
## Last modified: September 26, 2018
## 

# Run the source script
repo_dir <- "/panfs/roc/groups/6/smithkp/neyha001/Genomic_Selection/GenCorPrediction/"
source(file.path(repo_dir, "source_MSI.R"))



# # Run on a local machine
# repo_dir <- getwd()
# source(file.path(repo_dir, "source.R"))
# # Additional libraries
# library(pbsim)
# library(pbsimData)
# 
# # Load the two-row simulation genotypes
# load(file.path(gdrive_dir, "BarleyLab/Projects/SideProjects/Resources/s2_cap_simulation_data.RData"))



# Number of cores
n_cores <- 8 # Local machine
n_cores <- detectCores()



## Fixed parameters
tp_size <- 600
tp_select <- 30 # Number of TP individuals to choose as potential parents
nCandidates <- 1200 # Total number of progeny to generate in the next cycle
L <- 100
n_iter <- 50
n_env <- 3
n_rep <- 1
k_sp <- 1.75

## Mutable parameters
# Selection intensities
ints <- seq(0.01, 0.20, by = 0.01)
# Number of populations
nPop <- c(5, 20, 60, 120)
nPop <- setNames(nPop, paste0("nPop", nPop))
# Size of populations
popSize <- nCandidates / nPop
# Trait heritabilities
trait1_h2_list <- c(0.6)
trait2_h2_list <- c(0.6, 0.3)




## Outline the parameters to perturb
gencor_list <- c(-0.5, 0.5)
probcor_list <- data_frame(arch = c("pleio", "close_link", "loose_link"),
                           input = list(cbind(0, 1), cbind(5, 1), rbind(c(25, 0), c(30, 1)) ))

# Create a data.frame of parameters
param_df <- crossing(trait1_h2 = trait1_h2_list, trait2_h2 = trait2_h2_list, gencor = gencor_list, probcor = probcor_list, iter = seq(n_iter))

map_sim <- s2_snp_info %>%
  split(.$chrom) %>%
  map(~setNames(.$cM_pos, .$rs)) %>%
  map(~structure(., class = "A"))  %>%
  structure(., class = "map") %>%
  # Jitter
  qtl::jittermap(.) %>%
  `names<-`(., seq_along(.))

# Create the base genome - this is fixed for all simulations
genome <- sim_genome(map = map_sim)









# Split the parameter df
param_df_split <- param_df %>%
  assign_cores(n_cores) %>%
  split(.$core)

# Parallelize
simulation_out <- mclapply(X = param_df_split, FUN = function(core_df) {

  # ## For local machine
  # i <- 1
  # core_df <- param_df_split[[i]]
  # ##

  # Create a results list
  results_out <- vector("list", nrow(core_df))

  # Iterate over the rows of the param_df
  for (i in seq_along(results_out)) {

    gencor <- core_df$gencor[i]
    probcor <- core_df$input[[i]]
    trait1_h2 <- core_df$trait1_h2[i]
    trait2_h2 <- core_df$trait2_h2[i]
    
    # Set the maximum number of QTL
    max_qtl <- ifelse(probcor[1] != 0, L, 1.5 * L)

    # Simulate QTL
    qtl_model <- replicate(n = 2, matrix(NA, ncol = 4, nrow = L), simplify = FALSE)
    genome1 <- sim_multi_gen_model(genome = genome, qtl.model = qtl_model, add.dist = "geometric", max.qtl = max_qtl,
                                   corr = gencor, prob.corr = probcor)
    
    

    ## Adjust the genetic architecture - only if pleiotropy is not present
    if (probcor[1] != 0) {
      genome1 <- adj_multi_gen_model(genome = genome1, geno = s2_cap_genos, gencor = gencor)
    }

    
    # Create the TP by random selection
    tp1 <- create_pop(genome = genome1, geno = s2_cap_genos[sort(sample(nrow(s2_cap_genos), size = tp_size)),]) %>%
      # Phenotype the base population
      sim_phenoval(pop = ., h2 = c(trait1_h2, trait2_h2), n.env = n_env, n.rep = n_rep) %>%
      # Add marker effects
      pred_mar_eff(genome = genome1, training.pop = ., method = "RRBLUP")

    # Measure the genetic variance, genetic covariance, and genetic correlation in the tp1
    tp_summ <- tp1$geno_val %>%
      mutate(cor = cor(trait1, trait2)) %>% 
      gather(trait, value, trait1, trait2) %>% 
      group_by(trait) %>% 
      summarize_at(vars(cor, value), funs(mean, var)) %>% 
      select(trait, mean = value_mean, var = value_var, cor = cor_mean)
    
    ## Predict the genotypic values for the TP, then select the best from the tp
    # First create a set of weights based on the observed phenotypic variance
    # weights <- 1 / sqrt(diag(var(tp1$pheno_val$pheno_mean[,-1])))
    weights <- c(0.5, 0.5)

    par_pop_all <- pred_geno_val(genome = genome1, training.pop = tp1, candidate.pop = tp1) %>%
      {.$pred_val} %>% 
      mutate_at(vars(contains("trait")), funs(scale = as.numeric(scale(.)))) %>% 
      mutate(index = as.numeric(cbind(trait1_scale, trait2_scale) %*% weights))
    
    # Subset
    par_pop <- subset_pop(pop = tp1, individual = par_pop_all$ind[order(par_pop_all$index, decreasing = TRUE)[seq(tp_select)]])
    

    # Get the PGVs
    par_pop_pgv <- par_pop_all %>%
      select(ind, trait1, trait2) %>%
      gather(trait, pgv, -ind) %>%
      mutate(ind = as.character(ind))
    
    
    # Measure the genetic variance, co-variance, and correlation in the selected tp
    par_pop_summ <- par_pop$geno_val %>%
      mutate(cor = cor(trait1, trait2)) %>% 
      gather(trait, value, trait1, trait2) %>% 
      group_by(trait) %>% 
      summarize_at(vars(cor, value), funs(mean, var)) %>% 
      select(trait, mean = value_mean, var = value_var, cor = cor_mean)

    
    ## Create a crossing block with all of the possible non-reciprocal combinations of the TP
    crossing_block <- sim_crossing_block(parents = indnames(par_pop), n.crosses = choose(nind(par_pop), 2))

    # Select the best
    # crossing_block_use <- crossing_block_pgv[order(crossing_block_pgv$index, decreasing = T)[1:cross_preselect],]
    crossing_block_use <- crossing_block

    
    ## Predict genetic variance and correlation
    pred_out <- pred_genvar(genome = genome1, pedigree = sim_pedigree(n.ind = 5), training.pop = tp1, founder.pop = par_pop,
                            crossing.block = select(crossing_block_use, contains("parent"))) %>%
      mutate(pred_musp = pred_mu + (k_sp * sqrt(pred_varG))) %>%
      group_by(parent1, parent2) %>%
      mutate(pred_muspC = rev(pred_mu + (pred_corG * k_sp * sqrt(pred_varG)))) %>%
      group_by(trait) %>% 
      mutate_at(vars(contains("mu")), ~as.numeric(scale(.))) %>% # Create indices for both mu and mu_sp
      group_by(parent1, parent2) %>%
      mutate(mu_index = mean(pred_mu)) %>%
      ungroup() %>%
      mutate(mu_sp_index = (pred_musp + pred_muspC) / 2)


    ### For each selection method, select the top N crosses, where N is equal to nPop[i]
    ## Select the crosses that maximize the correlated response of trait 2
    cb_select_muspC <- pred_out %>%
      filter(trait == "trait1") %>%
      arrange(desc(mu_sp_index)) %>%
      {map(.x = nPop, function(n) slice(., 1:n) %>% select(contains("parent")) )}

    ## Now select crosses with the most favorable correlation
    cb_select_corG <- pred_out %>%
      filter(trait == "trait1") %>%
      arrange(desc(pred_corG)) %>%
      {map(.x = nPop, function(n) slice(., 1:n) %>% select(contains("parent")) )}

    ## Now select crosses based on an index of the predicted mean PGV
    cb_select_mean <- pred_out %>%
      filter(trait == "trait1") %>%
      arrange(desc(mu_index)) %>%
      {map(.x = nPop, function(n) slice(., 1:n) %>% select(contains("parent")) )}

    ## Now select random crosses
    cb_select_rand <- pred_out %>%
      select(contains("parent")) %>%
      {map(.x = nPop, function(n) sample_n(., size = n) )}
      

    ## For each crossing block, make the crosses
    cycle1_list <- list(muspC = cb_select_muspC, corG = cb_select_corG, mean = cb_select_mean, rand = cb_select_rand) %>%
      map(., ~map(., ~{
        cb <- .
        ped <- sim_pedigree(n.ind = nCandidates / nrow(cb))
        sim_family_cb(genome = genome1, pedigree = ped, founder.pop = par_pop, crossing.block = cb) }) )
    
    
    ## For each population, measure the mean, variance, and correlation
    cycle1_list_summ <- cycle1_list %>%
      map(., ~map(., "geno_val") %>%
            map2_df(.x = ., .y = names(.), ~mutate(.x, cor = cor(trait1, trait2)) %>% 
                   gather(trait, value, trait1, trait2) %>% 
                   group_by(trait) %>% 
                   summarize_at(vars(cor, value), funs(mean, var)) %>% 
                   mutate(nPop = .y) %>%
                   select(nPop, trait, mean = value_mean, var = value_var, cor = cor_mean) ) ) %>%
      map2_df(.x = ., .y = names(.), ~mutate(.x, selection = .y))
    
    ## For each population, predict genotypic values and create an index for selection
    ## Then use the index to select members of the population based on varying intensities
    ## Then calculate summaries
    cycle1_list_select <- cycle1_list %>%
      map(., ~map2_df(.x = ., .y = names(.), ~{
        cand_pop <- .x
        # Predict, then calculate index of PGVs
        pop1 <- pred_geno_val(genome = genome1, training.pop = tp1, candidate.pop = cand_pop, method = "RRBLUP")
        # Calculate an index
        pop_index <- pop1$pred_val %>%
          mutate_at(vars(contains("trait")), funs(scale = as.numeric(scale(.)))) %>% 
          mutate(index = as.numeric(cbind(trait1_scale, trait2_scale) %*% weights))
        
        # Select from the whole population
        pop_index_select_all <- pop_index %>%
          arrange(desc(index)) %>% 
          {map(ints, function(isp) slice(., 1:(isp * nrow(.))))}
        
        # # Select from each family
        # pop_index_select_family <- pop_index %>% 
        #   mutate(family = str_extract(ind, "[0-9]{4}")) %>% 
        #   group_by(family) %>% 
        #   arrange(desc(index)) %>% 
        #   {map(ints, function(isp) slice(., 1:(isp * nrow(pop_index))))} %>%
        #   map(ungroup)
        # 
        
        
        
        # Use the different selection intensities to subset the population
        # Then extract the genotypic values
        pop_index_select_all %>% 
          map("ind") %>% 
          map(as.character) %>%
          map(~subset_pop(pop = pop1, individual = .)) %>%
          map("geno_val") %>%
          map(., ~mutate(., cor = cor(trait1, trait2)) %>% 
                gather(trait, value, trait1, trait2) %>% 
                group_by(trait) %>% 
                summarize_at(vars(cor, value), funs(mean, var)) %>% 
                select(trait, mean = value_mean, var = value_var, cor = cor_mean) ) %>%
          map2_df(.x = ., .y = ints, ~mutate(.x, intensity = .y)) %>%
          mutate(nPop = .y)
        
    
      }) ) %>% map2_df(.x = ., .y = names(.), ~mutate(.x, selection = .y))
           
    
    ## Convert mean to response units (base population units)
    cycle1_list_response <- cycle1_list_select %>%
      left_join(., rename_at(tp_summ, vars(-trait), funs(paste0("base_", .))), by = "trait") %>%
      mutate(response = (mean - base_mean) / sqrt(base_var),
             stand_sd = sqrt(var) / sqrt(base_var)) %>%
      select(selection, nPop, intensity, trait, mean, var, cor, response, stand_sd)


    ## Add the response and other results
    results_out[[i]] <- list(
      response = cycle1_list_response,
      metadata = list(tp_summ = tp_summ, par_summ = par_pop_summ, cand_summ = cycle1_list_summ)
    )

  }

  # Add the results to the core_df, remove core
  core_df %>%
    mutate(results = results_out) %>%
    select(-core)

}, mc.cores = n_cores)

# Bind and save
popvar_gencor_cycle1_selection_simulation_out <- bind_rows(simulation_out)

# Save
save_file <- file.path(result_dir, "popvar_gencor_selection_simulation_results.RData")
save("popvar_gencor_cycle1_selection_simulation_out", file = save_file)

