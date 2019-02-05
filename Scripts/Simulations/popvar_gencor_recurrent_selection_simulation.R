## GenCorPrediction simulations of recurrent correlated response to selection
## by taking advantage of the predicted genetic correlation
## 
## 
## Author: Jeff Neyhart
## Last modified: September 26, 2018
## 

# Run the source script
repo_dir <- "/panfs/roc/groups/6/smithkp/neyha001/Genomic_Selection/GenCorPrediction/"
source(file.path(repo_dir, "source_MSI.R"))

# Load the two-row simulation genotypes
load(file.path(geno_dir, "s2_cap_simulation_data.RData"))

## Check if the results are present - if so only simulate the missing combinations
check_results <- F




# # Run on a local machine
# repo_dir <- getwd()
# source(file.path(repo_dir, "source.R"))
# # Additional libraries
# library(pbsim)
# 
# # Load the two-row simulation genotypes
# load(file.path(gdrive_dir, "BarleyLab/Projects/SideProjects/Resources/s2_cap_simulation_data.RData"))


# Filter for breeding programs relevant to my data
s2_cap_genos <- s2_cap_genos[str_detect(string = row.names(s2_cap_genos), pattern = "AB|BA|WA|N2|MT"),]
s2_cap_genos <- s2_cap_genos[,!colMeans(s2_cap_genos) %in% c(0, 2)]
s2_snp_info <- subset(s2_snp_info, rs %in% colnames(s2_cap_genos))




# Number of cores
n_cores <- 8 # Local machine
n_cores <- detectCores()








##### Recurrent selection experiment

## Fixed parameters
tp_size <- 600
tp_select <- 30
par_select <- 50
n_cross <- 20
i_sp <- 0.05
n_progeny <- (par_select / i_sp) / n_cross

L <- 100
n_iter <- 250
n_env <- 3
n_rep <- 1

# Selection intensity and number of cycles
# k_sp <- 1.76
k_sp <- 2.06
n_cycles <- 10


# Pedigree for later family development
ped <- sim_pedigree(n.ind = n_progeny, n.selfgen = Inf)


## Outline the parameters to perturb
trait1_h2_list <- 0.6
trait2_h2_list <- c(0.3, 0.6)
gencor_list <- c(-0.5, 0, 0.5)
selection_list <- c("mean", "muspC", "rand")
probcor_list <- data_frame(arch = c("pleio", "close_link", "loose_link"),
                           input = list(cbind(0, 1), cbind(5, 1), rbind(c(25, 0), c(35, 1)) ))

# Create a data.frame of parameters
param_df <- crossing(trait1_h2 = trait1_h2_list, trait2_h2 = trait2_h2_list, gencor = gencor_list, selection = selection_list,
                     probcor = probcor_list, iter = seq(n_iter))

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




## Check the results file
save_file <- file.path(result_dir, "popvar_gencor_recurrent_selection_simulation_results.RData")

# If it exists, load it and create the missing combinations
if (file.exists(save_file) & check_results) {
  load(save_file)
  
  ## Determine missing combinations
  missing_cases <- popvar_gencor_selection_simulation_out %>%
    select(-input, -results) %>%
    mutate_all(as.factor) %>%
    anti_join(x = complete_(., names(.)), y = .)  %>%
    mutate_all(as.character) %>%
    mutate_all(parse_guess)
  
  # Build a new parameter set
  param_df <- left_join(missing_cases, param_df)
  
}




# Split the parameter df
param_df_split <- param_df %>%
  assign_cores(n_cores) %>%
  split(.$core)


# Parallelize
simulation_out <- mclapply(X = param_df_split, FUN = function(core_df) {
  
  # # ## For local machine
  # i <- 1
  # core_df <- param_df_split[[i]]
  # # i = 3
  # ##
  
  # Create a results list
  results_out <- vector("list", nrow(core_df))
  
  
  # Iterate over the rows of the param_df
  for (i in seq_along(results_out)) {
    
    trait1_h2 <- core_df$trait1_h2[i]
    trait2_h2 <- core_df$trait2_h2[i]
    gencor <- core_df$gencor[i]
    selection <- core_df$selection[i]
    probcor <- core_df$input[[i]]
    
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
      # Predict marker effects
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
    
    ## Measure the frequency of favorable haplotypes (i.e. positive for trait1 and trait2)
    # First extract the allele effects and determine their sign
    qtl_pairs <- genome1$gen_model %>%
      map(~filter(., add_eff != 0))
    
    effect_pairs <- qtl_pairs %>%
      map("add_eff") %>% 
      pmap(c)
    # Get the pairs of the QTL
    qtl_pair_names <- qtl_pairs %>%
      map("qtl_name") %>% 
      pmap(c)
    
    fav_hap <- map(effect_pairs, ~sign(.) + 1) # These are the favorable haplotypes
    antag_hap1 <- map(fav_hap, ~c(2 - .[1], .[2])) # These are the antagonistic haplotypes
    antag_hap2 <- map(antag_hap1, ~2 - .) # These are the antagonistic haplotypes
    unfav_hap <- map(fav_hap, ~2 - .) # These are the unfavorable haplotypes
    
    ## List of haplotypes
    haplo_list <- list(favorable = fav_hap, antagonistic1 = antag_hap1, antagonistic2 = antag_hap2, unfavorable = unfav_hap)
    
    geno_mat <- do.call("cbind", tp1$geno)
    # Get a list of loci and the genotypes
    loci_geno_list <- map(qtl_pair_names, ~geno_mat[,.])
    
    ## Get the frequency of each haplotype
    tp_hap_freq <- haplo_list %>%
      map(~map2_dbl(.x = ., .y = loci_geno_list, ~mean(.y[,1] == .x[1] & .y[,2] == .x[2])))

    tp_haplotype_LD <- loci_geno_list %>% 
      # map(colnames) %>% map(~calc_LD(genome = genome1, pop = tp1, measure = "D", loci = .)) %>% map_dbl(~.[1,2])
      map_dbl(~cor(.)[1,2]) %>% ifelse(is.na(.), 0, .)
    
      
    # # Calculate the LD between these haplotypes
    # tp_LD <- cor(geno_mat[, unlist(qtl_pair_names)])[qtl_pairs[[1]]$qtl_name, qtl_pairs[[2]]$qtl_name]^2
    # tp_haplotype_LD <- 
    
    
    
    ## Start the recurrent selection
    recurrent_selection_out <- vector("list", n_cycles + 1)
    
    # Designate the tp as the first set of selection candidates
    candidates <- tp1

    # Iterate over cycles
    for (r in seq(n_cycles)) {
      
      par_pop_all <- pred_geno_val(genome = genome1, training.pop = tp1, candidate.pop = candidates) %>%
        {.$pred_val} %>% 
        mutate_at(vars(contains("trait")), funs(scale = as.numeric(scale(.)))) %>% 
        mutate(index = as.numeric(cbind(trait1_scale, trait2_scale) %*% weights))
      
      # Subset
      if (r == 1) {
        par_pop <- subset_pop(pop = candidates, individual = par_pop_all$ind[order(par_pop_all$index, decreasing = TRUE)[seq(tp_select)]])
      } else {
        par_pop <- subset_pop(pop = candidates, individual = par_pop_all$ind[order(par_pop_all$index, decreasing = TRUE)[seq(par_select)]])
      }
        
        
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
      
      # Calculate haplotype frequencies
      geno_mat <- do.call("cbind", par_pop$geno)
      # Get a list of loci and the genotypes
      loci_geno_list <- map(qtl_pair_names, ~geno_mat[,.])
      
      ## Get the frequency of each haplotype
      par_hap_freq <- haplo_list %>%
        map(~map2_dbl(.x = ., .y = loci_geno_list, ~mean(.y[,1] == .x[1] & .y[,2] == .x[2])))
      
      
      # Calculate the LD between these haplotypes
      par_haplotype_LD <- loci_geno_list %>% 
        # map(colnames) %>% map(~calc_LD(genome = genome1, pop = par_pop, measure = "D", loci = .)) %>% map_dbl(~.[1,2])
        map_dbl(~cor(.)[1,2]) %>% ifelse(is.na(.), 0, .)
      
      
      ## Create a crossing block with all possible crosses
      crossing_block_use <- sim_crossing_block(parents = indnames(par_pop), n.crosses = choose(nind(par_pop), 2))
      
      
      ## Split the stream based on parental selection
      # 1. mean - choose parents that maximize the predicted mean
      if (selection == "mean") {
        
        pred_out <- crossing_block_use %>%
          left_join(., par_pop_pgv, by = c("parent1" = "ind")) %>%
          left_join(., par_pop_pgv, by = c("parent2" = "ind", "trait")) %>%
          mutate(pred_mu = (pgv.x + pgv.y) / 2) %>% 
          group_by(trait) %>% 
          mutate(pred_mu = as.numeric(scale(pred_mu))) %>%
          group_by(parent1, parent2) %>%
          mutate(index = mean(pred_mu)) %>%
          select(-contains("pgv")) %>%
          ungroup()
        
        # Select on the index
        cb_select <- pred_out %>%
          filter(trait == "trait1") %>%
          arrange(desc(index)) %>%
          slice(1:n_cross) %>%
          select(contains("parent"))
        
        pred_variance <- NULL
        
      } else if (selection == "muspC") {
        
        ## Predict genetic variance and correlation
        pred_out <- pred_genvar(genome = genome1, pedigree = ped, training.pop = tp1, founder.pop = par_pop,
                                crossing.block = crossing_block_use) %>%
          mutate(pred_musp = pred_mu + (k_sp * sqrt(pred_varG))) %>%
          group_by(parent1, parent2) %>%
          mutate(pred_muspC = rev(pred_mu + (pred_corG * k_sp * sqrt(pred_varG)))) %>%
          group_by(trait) %>% 
          mutate_at(vars(contains("musp")), ~as.numeric(scale(.))) %>%
          ungroup() %>%
          mutate(index = (pred_musp + pred_muspC) / 2)
        
        # Select on the index
        cb_select <- pred_out %>%
          filter(trait == "trait1") %>%
          arrange(desc(index)) %>%
          slice(1:n_cross) %>%
          select(contains("parent"))
        
        ## Summarize variance of mean, sd, and corG
        pred_variance <- pred_out %>% 
          mutate(pred_sdG = sqrt(pred_varG)) %>% 
          group_by(trait) %>% 
          summarize_at(vars(pred_mu, pred_corG, pred_sdG), var)
        
      } else if (selection == "rand") {
        
        # Randomly select crosses
        cb_select <- sample_n(tbl = crossing_block_use, size = n_cross)
        
        pred_variance <- NULL
        
        
      }
      
      
      ## Create the crosses
      candidates <- sim_family_cb(genome = genome1, pedigree = ped, founder.pop = par_pop, crossing.block = cb_select, cycle.num = r) %>%
        pred_geno_val(genome = genome1, training.pop = tp1, candidate.pop = .)
      
      
      candidates_summ <- candidates$geno_val %>% 
        mutate(cor = cor(trait1, trait2)) %>% 
        gather(trait, value, trait1, trait2) %>% 
        group_by(trait) %>% 
        summarize_at(vars(cor, value), funs(mean, var)) %>% 
        select(trait, mean = value_mean, var = value_var, cor = cor_mean)
      
      
      # Calculate haplotype frequencies
      geno_mat <- do.call("cbind", candidates$geno)
      # Get a list of loci and the genotypes
      loci_geno_list <- map(qtl_pair_names, ~geno_mat[,.])
      ## Get the frequency of each haplotype
      cand_hap_freq <- haplo_list %>%
        map(~map2_dbl(.x = ., .y = loci_geno_list, ~mean(.y[,1] == .x[1] & .y[,2] == .x[2])))
      
      
      # Calculate the LD between these haplotypes
      cand_haplotype_LD <- loci_geno_list %>% 
        # map(colnames) %>% map(~calc_LD(genome = genome1, pop = candidates, measure = "D", loci = .)) %>% map_dbl(~.[1,2])
        map_dbl(~cor(.)[1,2]) %>% ifelse(is.na(.), 0, .)
      
      
      # Summarize the haplotype frequencies
      haplo_freq <- cbind(data.frame(population = c("parents", "candidates"), stringsAsFactors = FALSE), 
                          rbind(map_dbl(par_hap_freq, mean), map_dbl(cand_hap_freq, mean)),
                          hap_LD = c(mean(par_haplotype_LD, na.rm = T), mean(cand_haplotype_LD, na.rm = T)))
      
      ## Add the parent and candidate summary to the list
      recurrent_selection_out[[r]] <- list(cycle = r, parents = par_pop_summ, candidates = candidates_summ, haplo_freq = haplo_freq,
                                           pred_variance = pred_variance)
      
    }
    
    ## Make a final selection
    par_pop_all <- pred_geno_val(genome = genome1, training.pop = tp1, candidate.pop = candidates) %>%
    {.$pred_val} %>% 
      mutate_at(vars(contains("trait")), funs(scale = as.numeric(scale(.)))) %>% 
      mutate(index = as.numeric(cbind(trait1_scale, trait2_scale) %*% weights))
    
    # Subset
    par_pop <- subset_pop(pop = candidates, individual = par_pop_all$ind[order(par_pop_all$index, decreasing = TRUE)[seq(tp_select)]])
    
    
    # Measure the genetic variance, co-variance, and correlation in the selected tp
    par_pop_summ <- par_pop$geno_val %>%
      mutate(cor = cor(trait1, trait2)) %>% 
      gather(trait, value, trait1, trait2) %>% 
      group_by(trait) %>% 
      summarize_at(vars(cor, value), funs(mean, var)) %>% 
      select(trait, mean = value_mean, var = value_var, cor = cor_mean)
    
    # Calculate haplotype frequencies
    geno_mat <- do.call("cbind", par_pop$geno)
    # Get a list of loci and the genotypes
    loci_geno_list <- map(qtl_pair_names, ~geno_mat[,.])
    
    par_hap_freq <- haplo_list %>%
      map(~map2_dbl(.x = ., .y = loci_geno_list, ~mean(.y[,1] == .x[1] & .y[,2] == .x[2])))
    
    
    # Calculate the LD between these haplotypes
    par_haplotype_LD <- loci_geno_list %>% 
      # map(colnames) %>% map(~calc_LD(genome = genome1, pop = par_pop, measure = "D", loci = .)) %>% map_dbl(~.[1,2])
    map_dbl(~cor(.)[1,2]) %>% ifelse(is.na(.), 0, .)
    
    
    # Summarize the haplotype frequencies
    haplo_freq <- cbind(data.frame(population = c("parents"), stringsAsFactors = FALSE), 
                        rbind(map_dbl(par_hap_freq, mean)),
                        hap_LD = c(mean(par_haplotype_LD, na.rm = T)))
    
    
    ## Add the parent summary to the list
    recurrent_selection_out[[r + 1]] <- list(cycle = r + 1, parents = par_pop_summ, haplo_freq = haplo_freq)
    
    
    # Tidy
    recurrent_selection_out1 <- recurrent_selection_out %>%
      map_df(~t(.) %>% as_data_frame) %>% 
      unnest(cycle)
  
    haplo_freq_tidy <- recurrent_selection_out1 %>% 
      select(cycle, haplo_freq) %>% 
      unnest()
    
    
    
    ## Tidy everything
    recurrent_selection_tidy <- recurrent_selection_out1 %>% 
      select(-haplo_freq) %>%
      gather(population, summary, -cycle) %>% 
      filter(!map_lgl(summary, is.null)) %>% 
      unnest() %>%
      full_join(., haplo_freq_tidy, by = c("cycle", "population"))
    
    ## Add the response and other results
    results_out[[i]] <- bind_rows(
      cbind(tp_summ, cycle = 0, population = "tp", t(map_dbl(tp_hap_freq, mean)),
                                             hap_LD = mean(tp_haplotype_LD, na.rm = T)), 
      recurrent_selection_tidy) %>%
      select(cycle, trait, population, mean, var, cor, names(.))
    
  }
  
  # Add the results to the core_df, remove core
  core_df %>%
    mutate(results = results_out) %>%
    select(-core)
  
}, mc.cores = n_cores)

# Bind and save
popvar_gencor_selection_simulation_out <- bind_rows(simulation_out)



# Save as separate file for missing data if check_results == T
if (check_results) {
  save_file <- file.path(result_dir, "popvar_gencor_recurrent_selection_simulation_results_missing.RData")
  save("popvar_gencor_selection_simulation_out", file = save_file)
  
  
} else {
  save_file <- file.path(result_dir, "popvar_gencor_recurrent_selection_simulation_results.RData")
  save("popvar_gencor_selection_simulation_out", file = save_file)
  
}



















