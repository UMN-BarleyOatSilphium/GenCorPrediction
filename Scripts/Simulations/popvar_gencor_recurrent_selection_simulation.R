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
tp_select <- 25
cross_select <- 25
i_sp <- 0.05
n_progeny <- (cross_select / i_sp) / cross_select

L <- 100
n_iter <- 100
n_env <- 3
n_rep <- 1

# Selection intensity and number of cycles
# k_sp <- 1.76
k_sp <- 2.06
n_cycles <- 10


## Outline the parameters to perturb
trait1_h2_list <- 0.6
trait2_h2_list <- c(0.3, 0.6)
gencor_list <- c(-0.5, 0.5)
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



# Split the parameter df
param_df_split <- param_df %>%
  assign_cores(n_cores) %>%
  split(.$core)

# Pedigree for later family development
ped <- sim_pedigree(n.ind = n_progeny, n.selfgen = Inf)

# Parallelize
simulation_out <- mclapply(X = param_df_split, FUN = function(core_df) {
  
  # ## For local machine
  # i <- 1
  # core_df <- param_df_split[[i]]
  # i <- 80
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
    
    # Simulate QTL
    qtl_model <- replicate(n = 2, matrix(NA, ncol = 4, nrow = L), simplify = FALSE)
    genome1 <- sim_multi_gen_model(genome = genome, qtl.model = qtl_model, add.dist = "geometric", max.qtl = L,
                                   corr = gencor, prob.corr = probcor)
    
    ## Adjust the genetic architecture - only if pleiotropy is not present
    if (probcor[1] != 0) {
      genome1 <- adj_multi_gen_model(genome = genome1, geno = s2_cap_genos, gencor = gencor)
    }
    
    # Create the TP by random selection
    tp1 <- create_pop(genome = genome1, geno = s2_cap_genos[sort(sample(nrow(s2_cap_genos), size = tp_size)),]) %>%
      # Phenotype the base population
      sim_phenoval(pop = ., h2 = c(trait1_h2, trait2_h2), n.env = n_env, n.rep = n_rep)
    
    # Measure the genetic variance, genetic covariance, and genetic correlation in the tp1
    tp_summ <- tp1$geno_val %>%
      mutate(cor = cor(trait1, trait2)) %>% 
      gather(trait, value, trait1, trait2) %>% 
      group_by(trait) %>% 
      summarize_at(vars(cor, value), funs(mean, var)) %>% 
      select(trait, mean = value_mean, var = value_var, cor = cor_mean)
    
    ## Predict the genotypic values for the TP, then select the best from the tp
    # First create a set of weights based on the observed phenotypic variance
    weights <- c(0.5, 0.5)
    
    ## Measure the frequency of favorable haplotypes (i.e. positive for trait1 and trait2)
    # First extract the allele effects and determine their sign
    effect_pairs <- map(genome1$gen_model, "add_eff") %>% pmap(c)
    fav_hap <- map(effect_pairs, ~sign(.) + 1) # These are the favorable haplotypes
    neg_hap <- map(fav_hap, ~2 - .) # These are the negative haplotypes
    
    geno_mat <- do.call("cbind", tp1$geno)
    # Get a list of loci and the genotypes
    loci_geno_list <- map(genome1$gen_model, "qtl_name") %>% pmap(c) %>% map(~geno_mat[,.])
    # Get the frequency of the favorable haplotype
    tp_fav_hap_freq <- list(fav_hap, loci_geno_list) %>% pmap_dbl(~mean(.y[,1] == .x[1] & .y[,2] == .x[2]))
    tp_neg_hap_freq <- list(neg_hap, loci_geno_list) %>% pmap_dbl(~mean(.y[,1] == .x[1] & .y[,2] == .x[2]))
    
    ## Start the recurrent selection
    recurrent_selection_out <- vector("list", n_cycles + 1)
    
    # Designate the tp as the first set of selection candidates
    candidates <- tp1
    
    # Iterate over cycles
    for (r in seq(n_cycles)) {
      
      par_pop <- pred_geno_val(genome = genome1, training.pop = tp1, candidate.pop = candidates) %>%
        select_pop(pop = ., intensity = tp_select, index = weights, type = "genomic")
      
      # Get the PGVs
      par_pop_pgv <- par_pop$pred_val %>%
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
      loci_geno_list <- map(genome1$gen_model, "qtl_name") %>% pmap(c) %>% map(~geno_mat[,.])
      # Get the frequency of the favorable haplotype
      par_fav_hap_freq <- list(fav_hap, loci_geno_list) %>% pmap_dbl(~mean(.y[,1] == .x[1] & .y[,2] == .x[2]))
      par_neg_hap_freq <- list(neg_hap, loci_geno_list) %>% pmap_dbl(~mean(.y[,1] == .x[1] & .y[,2] == .x[2]))
      
      ## Create a crossing block with all of the possible non-reciprocal combinations of the TP
      crossing_block_use <- sim_crossing_block(parents = indnames(par_pop), n.crosses = choose(nind(par_pop), 2))
      
      
      ## Split the stream based on parental selection
      # 1. mean - choose parents that maximize the predicted mean
      if (selection == "mean") {
        
        pred_out <- crossing_block_use %>%
          left_join(., par_pop_pgv, by = c("parent1" = "ind")) %>%
          left_join(., par_pop_pgv, by = c("parent2" = "ind", "trait")) %>%
          mutate(pred_mu = (pgv.x + pgv.y) / 2) %>% 
          group_by(parent1, parent2) %>%
          mutate(index = mean(pred_mu)) %>%
          select(-contains("pgv")) %>%
          ungroup()
        
        # Select on the index
        cb_select <- pred_out %>%
          filter(trait == "trait1") %>%
          arrange(desc(index)) %>%
          slice(1:cross_select) %>%
          select(contains("parent"))
        
      } else if (selection == "muspC") {
        
        ## Predict genetic variance and correlation
        pred_out <- pred_genvar(genome = genome1, pedigree = ped, training.pop = tp1, founder.pop = par_pop,
                                crossing.block = crossing_block_use) %>%
          mutate(pred_musp = pred_mu + (k_sp * sqrt(pred_varG))) %>%
          group_by(parent1, parent2) %>%
          mutate(pred_muspC = rev(pred_mu + (pred_corG * k_sp * sqrt(pred_varG)))) %>%
          ungroup() %>%
          mutate(index = (pred_musp * weights[1]) + (pred_muspC * weights[2])) 
        
        # Select on the index
        cb_select <- pred_out %>%
          filter(trait == "trait1") %>%
          arrange(desc(index)) %>%
          slice(1:cross_select) %>%
          select(contains("parent"))
        
      } else if (selection == "rand") {
        
        # Randomly select crosses
        cb_select <- sample_n(tbl = crossing_block_use, size = cross_select)
        
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
      loci_geno_list <- map(genome1$gen_model, "qtl_name") %>% pmap(c) %>% map(~geno_mat[,.])
      # Get the frequency of the favorable haplotype
      cand_fav_hap_freq <- list(fav_hap, loci_geno_list) %>% pmap_dbl(~mean(.y[,1] == .x[1] & .y[,2] == .x[2]))
      cand_neg_hap_freq <- list(neg_hap, loci_geno_list) %>% pmap_dbl(~mean(.y[,1] == .x[1] & .y[,2] == .x[2]))
      
      # Summarize the haplotype frequencies
      haplo_freq <- data.frame(population = c("parents", "candidates"), pos_freq = c(mean(par_fav_hap_freq), mean(cand_fav_hap_freq)), 
                               neg_freq = c(mean(par_neg_hap_freq), mean(cand_neg_hap_freq)), stringsAsFactors = FALSE)
      
      ## Add the parent and candidate summary to the list
      recurrent_selection_out[[r]] <- list(cycle = r, parents = par_pop_summ, candidates = candidates_summ, haplo_freq = haplo_freq)
      
    }
    
    ## Make a final selection
    par_pop <- pred_geno_val(genome = genome1, training.pop = tp1, candidate.pop = candidates) %>%
      select_pop(pop = ., intensity = tp_select, index = weights, type = "genomic")
    
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
    loci_geno_list <- map(genome1$gen_model, "qtl_name") %>% pmap(c) %>% map(~geno_mat[,.])
    # Get the frequency of the favorable haplotype
    par_fav_hap_freq <- list(fav_hap, loci_geno_list) %>% pmap_dbl(~mean(.y[,1] == .x[1] & .y[,2] == .x[2]))
    par_neg_hap_freq <- list(neg_hap, loci_geno_list) %>% pmap_dbl(~mean(.y[,1] == .x[1] & .y[,2] == .x[2]))
    
    haplo_freq <- data.frame(population = c("parents"), pos_freq = c(mean(par_fav_hap_freq)), neg_freq = c(mean(par_neg_hap_freq)), stringsAsFactors = FALSE)
    
    
    ## Add the parent summary to the list
    recurrent_selection_out[[r + 1]] <- list(cycle = r + 1, parents = par_pop_summ, haplo_freq = haplo_freq)
    
    
    # Tidy
    recurrent_selection_out1 <- recurrent_selection_out %>%
      map_df(~t(.) %>% as_data_frame) %>% 
      unnest(cycle)
  
    haplo_freq_tidy <- recurrent_selection_out1 %>% 
      select(cycle, haplo_freq) %>% 
      unnest()
    
    # Adjust if the architecture is pleiotropic
    if (probcor[,1] == 0) {
      haplo_freq_tidy <- haplo_freq_tidy %>% mutate(neg_freq = 1 - pos_freq)
      tp_neg_hap_freq <- 1 - tp_fav_hap_freq
    }
    
    recurrent_selection_tidy <-recurrent_selection_out1 %>% 
      select(-haplo_freq) %>%
      gather(population, summary, -cycle) %>% 
      filter(!map_lgl(summary, is.null)) %>% 
      unnest() %>%
      full_join(., haplo_freq_tidy, by = c("cycle", "population"))
    
    ## Add the response and other results
    results_out[[i]] <- bind_rows(add_column(tp_summ, cycle = 0, population = "tp", pos_freq = mean(tp_fav_hap_freq), neg_freq = mean(tp_neg_hap_freq), .before = "trait"), 
                                  recurrent_selection_tidy)
    
  }
  
  # Add the results to the core_df, remove core
  core_df %>%
    mutate(results = results_out) %>%
    select(-core)
  
}, mc.cores = n_cores)

# Bind and save
popvar_gencor_selection_simulation_out <- bind_rows(simulation_out)

# Save
save_file <- file.path(result_dir, "popvar_gencor_recurrent_selection_simulation_results.RData")
save("popvar_gencor_selection_simulation_out", file = save_file)

















