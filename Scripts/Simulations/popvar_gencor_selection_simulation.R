## PopVarVal simulations of correlated response to selection
## by taking advantage of the predicted genetic correlation
## 
## 
## Author: Jeff Neyhart
## Last modified: September 7, 2018
## 

# Run the source script
repo_dir <- "/panfs/roc/groups/6/smithkp/neyha001/Genomic_Selection/PopVarVal/"
source(file.path(repo_dir, "source_MSI.R"))

# Load the two-row simulation genotypes
load(file.path(geno_dir, "s2_cap_simulation_data.RData"))




# # Run on a local machine
# repo_dir <- getwd()
# source(file.path(repo_dir, "source.R"))
# # Additional libraries
# library(pbsim)
# library(PopVar)
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



## Fixed parameters
tp_size <- 600
tp_select <- 25
# cross_preselect <- 150
cross_select <- 25
n_progeny <- 100
k_sp <- 1.76
k_sp1 <- 1.12

L <- 100
n_iter <- 50
n_env <- 3
n_rep <- 1
# Selection intensities
ints <- seq(0.01, 0.20, by = 0.01)


## Outline the parameters to perturb
trait1_h2_list <- trait2_h2_list <- c(0.3, 0.6, 1)
gencor_list <- c(-0.5, 0, 0.5)
probcor_list <- data_frame(arch = c("pleio", "close_link", "loose_link"),
                           input = list(cbind(0, 1), cbind(5, 1), cbind(30, 1) ))

# Create a data.frame of parameters
param_df <- crossing(trait1_h2 = trait1_h2_list, trait2_h2 = trait2_h2_list, gencor = gencor_list,
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

# Parallelize
simulation_out <- mclapply(X = param_df_split, FUN = function(core_df) {

  # ## For local machine
  # i <- 1
  # core_df <- param_df_split[[i]]
  # i <- 50
  # ##

  # Create a results list
  results_out <- vector("list", nrow(core_df))

  # Iterate over the rows of the param_df
  for (i in seq_along(results_out)) {

    trait1_h2 <- core_df$trait1_h2[i]
    trait2_h2 <- core_df$trait2_h2[i]
    gencor <- core_df$gencor[i]
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
    tp1_vcov <- var(tp1$geno_val[,-1])
    tp1_genvar <- diag(tp1_vcov)
    tp1_gencovar <- tp1_vcov[1,2]
    tp_cor <- cor(tp1$geno_val[,-1])[1,2]
    
    # Measure the phenotypic correlation in the TP
    tp_pheno_cor <- cor(tp1$pheno_val$pheno_mean[,-1])[1,2]
 
    ## Predict the genotypic values for the TP, then select the best from the tp
    # First create a set of weights based on the observed phenotypic variance
    weights <- c(0.5, 0.5)
    
    par_pop <- pred_geno_val(genome = genome1, training.pop = tp1, candidate.pop = tp1) %>%
      # select_pop(pop = ., intensity = sqrt(tp_select/nind(tp1)), index = c(1, 0), type = "genomic") %>%
      # select_pop(pop = ., intensity = sqrt(tp_select/nind(tp1)), index = c(0, 1), type = "genomic")
      select_pop(pop = ., intensity = tp_select, index = weights, type = "genomic")
    
    # Get the PGVs
    par_pop_pgv <- par_pop$pred_val %>%
      gather(trait, pgv, -ind) %>%
      mutate(ind = as.character(ind))

    # ## Plot the GVs of the tp and selections
    # tp1$geno_val %>%
    #   mutate(selected = ind %in% indnames(par_pop)) %>%
    #   qplot(x = trait1, y = trait2, data = ., color = selected)



    # Measure the genetic variance, co-variance, and correlation in the selected tp
    tp_select_vcov <- var(par_pop$geno_val[,-1])
    tp_select_var <- diag(tp_select_vcov)
    tp_select_covar <- tp_select_vcov[1,2]
    tp_select_cor <- cor(par_pop$geno_val[,-1])[1,2]

    ## Create a crossing block with all of the possible non-reciprocal combinations of the TP
    crossing_block <- sim_crossing_block(parents = indnames(par_pop), n.crosses = choose(nind(par_pop), 2))

    # Select the best
    # crossing_block_use <- crossing_block_pgv[order(crossing_block_pgv$index, decreasing = T)[1:cross_preselect],]
    crossing_block_use <- crossing_block
    
    # Pedigree for later family development
    ped <- sim_pedigree(n.ind = n_progeny, n.selfgen = Inf)


    ## Predict genetic variance and correlation
    pred_out <- pred_genvar(genome = genome1, pedigree = ped, training.pop = tp1, founder.pop = par_pop,
                            crossing.block = select(crossing_block_use, contains("parent"))) %>%
      mutate(pred_musp = pred_mu + (k_sp * sqrt(pred_varG))) %>%
      group_by(parent1, parent2) %>%
      # mutate(pred_muspC = rev(pred_mu + (k_sp * pred_corG * pred_corG))) %>%
      mutate(pred_muspC = rev(pred_mu + (k_sp1 * pred_corG * pred_corG))) %>%
      ungroup()



    # ## Get the expected correlation
    # expected_out <- calc_exp_genvar(genome = genome1, pedigree = ped, founder.pop = par_pop,
    #                                 crossing.block = select(crossing_block_use, contains("parent"))) %>%
    #   mutate(exp_musp = exp_mu + (pred_ksp * sqrt(exp_varG)))


    # Weight the mu_sp of trait1 and correlated mu_sp of trait2 to select crosses
    
    

    ## Select the crosses that maximize the sum of the correlated responses
    ## Alternatively, select the cross that gives the largest sum of the trait1 mu_sp and
    ## predicted correlated reponse in trait2
    cb_select_muspC <- pred_out %>% 
      filter(trait == "trait1") %>% 
      mutate(index = (pred_musp * weights[1]) + (pred_muspC * weights[2])) %>% 
      arrange(desc(index)) %>% 
      slice(1:cross_select) %>%
      select(contains("parent"))
    
    ## Now select crosses with the most favorable correlation
    cb_select_corG <- pred_out %>% 
      filter(trait == "trait1") %>% 
      arrange(desc(pred_corG)) %>% 
      slice(1:cross_select) %>%
      select(contains("parent"))

    ## Now select crosses based on an index of the predicted mean PGV
    cb_select_mean <- pred_out %>% 
      select(parent1:trait, pred_mu) %>% 
      spread(trait, pred_mu) %>% 
      mutate(mu_index = (trait1 * weights[1]) + (trait2 * weights[2])) %>%
      arrange(desc(mu_index)) %>% 
      slice(1:cross_select) %>%
      select(contains("parent"))

    # ## Now select on an index of the superior progeny means (this effectively ignores the correlation)
    # ## Select the crosses that maximize the correlated superior progeny mean of trait2
    # cb_select_musp <- pred_out %>% 
    #   select(parent1:trait, pred_musp) %>% 
    #   spread(trait, pred_musp) %>% 
    #   mutate(musp_index = (trait1 * weights[1]) + (trait2 * weights[2])) %>%
    #   arrange(desc(musp_index)) %>% 
    #   slice(1:cross_select) %>%
    #   select(contains("parent"))

    ## Now select random crosses
    cb_select_rand <- pred_out %>%
      select(contains("parent")) %>%
      sample_n(cross_select)


    ## Combine crossing blocks
    cb_list <- list(muspC = cb_select_muspC, corG = cb_select_corG, mean = cb_select_mean, rand = cb_select_rand)

    ### Create all of the crosses
    cycle1_list <- cb_list %>%
      map(~sim_family_cb(genome = genome1, pedigree = ped, founder.pop = par_pop, crossing.block = .) %>%
            pred_geno_val(genome = genome1, training.pop = tp1, candidate.pop = .))

    
    # Measure the genetic variance, covariance, and correlation in the C1
    c1_vcov <- cycle1_list %>%
      map("geno_val") %>% 
      map(~var(.[,-1]))
    
    c1_var <- map(c1_vcov, ~diag(.))
    c1_covar <- map(c1_vcov, ~.[1,2])
    c1_cor <- map(c1_vcov, ~.[1,2] / prod(sqrt(diag(.))))
    

    ## Select progeny on a range of selection intensities
    selections <- data_frame(intensity = ints)
    selections$summ <- selections$intensity %>% map(~{
      int <- .
      cycle1_list %>%
        map(~select_pop(pop = ., intensity = int, index = weights, type = "genomic")) %>%
        map("geno_val")
    })
  
    
    ## Summarize the mean and cor of each/both traits
    response <- selections %>%
      mutate(muspC = map(summ, "muspC"), corG = map(summ, "corG"), mean = map(summ, "mean"), rand = map(summ, "rand")) %>%
      select(-summ) %>% 
      gather(selection, out, -intensity) %>%
      unnest() %>%
      group_by(intensity, selection) %>% 
      mutate(corG = cor(trait1, trait2)) %>% 
      gather(trait, value, trait1, trait2) %>% 
      group_by(intensity, selection, trait, corG) %>% 
      summarize_at(vars(value), funs(mean, var)) %>%
      ungroup()
    
    ## Summarize the response
    response_summary <- response %>%
      # First add the mean from the "mean" selection
      left_join(., rename(subset(response, selection == "mean", c(intensity, trait, mean)), mean_control = mean), by = c("intensity", "trait")) %>%
      # Add the genetic variance from the base population
      left_join(., gather(as.data.frame(t(tp1_genvar)), trait, base_sdG), by = "trait") %>%
      # Add the mean from the parent population
      left_join(., gather(summarize_at(par_pop$geno_val, vars(contains("trait")), mean), trait, par_mean), by  = "trait") %>%
      # Calculate the standardized response (mean of progeny minus mean of parents)
      # Calcualte the relative response (mean of progeny minus mean of progeny under different selection method)
      mutate(relative_response = (mean - mean_control) / base_sdG,
             stand_response = (mean - par_mean) / base_sdG) %>%
      # Calculate an index for the responses
      group_by(intensity, selection) %>%
      mutate_at(vars(contains("response")), funs(index = sum(. * weights))) %>%
      ungroup() %>%
      select(-mean_control:-par_mean)
      
    ## How many of the crosses are the same?
    common_crosses <- combn(x = cb_list, m = 2, simplify = FALSE) %>%
      map_df(~data.frame(selection1 = names(.)[1], selection2 = names(.)[2], common_cross = nrow(reduce(., intersect)),
                         stringsAsFactors = FALSE))

    # Other data
    tp_meta <- crossing(type = c("tp_base", "tp_select"), variable = c("var", "cov", "cor")) %>% 
      mutate(value = list(tp_cor, tp1_gencovar, tp1_genvar, tp_select_cor, tp_select_covar, tp_select_var))
    
    c1_meta <- data.frame(type = "c1_all", variable = c("cor", "cov", "var"), row.names = NULL, stringsAsFactors = FALSE) %>%
      tbl_df() %>% mutate(value = list(c1_cor, c1_covar, c1_var))

    ## Add the response and other results
    results_out[[i]] <- list(
      response = response_summary,
      common_cross = common_crosses,
      correlations = bind_rows(tp_meta, c1_meta)
    )

  }

  # Add the results to the core_df, remove core
  core_df %>%
    mutate(results = results_out) %>%
    select(-core)

}, mc.cores = n_cores)

# Bind and save
popvar_gencor_selection_simulation_out <- bind_rows(simulation_out)

# Save
save_file <- file.path(result_dir, "popvar_gencor_selection_simulation_results.RData")
save("popvar_gencor_selection_simulation_out", file = save_file)















