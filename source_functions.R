## GenCorPred accessory functions
## 
## A script with useful functions used in the prediction analysis of predicting genetic
## correlations
## 



## Other/utility functions
# A function to assign cores to a data.frame
assign_cores <- function(df, n_core) {
  df$core <- sort(rep(seq(n_core), length.out = nrow(df)))
  return(df)
}


# A function to calculate heritability, BLUEs, and variance components from tidy phenotype data
summarize_pheno <- function(data, blue.model = c("lmer", "sommer")) {
  
  # Make sure necessary columns are in 'data'
  needed_cols <- c("trial", "environment", "location", "year", "line_name", "value", "std.error")
  stopifnot(all(needed_cols %in% names(data)))
  
  blue.model <- match.arg(blue.model)
  
  # If the number of trials/environment is > 1, fit a model to get the genotype mean
  # for a trait-environment combination
  mto_trial <- group_by(data, environment) %>% 
    summarize(n_trial = n_distinct(trial)) %>%
    filter(n_trial > 1)
  
  if (nrow(mto_trial) > 0) {
    
    env_mean <- data %>%
      filter(environment %in% mto_trial$environment) %>%
      group_by(environment) %>%
      do({
        data1 <- .
        fit1 <- lm(value ~ -1 + line_name + trial, data = data1)
        
        # Tidy
        tidy(fit1) %>% 
          select(term, estimate, std.error) %>% 
          filter(str_detect(term, "line_name")) %>% 
          mutate(term = str_replace(term, "line_name", "")) %>% 
          rename(line_name = term, value = estimate)
        
      })
    
    # Combine these results with the original data
    data1 <- bind_rows(
      data %>% filter(environment %in% mto_trial$environment) %>% distinct(environment, location, year, trait) %>% left_join(., env_mean, by = "environment"),
      data %>% filter(!environment %in% mto_trial$environment)
    ) %>%
      select(trial, names(.)) %>%
      arrange(environment)
    
  } else {
    data1 <- data
    
  }
  
  control <- lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.nRE = "ignore")
  wts <- data1$std.error^2
  
  # If the number of environments is < 2, drop relevant random effects
  if (n_distinct(data1$environment) < 2) {
    formula <- value ~ (1|line_name)
    exp <- "line_name / (line_name + (Residual / (n_e * n_r)))"
    
    ## Drop terms
    fit_noge <- lmer(formula = formula, data = data1, control = control, weights = wts)
    fit_nog <- lm(formula = value ~ 1, data = data1)
    
  } else {
    formula <- value ~ (1|line_name) + environment + (1|line_name:environment)
    exp <- "line_name / (line_name + (line_name:environment / n_e) + (Residual / (n_e * n_r)))"
    
    ## Drop terms
    fit_noge <- lmer(formula = value ~ (1|line_name) + environment, data = data1, control = control, weights = wts)
    fit_nog <- lmer(formula = value ~  environment + (1|line_name:environment), data = data1, control = control, weights = wts)
    
    
  }

  fit <- lmer(formula = formula, data = data1, control = control, weights = wts)
  
  plot_table <- xtabs(~line_name + environment, data1)
  
  # Get the harmonic mean of the number of environments / reps
  n_e <- plot_table %>% 
    ifelse(test = . > 1, 1, .) %>% 
    rowSums() %>% 
    harm_mean()
  
  n_r <- plot_table %>% 
    harm_mean()
  
  # Estimate heritability
  h2 <- herit(object = fit, exp = exp, n_r = n_r, n_e = n_e)
  
  

  # Calculate significance
  ge_sig <- lr_test(fit, fit_noge)
  g_sig <- lr_test(fit, fit_nog)
  sig_test <- bind_rows(g_sig, ge_sig) %>% 
    mutate(term = c("g", "ge")) %>% 
    select(term, names(.), -full_model)
  
  
  ## Split on whether to use lmer or sommer
  if (blue.model == "lmer") {
  
    ## Modify formula so line_name is fixed, then fit the model
    new_form <- tail(as.character(formula), 1) %>% 
      str_replace(string = ., pattern = "\\(1 \\| line_name\\)", "line_name") %>% str_c("value ~ -1 + ", .) %>% 
      as.formula()
    
    if (any(str_detect(new_form, "\\("))) {
      ## Now refit the model, but change genotype from random to fixed
      fit_blue <- lmer(formula = new_form, data = data1, control = control, weights = wts)
      
    } else {
      ## Now refit the model, but change genotype from random to fixed
      fit_blue <- lm(formula = new_form, data = data1)
      
    }
    
  
    
    # Tidy
    tidy_blue <- tidy(fit_blue) %>% 
      filter(str_detect(term, "line_name"), !str_detect(term, "sd")) %>%
      mutate(line_name = str_replace(term, "line_name", "")) %>% 
      select(line_name, value = estimate)
    
  } else if (blue.model == "sommer") {
    
    
    stopifnot(n_distinct(data$environment) > 1)
    
    ## Use sommer to calculate the genotype BLUEs
    mf <- model.frame(value ~ line_name + environment, data1)
    y <- model.response(mf)
    X <- model.matrix(~ -1 + line_name + environment, mf)
    Z <- model.matrix(~ -1 + line_name:environment, mf)
    K <- diag(ncol(Z))
    R <- diag(wts)
    
    fit_blue <- sommer::mmer(Y = y, X = X, Z = list(ge = list(Z = Z, K = diag(ncol(Z)))), R = list(unit = R))
    
    # Tidy
    tidy_blue <- fit_blue$beta.hat %>% 
      as.data.frame() %>% 
      rownames_to_column("term") %>% 
      rename(estimate = T1) %>%
      filter(str_detect(term, "line_name")) %>% 
      mutate(line_name = str_replace(term, "line_name", "")) %>% 
      select(line_name, value = estimate)
    
    
  }
  
  # Return all this nonsense
  data_frame(BLUE = list(tidy_blue), n_e = n_distinct(data$environment), h2 = list(h2), sig_test = list(sig_test))
  
}


# A function to calculate genetic variance
calc_varG <- function(data, method = c("lmer", "sommer")) {
  
  # Check the data input
  data <- droplevels(as.data.frame(data))
  method <- match.arg(method)
  
  # Check column names for the required columns
  needed_cols <- c("environment", "location", "year", "line_name", "value", "std.error", "family")
  stopifnot(all(needed_cols %in% names(data)))
  
  
  # Number of lines in the family
  n_lines <- n_distinct(data$line_name)
  n_env <- n_distinct(data$environment)
  
  plot_table <- xtabs(~line_name + environment, data)
  
  # Split based on the number of environments
  if (n_env > 1) {
  
    # Get the harmonic mean of the number of environments / reps
    n_e <- plot_table %>% 
      ifelse(test = . > 1, 1, .) %>% 
      rowSums() %>% 
      harm_mean()
    
    n_r <- plot_table %>% 
      harm_mean()
    
    wts <- data$std.error^2
    
    
    # Split flow based on method
    if (method == "lmer") {
      
      control <- lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.nRE = "ignore")
      formula <- value ~ (1|line_name) + environment + (1|line_name:environment)
      
      fit <- lmer(formula = formula, data = data, control = control, weights = wts, contrasts = list(environment = "contr.sum"))
      
      # Estimate heritability
      h2 <- herit(object = fit, exp = "line_name / (line_name + (line_name:environment / n_e) + (Residual / (n_e * n_r)))",
                  n_e = n_e, n_r = n_r)
      
      ## Drop terms
      fit_noge <- lmer(formula = value ~ (1|line_name) + environment, data = data, control = control, weights = wts)
      fit_nog <- lmer(formula = value ~  environment + (1|line_name:environment), data = data, control = control, weights = wts)
      
      # Calculate significance
      ge_sig <- lr_test(fit, fit_noge)
      fam_sig <- lr_test(fit, fit_nog)
      sig_test <- bind_rows(fam_sig, ge_sig) %>% 
        mutate(full_model = c("family", "ge")) %>% 
        rename(term_red = full_model)
      
      family_mean <- fixef(fit)[[1]]
      
    } else if (method == "sommer") {
      
      # Create the model matrices
      mf <- model.frame(value ~ line_name + environment, data)
      y <- model.response(mf)
      X <- model.matrix(~ 1 + environment, mf, contrasts.arg = list(environment = "contr.sum"))
      
      Zg <- model.matrix(~ -1 + line_name, mf)
      Kg <- diag(ncol(Zg))
      Zge <- model.matrix(~ -1 + line_name:environment, mf)
      Kge <- diag(ncol(Zge))
      
      R <- solve(diag(wts))
      
      # Fit the model
      fit <- sommer::mmer(Y = y, X = X, Z = list(g = list(Z = Zg, K = Kg), ge = list(Z = Zge, K = Kge)))
      
      varG <- fit$var.comp$g[1]
      varGE <- fit$var.comp$ge[1]
      varR <- fit$var.comp$units[1]
      
      h2 <- varG / (varG + (varGE / n_e) + (varR / (n_e + n_r)))
      var_comp <- data_frame(source = c("line_name:environment", "line_name", "Residual"),
                             variance = c(varGE, varG, varR))
      
      h2 <- list(heritability = h2, var_comp = var_comp)
      
      
      
      ## Drop terms
      fit_noge <- sommer::mmer(Y = y, X = X, Z = list(g = list(Z = Zg, K = Kg)))
      fit_nog <- sommer::mmer(Y = y, X = X, Z = list(ge = list(Z = Zge, K = Kge)))
      
      # Calculate significance
      ge_sig <- data_frame(term_red = "ge", statistic = -2 * (fit_noge$LL - fit$LL)) %>%
        mutate(df = 1, p_value = pchisq(q = statistic, df = df, lower.tail = FALSE))
      fam_sig <- data_frame(term_red = "family", statistic = -2 * (fit_nog$LL - fit$LL)) %>%
        mutate(df = 1, p_value = pchisq(q = statistic, df = df, lower.tail = FALSE))
      sig_test <- bind_rows(fam_sig, ge_sig)
      
      family_mean <- fit$beta.hat[1]
      
    }
    
  } else {
    
    n_r <- plot_table %>% 
      harm_mean()
    
    wts <- data$std.error^2
    
    
    # Split flow based on method
    if (method == "lmer") {
      
      control <- lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.nRE = "ignore")
      formula <- value ~ (1|line_name)
      
      fit <- lmer(formula = formula, data = data, control = control, weights = wts)
      
      # Estimate heritability
      h2 <- herit(object = fit, exp = "line_name / (line_name + (Residual / (n_r)))",
                  n_r = n_r)
      
      ## Drop terms
      fit_noge <- fit
      fit_nog <- lm(formula = value ~ 1, data = data)
      
      # Calculate significance
      ge_sig <- data_frame(term_red = "ge", statistic = -2 * (as.numeric(logLik(fit_noge)) - as.numeric(logLik(fit)))) %>%
        mutate(df = 1, p_value = pchisq(q = statistic, df = df, lower.tail = FALSE))
      fam_sig <- data_frame(term_red = "family", statistic = -2 * (as.numeric(logLik(fit_nog)) - as.numeric(logLik(fit)))) %>%
        mutate(df = 1, p_value = pchisq(q = statistic, df = df, lower.tail = FALSE))
      sig_test <- bind_rows(fam_sig, ge_sig)
      
      family_mean <- fixef(fit)[[1]]
      
    } else if (method == "sommer") {
      
      # Create the model matrices
      mf <- model.frame(value ~ line_name, data)
      y <- model.response(mf)
      X <- model.matrix(~ 1, mf, contrasts.arg = list(environment = "contr.sum"))
      
      Zg <- model.matrix(~ -1 + line_name, mf)
      Kg <- diag(ncol(Zg))
      
      R <- solve(diag(wts))
      
      # Fit the model
      fit <- sommer::mmer(Y = y, X = X, Z = list(g = list(Z = Zg, K = Kg)))
      
      varG <- fit$var.comp$g[1]
      varR <- fit$var.comp$units[1]
      
      h2 <- varG / (varG + (varR / (n_r)))
      var_comp <- data_frame(source = c("line_name", "Residual"),
                             variance = c(varG, varR))
      
      h2 <- list(heritability = h2, var_comp = var_comp)
      
      
      
      ## Drop terms
      fit_nog <- sommer::mmer(Y = y, X = X, Z = list(g = list(Z = diag(length(y)), K = diag(length(y)))))

      # Calculate significance
      ge_sig <- data_frame(term_red = "ge", statistic = -2 * (fit$LL - fit$LL)) %>%
        mutate(df = 1, p_value = pchisq(q = statistic, df = df, lower.tail = FALSE))
      fam_sig <- data_frame(term_red = "family", statistic = -2 * (fit_nog$LL - fit$LL)) %>%
        mutate(df = 1, p_value = pchisq(q = statistic, df = df, lower.tail = FALSE))
      sig_test <- bind_rows(fam_sig, ge_sig)
      
      family_mean <- fit$beta.hat[1]
      
    }
    
    
  }
  
  # Return all this nonsense
  data_frame(family_mean = family_mean, fit = list(fit), h2 = list(h2), sig_test = list(sig_test))
  
}



# A function to calculate family mean and superior progeny mean
calc_mean <- function(data, i = 0.1) {
  
  # Check the data input
  data <- droplevels(as.data.frame(data))
  stopifnot(between(i, 0, 1))
  
  # Check column names for the required columns
  needed_cols <- c("environment", "location", "year", "line_name", "value", "std.error", "family")
  stopifnot(all(needed_cols %in% names(data)))
  
  # Convert some variables to factors
  data1 <- mutate_at(data, vars(line_name, environment), as.factor)
  
  
  # Number of lines in the family
  n_lines <- n_distinct(data1$line_name)
  n_env <- n_distinct(data1$environment)
  
  plot_table <- xtabs(~line_name + environment, data1)
  
  # Split based on the number of environments
  if (n_env > 1) {
    
    wts <- data1$std.error^2

    control <- lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.nRE = "ignore")
    formula <- value ~ 1 + line_name + environment + (1|line_name:environment)
    
    fit <- lmer(formula = formula, data = data1, control = control, weights = wts, 
                contrasts = list(environment = "contr.sum", line_name = "contr.sum"))
    
    # Get the fixed effect estimates
    coefs <- tidy(fit) %>% 
      filter(str_detect(term, "line_name|Intercept"), group == "fixed") %>%
      select(term, estimate) %>%
      add_row(term = "last_line", estimate = -sum(tail(.$estimate, -1))) %>%
      mutate(term = c("family_mean", levels(data1$line_name)),
             mean = c(estimate[1], estimate[-1] + estimate[1]))

    
  } else {
    
    formula <- value ~ 1 + line_name
    
    fit <- lm(formula = formula, data = data1, contrasts = list(line_name = "contr.sum"))
    
    # Get the fixed effect estimates
    coefs <- tidy(fit) %>% 
      filter(str_detect(term, "line_name|Intercept")) %>%
      select(term, estimate) %>%
      add_row(term = "last_line", estimate = -sum(tail(.$estimate, -1))) %>%
      mutate(term = c("family_mean", levels(data1$line_name)),
             mean = c(estimate[1], estimate[-1] + estimate[1]))
    
  }
  
  geno_means <- subset(coefs, term != "family_mean", mean, drop = T)
  
  # Calculate the superior progeny mean
  mu_sp <- mean(geno_means[geno_means <= quantile(geno_means, i)])
    
  # Return all this nonsense
  data_frame(means = list(coefs), family_mean = coefs$mean[1], mu_sp = mu_sp)
  
}




# A function to return a tidy output from PopVar
tidy.popvar <- function(x) {
  x$predictions %>% 
    map(as_data_frame) %>% 
    map(~mutate_all(., unlist)) %>%
    list(., names(.)) %>%
    pmap(.f = function(df, tr) {
      mutate(df, trait = tr) }) %>%
    bind_rows() %>%
    select(trait, names(.)) %>%
    mutate(trait = str_replace(trait, "_param.df", ""))
}
  
  

## Function to correct a genetic model to force the desired genetic correlation
adj_multi_gen_model <- function(genome, geno, gencor) {
  
  # Grab the genetic model
  gen_model <- genome$gen_model
  
  # Get the positions of QTL with non-zero effect
  eff_qtl <- lapply(gen_model, function(x) which(x$add_eff != 0))
  # Remove the QTL with non-zero effect
  gen_model1 <- lapply(gen_model, subset, add_eff != 0, drop = FALSE)
  
  # Get the QTL names
  qtl_names <- lapply(gen_model1, "[[", "qtl_name")
  
  ## Pull the genotypes for these QTL
  qtl_geno <- pull_genotype(genome = genome, geno = geno, loci = qtlnames(genome)) - 1
  # Calculate D between the QTL, then subset for trait1 vs trait2
  D <- cor(qtl_geno)[qtl_names[[1]], qtl_names[[2]]]
  
  
  # Grab the qtl effects for the first trait (these will be unadulterated)
  trait1_qtl_eff <- gen_model1[[1]][,"add_eff"]
  
  # Create a list with these effects, then randomize these effects for the second trait
  qtl_eff <- list(
    setNames(trait1_qtl_eff, qtl_names[[1]]),
    setNames(sample(trait1_qtl_eff), qtl_names[[2]])
  )
  
  # qtl_eff <- list(
  #   setNames(gen_model1[[1]][,"add_eff"], qtl_names[[1]]),
  #   setNames(gen_model1[[2]][,"add_eff"], qtl_names[[2]])
  # )
  
  ## Adjust the trait 2 effects by the desired correlation
  A <- reduce(map(qtl_eff, as.matrix), cbind)
  Sigma <- rbind(c(1, gencor), c(gencor, 1))
  
  qtl_eff_adj <- A %*% chol(Sigma)
  # qtl_eff_adj <- do.call(cbind, qtl_eff)
  

  # Revise the effects for trait2
  qtl_eff_adj[,2] <- qtl_eff_adj[,2] %*% D
  
  # Add these effects to gen_model1
  gen_model1[[1]]$add_eff <- qtl_eff_adj[,1]
  gen_model1[[2]]$add_eff <- qtl_eff_adj[,2]
  
  # Add this gen_model back into the genome
  gen_model[[1]][eff_qtl[[1]],] <- gen_model1[[1]]
  gen_model[[2]][eff_qtl[[2]],] <- gen_model1[[2]]
  
  # Add these effects back into the genome
  genome$gen_model[[1]] <- gen_model[[1]]
  genome$gen_model[[2]] <- gen_model[[2]]
  
  
  # Return the genome
  return(genome)
  
}


### pop.predict using bi-variate models for marker effects
pop.predict.bivariate <- function (G.in = NULL, y.in = NULL, map.in = NULL, crossing.table = NULL, 
                                   parents = NULL, tail.p = 0.1, nInd = 200, map.plot = F, min.maf = 0.01, 
                                   mkr.cutoff = 0.5, entry.cutoff = 0.5, remove.dups = T, impute = "EM", 
                                   nSim = 25, frac.train = 0.6, nCV.iter = 100, nFold = NULL, 
                                   nFold.reps = 1, nIter = 12000, burnIn = 3000, models = c("rrBLUP", 
                                                                                            "BayesA", "BayesB", "BayesC", "BL", "BRR"), return.raw = F) 
{
  if (is.null(G.in)) 
    stop("Must provide a genotype (G.in) file.")
  if (is.null(y.in)) 
    stop("Must provide a phenotype (y.in) file.")
  if (is.null(map.in)) 
    stop("Must provide a map (map.in) file.")
  if (!is.null(min.maf) & min.maf >= 1) 
    stop("min.maf must be within the range [0, 1)")
  if (!is.null(entry.cutoff) & entry.cutoff > 1) 
    stop("entry.cutoff must be within the range (0, 1]")
  if (!is.null(mkr.cutoff) & mkr.cutoff > 1) 
    stop("mkr.cutoff must be within the range (0, 1]")
  if (impute == "pass") {
    min.maf = 0
    mkr.cutoff = 1
    entry.cutoff = 1
  }
  G.entries <- as.character(G.in[-1, 1])
  entries.removed <- NULL
  entries.to.remove <- c()
  G.markers <- as.character(t(G.in[1, -1]))
  map.markers <- as.character(map.in[, 1])
  mkrs.removed <- NULL
  mkrs.to.remove <- c()
  if (!all(G.markers %in% map.markers) & all(map.markers %in% 
                                             G.markers)) 
    stop("Markers in Genotype matrix and genetic map do not completely match.")
  map <- map.in[order(map.in[, 2], map.in[, 3]), ]
  G.mat <- as.matrix(G.in[-1, -1])
  class(G.mat) <- "numeric"
  G.mat <- G.mat[, order(match(map.markers, G.markers))]
  if (impute != "pass" && !all(unique(G.mat[, 1]) %in% c(-1, 
                                                         0, 1, NA))) 
    stop("\nNon-imputed genotypes need to be coded as -1, 0, 1.\nIf imputed genotypes are passed set impute = 'pass'.")
  mkrs.to.remove <- c()
  if (min.maf > 0) {
    maf.list <- apply(G.mat, 2, maf.filt)
    mkrs.to.remove <- c(mkrs.to.remove, which(maf.list < 
                                                min.maf))
  }
  if (mkr.cutoff < 1) {
    mkrNA.list <- apply(G.mat, 2, function(M) {
      return(length(which(is.na(M)))/length(M))
    })
    mkrs.to.remove <- unique(c(mkrs.to.remove, which(mkrNA.list > 
                                                       mkr.cutoff)))
  }
  if (length(mkrs.to.remove > 0)) {
    G.mat <- G.mat[, -mkrs.to.remove]
    map <- map[-mkrs.to.remove, ]
    mkrs.removed <- map.markers[mkrs.to.remove]
    map.markers <- map.markers[-mkrs.to.remove]
  }
  entries.to.remove <- c()
  if (remove.dups) 
    entries.to.remove <- c(entries.to.remove, which(duplicated.array(G.mat)))
  if (entry.cutoff < 1) {
    entryNA.list <- apply(G.mat, 1, function(E) {
      return(length(which(is.na(E)))/length(E))
    })
    entries.to.remove <- unique(c(entries.to.remove, which(entryNA.list > 
                                                             entry.cutoff)))
  }
  if (length(entries.to.remove > 0)) {
    G.mat <- G.mat[-entries.to.remove, ]
    entries.removed <- G.entries[entries.to.remove]
    G.entries <- G.entries[-entries.to.remove]
    if (!is.null(crossing.table)) {
      cross.tabl.1col <- rbind(cbind(as.numeric(row.names(crossing.table)), 
                                     as.character(crossing.table[, 1])), cbind(as.numeric(row.names(crossing.table)), 
                                                                               as.character(crossing.table[, 2])))
      tab.rows.2remove <- as.numeric(unique(unlist(sapply(entries.removed, 
                                                          function(X) {
                                                            return(cross.tabl.1col[grep(X, cross.tabl.1col[, 
                                                                                                           2]), 1])
                                                          }))))
      if (length(tab.rows.2remove) > 0) 
        crossing.table <- crossing.table[-tab.rows.2remove, 
                                         ]
    }
    if (!is.null(parents)) {
      parents <- parents[!parents %in% entries.removed]
    }
  }
  y <- as.matrix(y.in[match(G.entries, as.character(y.in[, 
                                                         1])), ])
  y.entries <- as.character(y[, 1])
  traits <- as.character(colnames(y))[-1]
  nTraits <- length(traits)
  name4out <- sample(10000:99999, 1)
  t.map <- t(map)
  rownames(t.map) <- NULL
  map4out <- cbind(c("pheno", "", ""), t.map)
  write.table(map4out, paste("map.tmp_", name4out, ".csv", 
                             sep = ""), row.names = F, col.names = F, sep = ",")
  options(warn = -1)
  read.map.out <- capture.output(read.map <- qtl::read.cross(format = "csv", 
                                                             crosstype = "riself", file = paste("map.tmp_", name4out, 
                                                                                                ".csv", sep = ""), na.strings = "NA"))
  print(paste("Number of Markers Read in: ", unlist(strsplit(read.map.out[3], 
                                                             split = " "), recursive = T)[2], sep = ""), quote = F)
  unlink(paste("map.tmp_", name4out, ".csv", sep = ""))
  map_t1 <- qtl::pull.map(read.map)
  options(warn = 0)
  if (map.plot == T) 
    qtl::plot.map(map_t1)
  if (impute == "EM") 
    G.imp <- rrBLUP::A.mat(G.mat, min.MAF = 0, max.missing = 1, 
                           impute.method = "EM", return.imputed = T)$imputed
  if (impute == "mean") 
    G.imp <- rrBLUP::A.mat(G.mat, min.MAF = 0, max.missing = 1, 
                           impute.method = "mean", return.imputed = T)$imputed
  if (impute == "pass") 
    G.imp <- G.mat
  
  ## List of trait pairs
  trait_pairs <- combn(x = traits, m = 2, simplify = FALSE)
  nTraitPairs <- length(trait_pairs)
  mkr_effs_list <- beta.list <- list()
  
  param_df_list <- setNames(object = vector("list", nTraitPairs), nm =sapply(X = trait_pairs, paste0, collapse = ":"))
  
  
  for (t in 1:nTraitPairs) {
    
    trait_pair <- trait_pairs[[t]]
    y_notNAs <- !is.na(y[, trait_pair])
    y_TP <- apply(X = y[y_notNAs[,1], trait_pair], MARGIN = 2, FUN = as.numeric)
    TP.entries <- y.entries[y_notNAs[,1]]
    G_TP <- G.imp[y_notNAs[,1], ]
    
    ## Predict marker effects
    X <- matrix(data = 1, nrow = nrow(y_TP), ncol = 1)
    K <- diag(ncol(G_TP))
    
    # fit <- EMMREML::emmremlMultivariate(Y = t(y_TP), X = t(X), Z = t(G_TP), K = K)
    # 
    # # Extract model results
    # beta.list[[t]] <- beta <- fit$Bhat[,1]
    # mkr_effs_list[[t]] <- mkr_effs.mat <- t(fit$Gpred)
    
    fit <- sommer::mmer(Y = y_TP, X = X, Z = list(g = list(Z = G_TP, K = K)), silent = TRUE)
    
    # Extract model results
    beta.list[[t]] <- beta <- fit$beta.hat[1,]
    mkr_effs_list[[t]] <- mkr_effs.mat <- fit$u.hat$g

    if (!is.null(crossing.table)) {
      crossing.table <- as.matrix(crossing.table)
      crossing.mat <- PopVar:::par.position(crossing.table, par.entries = G.entries)$parent.position
      crosses.possible <- PopVar:::par.position(crossing.table, par.entries = G.entries)$crosses.possible
    }
    
    par.BVs <- G.imp %*% mkr_effs.mat
    # Rescale the breeding values
    par_BV_rescale <- t(t(par.BVs) + beta)
    
    ## Temporary data.frame
    traits <- names(beta)
    nTraits <- length(traits)
    df.tmp <- data.frame(cbind(crosses.possible[, 1], crosses.possible[,2], matrix(list(rep(NA, times = nSim)), nrow = nrow(crosses.possible), ncol = (8 + 3 * (nTraits - 1)))))
    
    names(df.tmp)[1:10] <- c("Par1", "Par2", "midPar.Pheno", 
                             "midPar.GEBV", "pred.mu", "pred.mu_sd", "pred.varG", 
                             "pred.varG_sd", "mu.sp_low", "mu.sp_high")
    for (n in 1:nTraits) {
      if (n == 1) 
        param.dfs <- list()
      param.dfs[[n]] <- as.matrix(df.tmp)
    }
    names(param.dfs) <- paste(traits, "_param.df", sep = "")
    cat("\n")
    cat(paste("\nBrewing", nSim, "populations of", nInd, "individuals for each cross... Please be patient", 
              sep = " "))
    cat("\n")
    prog.bar <- txtProgressBar(min = 1, max = (nrow(crossing.mat) * 
                                                 nSim), style = 3)
    p = 1
    M <- nInd
    for (s in 1:nSim) {
      sim.pop <- qtl::sim.cross(map_t1, type = "riself", n.ind = M, model = NULL)
      qtl::write.cross(sim.pop, "csv", paste("sim.pop.tmp_", name4out, sep = ""))
      pop.mat <- as.matrix(read.csv(paste("sim.pop.tmp_", name4out, ".csv", sep = ""), header = T))[3:(M + 2), 2:(ncol(G_TP) + 1)]
      unlink(paste("sim.pop.tmp_", name4out, ".csv", sep = ""))
      
      for (z in 1:nrow(crossing.mat)) {
        setTxtProgressBar(prog.bar, p)
        pop.mat2 <- matrix(NA, nrow = nrow(pop.mat), ncol = ncol(pop.mat))
        par1 <- G.imp[crossing.mat[z, 1], ]
        par2 <- G.imp[crossing.mat[z, 2], ]
        for (r in 1:M) {
          pop.mat2[r, which(pop.mat[r, ] == "A")] <- par1[which(pop.mat[r, ] == "A")]
          pop.mat2[r, which(pop.mat[r, ] == "B")] <- par2[which(pop.mat[r,] == "B")]
        }
        mkr.has.0 <- apply(pop.mat2, 2, function(X) {
          return(length(which(X == 0)))
        })
        replace.0.mat <- rbind(which(mkr.has.0 != 0), mkr.has.0[which(mkr.has.0 != 0)])
        if (ncol(replace.0.mat) > 0) {
          for (b in 1:ncol(replace.0.mat)) {
            pop.mat2[which(pop.mat2[, replace.0.mat[1, 
                                                    b]] == 0), replace.0.mat[1, b]] <- sample(c(1, 
                                                                                                -1), size = replace.0.mat[2, b], replace = T, 
                                                                                              prob = c(0.5, 0.5))
          }
        }
        prog_pred.mat <- pop.mat2 %*% mkr_effs.mat
        prog_pred_mat_rescale <- t(t(prog_pred.mat) + beta)
        
        for (n in 1:nTraits) {
          if (s == 1) {
            if (nTraits > 1) 
              colnames(param.dfs[[n]])[11:(10 + 3 * (nTraits - 
                                                       1))] <- c(paste("low.resp_", traits[-n], 
                                                                       sep = ""), paste("high.resp_", traits[-n], 
                                                                                        sep = ""), paste("cor_w/_", traits[-n], 
                                                                                                         sep = ""))
            param.dfs[[n]][[z, "midPar.Pheno"]][s] <- 0.5 * 
              (as.numeric(y[crossing.mat[z, 1], n + 1]) + 
                 as.numeric(y[crossing.mat[z, 2], n + 1]))
            param.dfs[[n]][[z, "midPar.GEBV"]][s] <- 0.5 * 
              (as.numeric(par.BVs[crossing.mat[z, 1], n]) + 
                 as.numeric(par.BVs[crossing.mat[z, 2], 
                                    n]))
          }
          param.dfs[[n]][[z, "pred.mu"]][s] <- mean(prog_pred.mat[, 
                                                                  n])
          param.dfs[[n]][[z, "pred.mu_sd"]][s] <- mean(prog_pred.mat[, 
                                                                     n])
          param.dfs[[n]][[z, "pred.varG"]][s] <- var(prog_pred.mat[, 
                                                                   n])
          param.dfs[[n]][[z, "pred.varG_sd"]][s] <- var(prog_pred.mat[, 
                                                                      n])
          param.dfs[[n]][[z, "mu.sp_low"]][s] <- PopVar:::tails(prog_pred.mat[, 
                                                                     n], tail.p = tail.p)[2]
          param.dfs[[n]][[z, "mu.sp_high"]][s] <- PopVar:::tails(prog_pred.mat[, 
                                                                      n], tail.p = tail.p)[1]
          if (nTraits > 1) {
            index <- 1
            for (n2 in (1:nTraits)[-n]) {
              param.dfs[[n]][[z, 10 + index]][s] <- mean(prog_pred.mat[, 
                                                                       n2][which(prog_pred.mat[, n] <= quantile(prog_pred.mat[, 
                                                                                                                              n], probs = tail.p))], na.rm = T)
              param.dfs[[n]][[z, 10 + (nTraits - 1) + index]][s] <- mean(prog_pred.mat[, 
                                                                                       n2][which(prog_pred.mat[, n] >= quantile(prog_pred.mat[, 
                                                                                                                                              n], probs = 1 - tail.p))], na.rm = T)
              param.dfs[[n]][[z, 10 + 2 * (nTraits - 1) + 
                                index]][s] <- cor(prog_pred.mat[, n], prog_pred.mat[, 
                                                                                    n2], use = "complete.obs")
              index <- index + 1
            }
          }
        }
        p <- p + 1
      }
    }
    preds.per.sim <- param.dfs
    for (n in 1:nTraits) {
      col.names <- colnames(param.dfs[[n]])
      for (c in 3:length(col.names)) {
        name.tmp <- col.names[c]
        if (name.tmp %in% c("pred.mu_sd", "pred.varG_sd")) 
          param.dfs[[n]][, c] <- sapply(param.dfs[[n]][, 
                                                       c], FUN = sd, na.rm = T)
        if (!name.tmp %in% c("pred.mu_sd", "pred.varG_sd")) 
          param.dfs[[n]][, c] <- sapply(param.dfs[[n]][, 
                                                       c], FUN = mean, na.rm = T)
      }
    }
    
    ## Add param_df to list
    param_df_list[[t]] <- param.dfs
    
  } # Close loop over trait pairs
    
  return(list(predictions = param_df_list))
}
  
  
  
  
  
  
  
  
  