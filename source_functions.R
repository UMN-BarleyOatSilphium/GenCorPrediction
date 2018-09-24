## S2MET Functions
## 
## A script with useful functions used in the prediction analysis of the S2MET
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
  
  ## Adjust the trait 2 effects by the desired correlation
  qtl_eff_adj <- reduce(map(qtl_eff, as.matrix), cbind) %*% chol(rbind(c(1, gencor), c(gencor, 1)))

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
  
    
 