## PopVarValidation
## Calculation and analysis of genetic correlation
##
## 
## Author: Jeff Neyhart
## Last modified: August 29, 2018
## 

# Load the source script
repo_dir <- getwd()
source(file.path(repo_dir, "source.R"))

# Load the pbr library
library(pbr)
library(broom)
library(ggridges)
library(cowplot)
library(modelr)
library(EMMREML)


# Load the PVV BLUEs
load(file.path(data_dir, "PVV_BLUE.RData"))
# Load the S2 BLUEs and filter
load(file.path(gdrive_dir, "BarleyLab/Breeding/PhenotypicData/Final/MasterPhenotypes/S2_tidy_BLUE.RData"))

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
    
    ## Return variance components
    varcomp <- fit$Vg
    data.frame(trait1 = colnames(Y)[1], trait2 = colnames(Y)[2], correlation = varcomp[1,2] / prod(sqrt(diag(varcomp))),
               row.names = NULL, stringsAsFactors = FALSE)
    
  })


vp_family_tomodel <- s2_tidy_BLUE %>%
  filter(trait %in% traits, line_name %in% c(pot_pars, exper), str_detect(trial, "PVV"))


boot_rep <- 1000
alpha <- 0.05
i <- 0.1


## First if the same line is measured in two different trials, calculate an environment mean
vp_family_tomodel1 <- vp_family_tomodel %>%
  filter(line_name %in% exper) %>%
  group_by(trait, environment) %>%
  do({
    df <- .
    
    if (n_distinct(df$trial) > 1) {
      
      fit <- lm(value ~ -1 + line_name + trial, data = df)
      fit_tidy <- tidy(fit) %>% filter(str_detect(term, "line_name")) %>% 
        mutate(line_name = str_replace(term, "line_name", "")) %>% 
        select(line_name, value = estimate, std.error)
      
      df %>% 
        mutate(trial = NA) %>% 
        distinct(trial, environment, location, year, trait) %>% 
        cbind(., fit_tidy)
      
    } else {
      
      df
    }
    
  }) %>% ungroup()


## Calculate genetic correlation via REML by fitting a model with G and GE

# ## Use the stage-one BLUEs to calculate correlation
# vp_family_tomodel1 <- vp_family_tomodel1 %>% 
#   mutate(family = str_extract(line_name, "4[0-9]{3}"))

# Use the stage-two BLUEs to calculate correlation
vp_family_tomodel1 <- vp_BLUE %>%
  filter(line_name %in% exper) %>%
  mutate(family = str_extract(line_name, "4[0-9]{3}"))


# What families were measured for both traits in the pair
vp_family_tomodel2 <- trait_pairs %>%
  map(~{
    trs <- .
    families <- vp_family_tomodel1 %>% 
      distinct(trait, family) %>% 
      group_by(family) %>%
      filter(trait %in% trs) %>% 
      filter(n() == length(trs))
    
    left_join(families, vp_family_tomodel1, by = c("trait", "family")) %>% 
      ungroup()
  })

## Fit a model per family
## Use EMMREML


vp_family_corG1 <- vp_family_tomodel2 %>% 
  map_df(~{
    df <- .
    
    # df %>%
    #   select(family, environment, line_name, trait, value) %>%
    #   spread(trait, value) %>%
    #   group_by(family) %>%
    #   do({
    #     df2 <- .
    #     
    #     fit2 <- sommer::mmer2(fixed = Y ~ environment, 
    #                           random = ~us(trait):line_name + us(trait):line_name:environment, 
    #                           data = df2, rcov = ~ us(traits):units, silent = T)
    #     
    #     ## Return variance components
    #     varcomp <- fit2$var.comp$line_name
    #     data.frame(trait1 = colnames(Y)[1], trait2 = colnames(Y)[2], correlation = varcomp[1,2] / prod(sqrt(diag(varcomp))),
    #                row.names = NULL, stringsAsFactors = FALSE)
    #     
    #   })
    
    df %>%
      select(family, line_name, trait, value) %>%
      spread(trait, value) %>%
      group_by(family) %>%
      do({
        df2 <- .

        f <- as.formula(paste("cbind(", paste(tail(names(df2), 2), collapse = ", "), ") ~ line_name"))
        mf <- model.frame(f, df2)
        # Bind the traits
        Y <- model.response(mf)
        X <- model.matrix(~ 1, mf)

        # Random effects
        Zg <- model.matrix(~ 1 + line_name, mf)

        # Fit the model
        fit2 <- emmremlMultivariate(Y = t(Y), X = t(X), Z = t(Zg), K = t(diag(ncol(Zg))))
        
        ## Return variance components
        varcomp <- fit2$Vg
        data.frame(trait1 = colnames(Y)[1], trait2 = colnames(Y)[2], correlation = varcomp[1,2] / prod(sqrt(diag(varcomp))),
                   row.names = NULL, stringsAsFactors = FALSE)
        
      }) %>% ungroup()

    
  })



# Distributions
vp_family_corG1 %>%
  ggplot(aes(x = correlation)) +
  geom_histogram() +
  facet_grid(trait1 ~ trait2)

# Min, max, mean for each trait pair
vp_family_corG1 %>% 
  group_by(trait1, trait2) %>% 
  summarize_at(vars(correlation), funs(min, max, mean))

# trait1      trait2         min   max     mean
# 1 FHBSeverity HeadingDate -0.724 0.638 -0.0662 
# 2 FHBSeverity PlantHeight -0.672 0.635 -0.00100
# 3 HeadingDate PlantHeight -0.663 0.686 -0.106 

## Save the correlation calculations
save("vp_family_corG1", "tp_corG", file = file.path(result_dir, "correlation_analysis.RData"))











