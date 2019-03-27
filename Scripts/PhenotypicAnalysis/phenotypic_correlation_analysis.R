## GenCorPrediction
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
library(pbr) # Available from https://github.com/neyhartj/pbr
library(broom)
library(ggridges)
library(cowplot)
library(modelr)
library(EMMREML)



boot_rep <- 1000
alpha <- 0.05
i <- 0.1


# Load the PVV BLUEs
load(file.path(data_dir, "PVV_BLUE.RData")) # This is available from the GitHub repository
# # Load the S2 BLUEs and filter
# load(file.path(gdrive_dir, "BarleyLab/Breeding/PhenotypicData/Final/MasterPhenotypes/S2_tidy_BLUE.RData"))


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
    # fit <- sommer::mmer(Y = Y, X = X, Z = list(g = list(Z = Z, K = K)))
    
    
    ## Return variance components
    varcompG <- fit$Vg
    varcompE <- fit$Ve

    # varcomp <- fit$var.comp$g
    data_frame(trait1 = colnames(Y)[1], trait2 = colnames(Y)[2], corG = varcompG[1,2] / prod(sqrt(diag(varcompG))), corE = varcompE[1,2] / prod(sqrt(diag(varcompE))))
    
  })


# trait1      trait2 correlation
# 1 FHBSeverity HeadingDate  -0.9896474
# 2 FHBSeverity PlantHeight  -0.6121457
# 3 HeadingDate PlantHeight   0.3841620




# Use the stage-two BLUEs to calculate correlation
vp_family_tomodel1 <- vp_BLUE %>%
  filter(line_name %in% exper) %>%
  mutate(family = str_extract(line_name, "4[0-9]{3}"), std.error = NA)


# What families were measured for both traits in the pair
vp_family_tomodel2 <- trait_pairs %>%
  map(~{
    trs <- .
    
    if ("environment" %in% names(vp_family_tomodel1)) {
      distc_families <- vp_family_tomodel1 %>% 
        distinct(trait, family, environment) %>%
        group_by(family, environment)
      
    } else {
      distc_families <- vp_family_tomodel1 %>% 
        distinct(trait, family) %>%
        group_by(family)
      
    }
      
    families <- distc_families %>% 
      filter(trait %in% trs) %>% 
      filter(n() == length(trs))
    
    left_join(families, vp_family_tomodel1) %>% 
      ungroup()
  })


## Fit a model per family
## Use EMMREML


vp_family_corG1 <- vp_family_tomodel2 %>% 
  map_df(~{
    df <- .
    
    ## Modify the data.frame
    df1 <- df %>%
      dplyr::select(family, line_name, trait, value) %>%
      spread(trait, value) %>%
      # rename_at(vars(which(names(.) %in% traits)), funs(c("trait1", "trait2"))) %>%
      mutate_at(vars(family, line_name), as.factor) %>%
      as.data.frame()


    df1 %>%
      group_by(family) %>%
      do({
        df2 <- .

        f <- as.formula(paste("cbind(", paste(tail(names(df2), 2), collapse = ", "), ") ~ line_name"))
        mf <- model.frame(f, df2)
        # Bind the traits
        Y <- model.response(mf)
        X <- model.matrix(~ 1, mf)

        # Random effects
        Zg <- model.matrix(~ 1 + line_name, droplevels(mf))

        # Fit the model
        fit2 <- emmremlMultivariate(Y = t(Y), X = t(X), Z = t(Zg), K = t(diag(ncol(Zg))))
        
        # ## Simulate for bootstrapping
        # y_sim <- replicate(n = boot_rep, X %*% t(fit2$Bhat) + mvtnorm::rmvnorm(n = nrow(Y), mean = c(0, 0), sigma = fit2$Vg) + mvtnorm::rmvnorm(n = nrow(Y), mean = c(0, 0), sigma = fit2$Ve), simplify = FALSE)
        # y_sim_fit <- map(y_sim, ~emmremlMultivariate(Y = t(.), X = t(X), Z = t(Zg), K = t(diag(ncol(Zg)))))
        # 
        
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




## For each trait, find the best i% of lines
vp_family_best_lines <- vp_family_tomodel1 %>% 
  split(.$trait) %>% 
  map_df(~split(., .$family) %>% 
           map_df(~top_n(x = ., n = round(i * nrow(.)), wt = -value) %>% 
                    group_by(trait, family) %>% nest(line_name)))

## Now calculate the mean of the progeny for that trait
vp_family_musp <- vp_family_best_lines %>%
  mutate(musp = map2_dbl(.x = trait, .y = data, ~mean(subset(vp_family_tomodel1, trait == .x & line_name %in% .y$line_name, value, drop = T))))

## Now calculate the correlated superior progeny mean for the second trait
vp_family_muspC <- vp_family_musp %>%
  mutate(muspC = map2(.x = trait, .y = data, ~filter(vp_family_tomodel1, trait != .x, line_name %in% .y$line_name) %>% group_by(trait) %>% 
                        summarize(muspC = mean(value)) %>% rename(traitC = trait))) %>%
  unnest(muspC) %>%
  rename(trait1 = trait, trait2 = traitC)


## Save the correlation calculations
save("vp_family_corG1", "vp_family_muspC", "tp_corG", file = file.path(result_dir, "correlation_analysis.RData"))











