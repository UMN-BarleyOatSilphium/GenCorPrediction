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

s2_tp_BLUE <- s2_tidy_BLUE %>%
  filter(trait %in% traits, line_name %in% tp, year %in% 2014:2015, location %in% c("STP", "CRM")) %>%
  group_by(trait, environment, line_name) %>%
  summarize(value = mean(value)) %>%
  ungroup()
tp_tomodel <- map(trait_pairs, ~filter(s2_tp_BLUE, trait %in% ., line_name %in% tp_geno) %>% spread(trait, value))


# Model
tp_corG <- tp_tomodel %>% 
  map_df(~{
    df <- .
    

    ## Option 4 - sommer with GxE
    f <- as.formula(paste("cbind(", paste(tail(names(df), 2), collapse = ","), ") ~ line_name + environment"))
    mf <- model.frame(f, df)
    
    Y <- model.response(mf)
    X <- model.matrix(~ 1 + environment, mf) # Fixed effect of environment 
    Zg <- model.matrix(~ -1 + line_name, mf); colnames(Zg) <- levels(mf$line_name)
    Ze <- model.matrix(~ -1 + environment, mf); colnames(Ze) <- levels(mf$environment)
    # Zge <- model.matrix(~ -1 + line_name:environment, mf); colnames(Zge) <- levels(mf$line_name)
    
    Kg <- K
    Ke <- diag(ncol(Ze))
    Kge <- (Zg %*% Kg %*% t(Zg)) * tcrossprod(Ze)
    Zge <- diag(ncol(Kge))
    
    fit1 <- sommer::mmer(Y = Y, X = X, Z = list(corG = list(Z = Zg, K = Kg), corGE = list(Z = Zge, K = Kge)),
                         init = replicate(3, cov(Y), simplify = F))
    
    ## Calculate correlations
    cors <- map_df(fit1$var.comp, ~.[1,2] / prod(sqrt(diag(.))))
    
    ## Return df
    data.frame(trait1 = colnames(Y)[1], trait2 = colnames(Y)[2], cors, corP = cor(Y)[1,2], n = ncol(Zg), stringsAsFactors = FALSE)
    
    
  })

# trait1      trait2       corG       corE      corGE       units       corP
# 1 FHBSeverity HeadingDate -0.8419292 -0.2771820 -0.4395983  0.01424705 -0.2771820
# 2 FHBSeverity PlantHeight -0.4400359 -0.2594616 -1.2866123  0.05829586 -0.2594616
# 3 HeadingDate PlantHeight  0.4786676  0.7551931  0.7609902 -0.12331477  0.7551931

tp_corG <- as_data_frame(tp_corG)


## Signficance - use z transformation
tp_corG %>% 
  select(trait1, trait2, corG, n) %>%
  mutate(z = ztrans(corG), 
         se = 1 / sqrt(n - 3),
         upper = z + (1.96 * se), 
         lower = z - (1.96 * se)) %>% 
  mutate_at(vars(upper, lower), zexp) %>% 
  mutate(significant = map2_lgl(upper, lower, ~!between(0, .x, .y)), 
         significant = ifelse(significant, "*", ""))


# trait1      trait2        corG     n      z     se  upper  lower significant
# 1 FHBSeverity HeadingDate -0.842   175 -1.23  0.0762 -0.793 -0.880 *          
# 2 FHBSeverity PlantHeight -0.440   175 -0.472 0.0762 -0.312 -0.552 *          
# 3 HeadingDate PlantHeight  0.479   175  0.521 0.0762  0.585  0.356 * 









## Estimates of genetic correlation in each of the validation families

# Use the stage-two BLUEs to calculate correlation
vp_family_tomodel_experimental_BLUE <- vp_BLUE %>%
  filter(line_name %in% exper) %>%
  mutate(family = str_extract(line_name, "4[0-9]{3}"), std.error = NA)


## Calculate genetic correlation using REML and the BLUEs accounting for G and GE

# What families were measured for both traits in the pair
vp_family_tomodel2 <- trait_pairs %>%
  map(~{
    trs <- .
    
    if ("environment" %in% names(vp_family_tomodel_experimental_BLUE)) {
      distc_families <- vp_family_tomodel_experimental_BLUE %>% 
        distinct(trait, family, environment) %>%
        group_by(family, environment)
      
    } else {
      distc_families <- vp_family_tomodel_experimental_BLUE %>% 
        distinct(trait, family) %>%
        group_by(family)
      
    }
    
    families <- distc_families %>% 
      filter(trait %in% trs) %>% 
      filter(n() == length(trs))
    
    left_join(families, vp_family_tomodel_experimental_BLUE) %>% 
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
      mutate_at(vars(family, line_name), as.factor) %>%
      as.data.frame()
    
    # Fixed formula
    fixed_form <- as.formula(paste("cbind(", paste(sort(unique(df$trait)), collapse = ", "), ") ~ 1"))
    
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
        fit1 <- emmremlMultivariate(Y = t(Y), X = t(X), Z = t(Zg), K = t(diag(ncol(Zg))))
        
        ## Return variance components
        varcomp <- fit1$Vg
        data.frame(trait1 = colnames(Y)[1], trait2 = colnames(Y)[2], correlation = varcomp[1,2] / prod(sqrt(diag(varcomp))),
                   row.names = NULL, stringsAsFactors = FALSE)
        
      }) %>% ungroup()
    
    
  })



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
