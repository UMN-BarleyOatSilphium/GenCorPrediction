## GenCorPrediction - prediction accuracy and analysis 
##
## This script will look at the prediction accuracy of family correlation
## 
## Author: Jeff Neyhart
## Last modified: 22 July 2019
## 

# Load the source script
repo_dir <- getwd()
source(file.path(repo_dir, "source.R"))

library(cowplot)
library(gridExtra)

# Boot reps
boot_reps <- 1000
alpha <- 0.05

## Load the predictions and the results
# Predictions
load(file.path(result_dir, "prediction_results.RData"))
# Results
load(file.path(result_dir, "correlation_analysis.RData"))


## First subset the relevant columns
popvar_pred <- list(pred_results_realistic, pred_results_relevant) %>% 
  setNames(c("realistic", "relevant")) %>%
  map(~{
    filter(., trait %in% traits) %>%
      left_join(., cross_list, by = c("Par1" = "parent1", "Par2" = "parent2")) %>%
      select(parent1 = Par1, parent2 = Par2, family, trait, note, pred_mu = pred.mu, pred_varG = pred.varG, musp_high = mu.sp_high,
             musp_low = mu.sp_low, cor_HeadingDate = `cor_w/_HeadingDate`, cor_PlantHeight = `cor_w/_PlantHeight`,
             cor_FHBSeverity = `cor_w/_FHBSeverity`, muspC_HeadingDate = low.resp_HeadingDate, muspC_PlantHeight = low.resp_PlantHeight,
             muspC_FHBSeverity = low.resp_FHBSeverity)
  })

# Filter for the empirical crosses
# Then select only the relevant parameters and tidy it up
popvar_pred_cross <- popvar_pred$realistic %>%
  mutate(., family = str_c(4, family)) %>%
  filter(family %in% unique(vp_family_corG1$family)) %>%
  select(family, parent1, parent2, trait, note, musp_low, contains("cor"), contains("muspC"))

popvar_pred_cross_corG <- popvar_pred_cross %>%
  select(-musp_low, -contains("muspC")) %>%
  gather(trait2, prediction, contains("cor")) %>%
  mutate(trait2 = str_replace(trait2, "cor_", "")) %>%
  rename(trait1 = trait) %>%
  filter(!is.na(prediction))


# Combine the predictions with the estimates - remove NAs
popvar_pred_obs_corG <- left_join(popvar_pred_cross_corG, rename(vp_family_corG1, estimate = correlation)) %>%
  filter(!is.na(estimate))





### Measure prediction accuracy
# Calculate the correlation between predictions and observations
set.seed(242)
pred_acc <- popvar_pred_obs_corG %>% 
  group_by(trait1, trait2) %>% 
  do(cbind(bootstrap(x = .$prediction, y = .$estimate, fun = "cor", boot.reps = boot_reps, alpha = alpha), n_fam = length(.$prediction))) %>%
  rowwise() %>%
  mutate(annotation = ifelse(!between(0, ci_lower, ci_upper), "*", ""),
         trait_pair = str_c(trait1, " / ", trait2)) %>%
  ungroup()

# trait1      trait2      statistic    base    se     bias ci_lower ci_upper n_fam annotation trait_pair               
# 1 FHBSeverity HeadingDate cor        0.241  0.255 -0.0146   -0.300     0.669    14 ""         FHBSeverity / HeadingDate
# 2 FHBSeverity PlantHeight cor       -0.0119 0.296  0.0121   -0.530     0.589    14 ""         FHBSeverity / PlantHeight
# 3 HeadingDate PlantHeight cor        0.412  0.174 -0.00263   0.0239    0.711    26 *          HeadingDate / PlantHeight

