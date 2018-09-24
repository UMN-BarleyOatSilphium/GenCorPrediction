## PopVarValidation - prediction accuracy and analysis 
##
## This script will look at the prediction accuracy of family correlation
## 
## Author: Jeff Neyhart
## Last modified: September 19, 2018
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
load(file.path(result_dir, "vp_family_correlation_analysis.RData"))


## First subset the relevant columns
popvar_pred <- list(pred_results_realistic, pred_results_relevant) %>% 
  setNames(c("realistic", "relevant")) %>%
  map(~{
    filter(., trait %in% traits) %>%
      left_join(., cross_list, by = c("Par1" = "parent1", "Par2" = "parent2")) %>%
      select(parent1 = Par1, parent2 = Par2, family, trait, note, pred_mu = pred.mu, pred_varG = pred.varG, musp_high = mu.sp_high,
             musp_low = mu.sp_low, cor_HeadingDate = `cor_w/_HeadingDate`, cor_PlantHeight = `cor_w/_PlantHeight`,
             cor_FHBSeverity = `cor_w/_FHBSeverity`)
  })

# Filter for the empirical crosses
# Then select only the relevant parameters and tidy it up
popvar_pred_cross <- popvar_pred %>%
  map(~mutate(., family = str_c(4, family)) %>%
        filter(family %in% unique(vp_family_corG1$family)) %>%
        select(family, parent1, parent2, trait, note, contains("cor")) %>%
        gather(trait2, prediction, contains("cor")) %>%
        mutate(trait2 = str_replace(trait2, "cor_", "")) %>%
        rename(trait1 = trait) %>%
        filter(!is.na(prediction))) %>%
  list(., names(.)) %>%
  pmap(~{names(.x)[ncol(.x)] <- .y; .x}) %>% 
  reduce(left_join) %>%
  gather(tp_set, prediction, realistic, relevant)

## Edit the bootstrapping results
vp_family_results <- vp_family_boot %>% 
  mutate(parameter = case_when(parameter == "mu" ~ "family_mean", parameter == "varG" ~ "variance", TRUE ~ "mu_sp"),
         expectation = base + bias) %>%
  select(trait, family, parameter, estimate = base, expectation)

# Combine the predictions with the estimates - remove NAs
popvar_pred_obs <- left_join(popvar_pred_cross, rename(vp_family_corG1, estimate = correlation)) %>%
  filter(!is.na(estimate))








### Measure prediction accuracy
# Calculate the correlation between predictions and observations
set.seed(242)
pred_acc <- popvar_pred_obs %>% 
  group_by(trait1, trait2, tp_set) %>% 
  do(cbind(bootstrap(x = .$prediction, y = .$estimate, fun = "cor", boot.reps = boot_reps, alpha = alpha), n_fam = length(.$prediction))) %>%
  rowwise() %>%
  mutate(annotation = ifelse(!between(0, ci_lower, ci_upper), "*", ""),
         trait_pair = str_c(trait1, " / ", trait2)) %>%
  ungroup()

# ## Do predictions of the variance improve when using the unbiased expectation?
# pred_acc_exp <- popvar_pred_obs %>%
#   group_by(trait, parameter, tp_set) %>%
#   do(bootstrap(x = .$prediction, y = .$expectation, fun = "cor", boot.reps = boot_reps, alpha = alpha)) %>%
#   rowwise() %>%
#   mutate(annotation = ifelse(!between(0, ci_lower, ci_upper), "*", "")) %>%
#   ungroup()


# trait1      trait2      tp_set    statistic    base    se      bias ci_lower ci_upper n_fam annotation
# 1 FHBSeverity HeadingDate realistic cor        0.241  0.262 -0.000556  -0.287     0.692    14 ""        
# 2 FHBSeverity HeadingDate relevant  cor        0.191  0.263  0.000722  -0.353     0.658    14 ""        
# 3 FHBSeverity PlantHeight realistic cor       -0.0119 0.275  0.00431   -0.516     0.522    14 ""        
# 4 FHBSeverity PlantHeight relevant  cor        0.0781 0.239  0.00372   -0.408     0.515    14 ""        
# 5 HeadingDate PlantHeight realistic cor        0.412  0.174 -0.00650    0.0501    0.699    26 *         
# 6 HeadingDate PlantHeight relevant  cor        0.160  0.202  0.00515   -0.246     0.539    26 ""

## Just the realistic training set
filter(pred_acc, tp_set == "realistic")

# trait1      trait2      tp_set    statistic    base    se      bias ci_lower ci_upper n_fam annotation
# 1 FHBSeverity HeadingDate realistic cor        0.241  0.262 -0.000556  -0.287     0.692    14 ""        
# 2 FHBSeverity PlantHeight realistic cor       -0.0119 0.275  0.00431   -0.516     0.522    14 ""        
# 3 HeadingDate PlantHeight realistic cor        0.412  0.174 -0.00650    0.0501    0.699    26 * 


# ## Create and write a table
# pred_table <- pred_acc %>% 
#   mutate(annotation = str_c(round(base, 2), " (", round(ci_lower, 2), ", ", round(ci_upper, 2), ")")) %>% 
#   select(trait, n_fam, tp_set, parameter, annotation) %>% 
#   mutate(parameter = str_replace_all(parameter, param_replace), 
#          parameter = factor(parameter, levels = param_replace)) %>% 
#   spread(parameter, annotation) %>%
#   split(.$tp_set)
# 
# # for (df in pred_table) {
# #   write_csv(x = select(df, -tp_set), path = file.path(fig_dir, str_c("pred_accuracy_", unique(df$tp_set), ".csv")))
# # }




# Plot the same results
g_pred_acc <- popvar_pred_obs %>%
  mutate(trait_pair = str_c(trait1, " / ", trait2)) %>%
  split(.$tp_set) %>%
  map(~ggplot(., aes(x = prediction, y = estimate)) +
        geom_smooth(method = "lm", se = FALSE) + 
        geom_point() + 
        geom_text(data = mutate(subset(pred_acc, tp_set == unique(.$tp_set)), annotation = str_c("r = ", round(base, 2), annotation)),
                  aes(x = Inf, y = -Inf, label = annotation), size = 3, hjust = 1.2, vjust = -1) + 
        ylab("Observation") +
        xlab("Prediction") + 
        facet_wrap(~ trait_pair, ncol = 3, scales = "free") + 
        theme_acs() )



# Save the plots
for (i in seq_along(g_pred_acc)) {
  filename <- str_c(names(g_pred_acc)[i], "_corG_pred_acc.jpg")
  ggsave(filename = filename, plot = g_pred_acc[[i]], path = fig_dir, height = 2, width = 6, dpi = 1000)
}







## Analyze bias
# First calculate bias on a per-family basis
#### This needs to be rechecked! ####
popvar_bias <- popvar_pred_obs %>% 
  # mutate_at(vars(prediction, estimate), abs) %>%
  mutate(bias = (prediction - estimate) / abs(estimate)) %>%
  # filter(tp_set == "realistic") %>%
  select(family, trait1, trait2, tp_set, prediction:bias)

# Next summarize over traits
popvar_trait_bias <- popvar_bias %>% 
  group_by(tp_set, trait1, trait2) %>% 
  summarize(bias = mean(prediction - estimate) / mean(estimate))

# Plot
popvar_bias %>%
  ggplot(aes(x = family, y = bias)) + 
  geom_point() + 
  facet_wrap(~ trait1 + trait2, scales = "free")

