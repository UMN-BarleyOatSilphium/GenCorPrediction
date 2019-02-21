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

# Mean and range
popvar_pred_cross_corG %>% group_by(trait1, trait2) %>% summarize_at(vars(prediction), funs(mean, min, max))





popvar_pred_cross_muspC <- popvar_pred_cross %>%
  select(-contains("cor")) %>%
  gather(trait2, prediction, contains("muspC")) %>%
  mutate(trait2 = str_replace(trait2, "muspC_", "")) %>%
  rename(trait1 = trait) %>%
  filter(!is.na(prediction))


# ## Reformat the predictions using different models
# popvar_pred_model <- pred_results_model %>%
#   left_join(., cross_list, by = c("Par1" = "parent1", "Par2" = "parent2")) %>%
#   select(parent1 = Par1, parent2 = Par2, family, trait, model, pred_mu = pred.mu, pred_varG = pred.varG, musp_high = mu.sp_high,
#          musp_low = mu.sp_low, cor_HeadingDate = `cor_w/_HeadingDate`, cor_PlantHeight = `cor_w/_PlantHeight`,
#          cor_FHBSeverity = `cor_w/_FHBSeverity`) %>%
#   mutate(., family = str_c(4, family)) %>%
#   filter(family %in% unique(vp_family_corG1$family)) %>%
#   select(family, parent1, parent2, trait, model, contains("cor")) %>%
#   gather(trait2, prediction, contains("cor")) %>%
#   mutate(trait2 = str_replace(trait2, "cor_", "")) %>%
#   rename(trait1 = trait) %>%
#   filter(!is.na(prediction)) %>%
#   select(family:parent2, model, names(.))




# Combine the predictions with the estimates - remove NAs
popvar_pred_obs_corG <- left_join(popvar_pred_cross_corG, rename(vp_family_corG1, estimate = correlation)) %>%
  filter(!is.na(estimate))

# Mean and range of estimated correlations
popvar_pred_obs_corG %>% group_by(trait1, trait2) %>% summarize_at(vars(estimate), funs(mean, min, max))




popvar_pred_obs_muspC <- left_join(popvar_pred_cross_muspC, vp_family_muspC) %>%
  select(family:trait1, trait2, note, musp_estimate = musp, musp_prediction = musp_low,
         muspC_estimate = muspC, muspC_prediction = prediction) %>%
  filter(!is.na(musp_estimate))

# popvar_pred_obs_model <- left_join(popvar_pred_model, rename(vp_family_corG1, estimate = correlation)) %>%
#   filter(!is.na(estimate))








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






## Accuracy of superior progeny mean
set.seed(242)
pred_acc_muspC <- popvar_pred_obs_muspC %>% 
  group_by(trait1, trait2) %>% 
  gather(variable, value, contains("musp")) %>% 
  separate(variable, c("variable", "type"), sep = "_") %>% 
  spread(type, value) %>%
  group_by(trait1, trait2, variable) %>%
  do(cbind(bootstrap(x = .$prediction, y = .$estimate, fun = "cor", boot.reps = boot_reps, alpha = alpha), n_fam = length(.$prediction))) %>%
  rowwise() %>%
  mutate(annotation = ifelse(!between(0, ci_lower, ci_upper), "*", ""),
         trait_pair = str_c(trait1, " / ", trait2)) %>%
  ungroup()
  
  


## Compare models for prediction accuracy
set.seed(242)
pred_acc_model <- popvar_pred_obs_model %>% 
  group_by(trait1, trait2, model) %>% 
  do(cbind(bootstrap(x = .$prediction, y = .$estimate, fun = "cor", boot.reps = boot_reps, alpha = alpha), n_fam = length(.$prediction))) %>%
  rowwise() %>%
  mutate(annotation = ifelse(!between(0, ci_lower, ci_upper), "*", ""),
         trait_pair = str_c(trait1, " / ", trait2)) %>%
  ungroup()




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
        ylim(c(-0.75, 0.75)) +
        facet_grid(~ trait_pair, scales = "free_x") + 
        theme_acs() )



# Save the plots
for (i in seq_along(g_pred_acc)) {
  filename <- str_c(names(g_pred_acc)[i], "_corG_pred_acc.jpg")
  ggsave(filename = filename, plot = g_pred_acc[[i]], path = fig_dir, height = 2, width = 5, dpi = 1000)
}




## Plot the realistic results in presentation format
pred_acc_realistic_annotation <- popvar_pred_obs %>%
  filter(tp_set == "realistic") %>% 
  distinct(trait1, trait2, tp_set) %>% 
  left_join(., pred_acc) %>%
  mutate(annotation = str_c("r = ", round(base, 2), annotation))

g_pred_acc_realistic <- popvar_pred_obs %>%
  filter(tp_set == "realistic") %>%
  mutate(trait_pair = str_c(trait1, " / ", trait2)) %>%
  ggplot(aes(x = prediction, y = estimate)) +
  geom_smooth(method = "lm", se = FALSE) + 
  geom_point() + 
  geom_text(data = pred_acc_realistic_annotation, aes(x = Inf, y = -Inf, label = annotation), size = 3, hjust = 1.2, vjust = -1) + 
  ylab("Observation") +
  xlab("Prediction") + 
  ylim(c(-0.75, 0.75)) +
  facet_grid(~ trait_pair, scales = "free_x") + 
  theme_presentation2()


ggsave(filename = "realistic_corG_pred_acc_presentation.jpg", plot = g_pred_acc_realistic, path = fig_dir,
       height = 3.5, width = 8, dpi = 1000)



# Save the plots
for (i in seq_along(g_pred_acc)) {
  filename <- str_c(names(g_pred_acc)[i], "_corG_pred_acc.jpg")
  ggsave(filename = filename, plot = g_pred_acc[[i]], path = fig_dir, height = 2, width = 5, dpi = 1000)
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

