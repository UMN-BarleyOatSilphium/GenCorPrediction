## PopVarValidation - genetic correlation simulation analysis
## 
## Author: Jeff Neyhart
## Last modified: September 3, 2018
## 

# Load the source script
repo_dir <- getwd()
source(file.path(repo_dir, "source.R"))

# Other libraries
library(cowplot)

# Significance level
alpha <- 0.05

# Create a replacement vector
arch_replace <- c(loose_link = "Loose Linkage", close_link = "Close Linkage", pleio = "Pleiotropy")




# Load the simulation results
load(file.path(result_dir, "popvar_gencor_simulation_prediction_results.RData"))


## Are there any missing combinations?
popvar_prediction_simulation_out %>%
  distinct(trait1_h2, trait2_h2, nQTL, tp_size, gencor, arch, model, iter) %>%
  mutate_all(as.factor) %>% 
  mutate(obs = T) %>% 
  complete(trait1_h2, trait2_h2, nQTL, tp_size, gencor, arch, model, iter, fill = list(obs = F)) %>% 
  filter(!obs)

pred_sim_tidy <- popvar_prediction_simulation_out %>%
  bind_cols(., as_data_frame(transpose(popvar_prediction_simulation_out$results))) %>%
  select(-results) %>%
  mutate_at(vars(trait1_h2, trait2_h2, nQTL, tp_size, gencor, arch, model), as.factor) %>%
  mutate(arch = factor(str_replace_all(arch, arch_replace), levels = arch_replace))


## Tidy the results
sim_summary_tidy <- pred_sim_tidy %>% 
  unnest(summary)

sim_meta_tidy <- pred_sim_tidy %>% 
  unnest(other)



## What are the genetic correlations in the tp for each tp size and intented gencor?

# Plot
g_base_corG <- sim_meta_tidy %>% 
  filter(variable == "tp_gencor", tp_size == 600) %>%
  group_by(arch, gencor, nQTL) %>%
  mutate(n = n()) %>%
  ungroup() %>%
  {ggplot(data = ., aes(x = value, fill = gencor)) + 
  geom_density(alpha = 0.5) + 
  facet_grid(nQTL ~ arch) +
  xlim(c(-1, 1)) + 
  labs(caption = paste0("N_TP = ", unique(.$tp_size), ", n = ", unique(.$n))) +
  theme_acs()}

ggsave(filename = "gencor_base_corG.jpg", plot = g_base_corG, path = fig_dir, height = 4, width = 5, dpi = 1000)


## Summarize
sim_out_summ <- sim_summary_tidy %>%
  gather(variable, value, accuracy, bias) %>% 
  group_by(trait1_h2, trait2_h2, nQTL, tp_size, gencor, arch, model, trait, parameter, variable) %>% 
  summarize_at(vars(value), funs(mean(., na.rm = T), sd(., na.rm = T), n())) %>% 
  ungroup() %>%
  mutate(stat = qt(p = 1 - (alpha / 2), df = n - 1) * (sd / sqrt(n) ),
         lower = mean - stat, upper = mean + stat,
         round_mean = round(mean, 2))





### Fit a model for corG
fit <- lm(accuracy ~ trait1_h2 + trait2_h2 + nQTL + tp_size + gencor + arch + model + gencor:arch + trait1_h2:trait2_h2 + 
            model:arch + model:nQTL + model:nQTL:arch, 
          data = sim_summary_tidy, subset = parameter == "corG")
anova(fit)

# Effect plot
effs <- effects::allEffects(fit)
## Convert to df
effs_df <- map(effs, as.data.frame)

# Plot the nQTL:arch:model results
g_arch_accuracy <- effs_df$`nQTL:arch:model` %>%
  mutate(arch = factor(arch, levels = arch_replace)) %>%
  ggplot(aes(x = arch, y = fit, ymin = lower, ymax = upper, color = model, shape = nQTL)) + 
  geom_point(position = position_dodge(0.5)) + 
  geom_errorbar(position = position_dodge(0.5), width = 0.5) +
  scale_color_discrete(name = "Model") +
  xlab("Genetic architecture") +
  ylab("Prediction accuracy") +
  theme_acs() +
  theme(legend.position = c(0.15, 0.25), legend.key.height = unit(0.75, "line"))

ggsave(filename = "gencor_accuracy_architecture.jpg", plot = g_arch_accuracy, path = fig_dir,
       height = 3, width = 3.5, dpi = 1000)


## TP size breaks for the scale
tp_size_breaks <- parse_number(unique(sim_out_summ$tp_size))


## Plot using points
g_accuracy <- sim_out_summ %>% 
  filter(parameter == "corG", variable == "accuracy") %>%
  mutate(trait2_h2 = factor(trait2_h2, levels = rev(levels(trait2_h2))),
         tp_size = parse_number(tp_size)) %>%
  ggplot(aes(x = tp_size, y = mean, ymin = lower, ymax = upper, color = arch)) + 
  geom_point(size = 1) +
  geom_line() + 
  geom_errorbar(width = 0.25) +
  scale_color_discrete(name = NULL) +
  scale_x_continuous(breaks = tp_size_breaks) + 
  facet_grid(trait2_h2 + trait1_h2 ~ nQTL + gencor + model, labeller = label_both) +
  theme_acs() +
  theme(legend.position = c(0.85, 0.75), legend.key.height = unit(0.5, "lines"))



## Plot using points
g_accuracy1 <- sim_out_summ %>% 
  # filter(parameter == "corG", variable == "accuracy", gencor == 0.5, nQTL == 100) %>%
  filter(parameter == "corG", variable == "accuracy", gencor == 0.5, nQTL == 30) %>%
  mutate(trait2_h2 = factor(trait2_h2, levels = rev(levels(trait2_h2))),
         tp_size = parse_number(tp_size)) %>%
  {ggplot(data = ., aes(x = tp_size, y = mean, ymin = lower, ymax = upper, color = arch, lty = model)) + 
  geom_point(size = 1) +
  geom_line() + 
  geom_errorbar(width = 25) +
  scale_color_discrete(name = NULL) +
  scale_linetype_discrete(name = NULL) + 
  scale_x_continuous(breaks = tp_size_breaks) + 
  facet_grid(trait2_h2 ~ trait1_h2, labeller = label_both) +
  labs(caption = paste0("Base corG: ", unique(.$gencor), ", nQTL: ", unique(.$nQTL))) +
  theme_acs() +
  theme(legend.position = c(0.85, 0.75), legend.key.height = unit(0.5, "lines"))}

# Save
ggsave(filename = "gencor_simulation_accuracy_model30.jpg", plot = g_accuracy1, path = fig_dir,
       height = 5, width = 5, dpi = 1000)

## Plot using points
g_accuracy2 <- sim_out_summ %>% 
  filter(parameter == "corG", variable == "accuracy", gencor == 0.5, nQTL == 100,
         model == "RRBLUP") %>%
  mutate(trait2_h2 = factor(trait2_h2, levels = rev(levels(trait2_h2))),
         tp_size = parse_number(tp_size)) %>%
  ggplot(aes(x = tp_size, y = mean, ymin = lower, ymax = upper, color = arch)) + 
  geom_point(size = 1) +
  geom_line() + 
  geom_errorbar(width = 25) +
  scale_color_discrete(name = NULL) +
  scale_linetype_discrete(name = NULL) + 
  scale_x_continuous(breaks = tp_size_breaks) + 
  facet_grid(trait2_h2 ~ trait1_h2, labeller = label_both) +
  labs(caption = "Base corG: 0.5, nQTL: 100") +
  theme_acs() +
  theme(legend.position = c(0.85, 0.75), legend.key.height = unit(0.5, "lines"))

# Save
ggsave(filename = "gencor_simulation_accuracy.jpg", plot = g_accuracy2, path = fig_dir,
       height = 5, width = 5, dpi = 1000)





# ## Bias
# g_bias <- sim_out_summ %>% 
#   filter(parameter == "corG", variable == "bias") %>%
#   ggplot(aes(x = trait1_h2, y = trait2_h2, fill = mean, label = round_mean)) + 
#   geom_tile() + 
#   # geom_text(size = 2) + 
#   scale_fill_gradient2(low = "red", mid = "white", high = "blue") +
#   facet_grid(arch + gencor ~ tp_size + nQTL, labeller = label_both) +
#   theme_acs()
# 
# # Save
# ggsave(filename = "gen_cor_simulation_bias.jpg", plot = g_bias, path = fig_dir,
#        height = 3, width = 5, dpi = 1000)













### Genetic correlation genetic architecture simulation
# Load the simulation results
load(file.path(result_dir, "popvar_gencor_space_simulation_results.RData"))

# Mutate the architecture space combinations
sim_out1 <- popvar_corG_space_simulation_out %>% 
  mutate(probcor = map(probcor, ~`names<-`(as.data.frame(.), c("dLinkage", "pLinkage")) %>% tail(., 1))) %>% # The tail is used to remove the probabilities of pleiotropy))
  unnest(probcor)
  
## Are there any missing combinations?
sim_out1 %>%
  distinct(trait1_h2, trait2_h2, gencor, dLinkage, pLinkage, iter) %>%
  mutate_all(as.factor) %>% 
  mutate(obs = T) %>% 
  complete(trait1_h2, trait2_h2, gencor, dLinkage, pLinkage, iter, fill = list(obs = F)) %>% 
  filter(!(dLinkage == 0 & pLinkage != 0)) %>%
  filter(!obs)

## Good!

# Tidy the results and extract the probability of linkage and the degree of linkage
sim_results_tidy <- sim_out1 %>%
  ## Extract the results
  bind_cols(., as_data_frame(transpose(.$results))) %>%
  select(-results) %>%
  mutate_at(vars(trait1_h2:gencor, dLinkage, pLinkage), as.factor)
  

## Extract each dataset
correlation_tidy <- unnest(sim_results_tidy, other)
predictions_tidy <- unnest(sim_results_tidy, summary)


## Extract the training population genetic correlation
base_cor_summ <- correlation_tidy %>%
  group_by(trait1_h2, trait2_h2, nQTL, tp_size, gencor, dLinkage, pLinkage, variable) %>%
  summarize_at(vars(value), funs(mean(., na.rm = T), sd(., na.rm = T), n())) %>%
  ## Fill-in missing combinations when dLinkage == 0
  bind_rows(., 
            filter(., dLinkage == 0) %>%
            group_by(variable, add = T) %>% 
            complete(trait1_h2, trait2_h2, nQTL, tp_size, gencor, pLinkage, variable) %>%
            group_by(trait1_h2, trait2_h2, nQTL, tp_size, gencor, variable) %>% 
            mutate_at(vars(mean, sd), funs(mean), na.rm = T) %>%
            ungroup() %>%
            distinct(trait1_h2, trait2_h2, nQTL, tp_size, gencor, pLinkage, dLinkage, variable, mean, sd) %>%
            filter(pLinkage != 1)
  ) %>% ungroup()





## Plot
g_base_cor <- base_cor_summ %>%
  filter(variable == "tp_gencor") %>%
  ggplot(aes(x = pLinkage, y = dLinkage, fill = mean)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", high = "red") +
  facet_wrap(~ gencor, nrow = 1) +
  theme_acs()

# save
ggsave(filename = "gencor_arch_space_base_corG.jpg", plot = g_base_cor, path = fig_dir,
       height = 3, width = 6, dpi = 1000)



  
## Extract the prediction results
# Summarize
pred_results_summ <- predictions_tidy %>%
  gather(variable, value, accuracy, bias) %>%
  filter(!(variable == "bias" & abs(value) > 2)) %>%
  group_by(trait1_h2, trait2_h2, nQTL, tp_size, gencor, dLinkage, pLinkage, trait, parameter, variable) %>%
  summarize_at(vars(value), funs(mean(., na.rm = T), sd(., na.rm = T), n())) %>%  ## Fill-in missing combinations when dLinkage == 0
  bind_rows(., 
            filter(., dLinkage == 0) %>%
              group_by(variable, add = T) %>% 
              complete(trait1_h2, trait2_h2, nQTL, tp_size, gencor, pLinkage, parameter, variable) %>%
              group_by(trait1_h2, trait2_h2, nQTL, tp_size, gencor, parameter, variable) %>% 
              mutate_at(vars(mean, sd), funs(mean), na.rm = T) %>%
              ungroup() %>%
              distinct(trait1_h2, trait2_h2, nQTL, tp_size, gencor, pLinkage, dLinkage, parameter, variable, mean, sd) %>%
              filter(pLinkage != 1)
  ) %>% ungroup()


## Plot results for genetic correlation
g_pred_acc_corG <- pred_results_summ %>%
  filter(parameter == "corG", variable == "accuracy") %>%
  ggplot(aes(x = pLinkage, y = dLinkage, fill = mean)) +
  geom_tile() +
  scale_fill_gradient(limits = c(0.40, 0.72), low = "white", high = "green", name = "Prediction\naccuracy") +
  facet_grid(~ gencor) +
  ylab("Maximum distance between QTL (cM)") +
  xlab("Proportion of non-pleiotropic QTL") +
  theme_acs()


# bias
g_pred_bias_corG <- pred_results_summ %>%
  filter(parameter == "corG", variable == "bias") %>%
  ggplot(aes(x = pLinkage, y = dLinkage, fill = mean)) +
  geom_tile() +
  scale_fill_gradient2(low = "red", high = "blue", name = "Bias") +
  facet_grid(~ gencor) +
  ylab("Maximum distance between QTL (cM)") +
  xlab("Proportion of non-pleiotropic QTL") +
  theme_acs()


# Combine
g_pred_corG <- plot_grid(
  g_pred_acc_corG + theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()), 
  g_pred_bias_corG, 
  ncol = 1, align = "hv")

ggsave(filename = "gencor_arch_space_pred.jpg", plot = g_pred_corG, path = fig_dir, width = 6, height = 5, dpi = 1000)




# ## Fit a model for correlation
# fit <- lm(value ~ gencor + dLinkage + pLinkage + gencor:pLinkage, data = correlation_tidy, subset = variable == "tp_gencor" & dLinkage != 0)
# anova(fit)
# plot(effects::allEffects(fit))
# 
# ## Notes
# ## 1. It works
# 
predictions_tidy_tomodel <- predictions_tidy %>%
  mutate(pLinkage = parse_number(pLinkage)) %>%
  filter(parameter == "corG", dLinkage != 0)



# Treat pLinkage as numeric
fit <- lm(accuracy ~ gencor + dLinkage + pLinkage + gencor:pLinkage, data = predictions_tidy_tomodel)
anova(fit)
plot(effects::allEffects(fit))

## Note
## 1. Upward trend in prediction accuracy with decreasing pleiotropy, as expected.

effs_df <- map(effects::allEffects(fit), as.data.frame)

# Plot
g_pLinkage_accuracy <- effs_df$`gencor:pLinkage` %>%
  # mutate(pLinkage = parse_number(pLinkage)) %>%
  ggplot(aes(x = pLinkage, y = fit, ymin = lower, ymax = upper, fill = gencor)) +
  # geom_point() +
  geom_line() +
  geom_ribbon(alpha = 0.2) +
  scale_fill_discrete(name = "Genetic\ncorrelation") +
  ylab("Prediction accuracy") +
  xlab("Proportion of non-pleiotropic QTL") +
  theme_acs()

ggsave(filename = "gencor_plinkage_accuracy.jpg", plot = g_pLinkage_accuracy, path = fig_dir, width = 3, height = 3, dpi = 1000)





















### Genetic correlation selection simulation


# Load the simulation results
load(file.path(result_dir, "popvar_gencor_selection_simulation_results.RData"))

## Are there any missing trials?
popvar_gencor_selection_simulation_out %>%
  distinct(trait1_h2, trait2_h2, gencor, arch, iter) %>%
  mutate_all(as.factor) %>% 
  mutate(obs = T) %>% 
  complete(trait1_h2, trait2_h2, gencor, arch, iter, fill = list(obs = F)) %>% 
  filter(!obs) %>%
  nrow()



# Tidy
popvar_selection_sim_tidy <- popvar_gencor_selection_simulation_out %>% 
  bind_cols(., as_data_frame(transpose(popvar_gencor_selection_simulation_out$results))) %>% 
  select(-input, -results) %>%
  mutate_at(vars(trait1_h2, trait2_h2, gencor, arch), as.factor) %>%
  mutate(arch = factor(str_replace_all(arch, arch_replace), levels = arch_replace))

# Get a vector of the selection types
selections <- popvar_selection_sim_tidy %>% unnest(response) %>% pull(selection) %>% unique()


## Unnest the correlations
correlations_tidy <- popvar_selection_sim_tidy %>%
  unnest(correlations) %>%
  filter(variable %in% c("cor", "cov")) %>% 
  mutate(data = map(value, ~{
    if (is.list(.)) {
      data_frame(selection = names(.), value = unlist(.))
    } else {
      data_frame(selection = selections, value = .)
    } })) %>%
  unnest(data)

## Unnest the variance
variance_tidy <- popvar_selection_sim_tidy %>%
  unnest(correlations) %>%
  filter(variable == "var") %>% 
  mutate(data = map(value, ~{
    if (is.list(.)) {
      data_frame(selection = names(.), trait1 = map_dbl(., "trait1"), trait2 = map_dbl(., "trait2"))
    } else {
      data_frame(selection = selections, trait1 = .["trait1"], trait2 = .["trait2"])
    } })) %>%
  unnest(data) %>%
  gather(trait, value, trait1, trait2)



# Split into two parts
response_tidy <- popvar_selection_sim_tidy %>%
  unnest(response) %>%
  mutate(intensity = round(intensity, 2))

# Response for each trait
response_tidy_part1 <- response_tidy %>%
  select(trait1_h2:trait, relative_response, stand_response) %>% 
  gather(variable, value, relative_response, stand_response)

# Response for the trait index
response_tidy_part2 <- response_tidy %>% 
  select(trait1_h2:trait, contains("index")) %>% 
  filter(trait == "trait1") %>% 
  gather(variable, value, relative_response_index, stand_response_index)

# Genetic correlation and variance in the selected C1
response_tidy_corG <- response_tidy %>%
  select(trait1_h2:trait, corG, var) %>%
  group_by(trait1_h2, trait2_h2, gencor, arch, iter, intensity, selection) %>%
  mutate(cov = corG[1] * prod(sqrt(var))) %>%
  ungroup() %>%
  rename(cor = corG) %>%
  gather(variable, value, cor, cov, var)


## Combine with the previous tidy correlations
correlations_tidy1 <- bind_rows(
  correlations_tidy,
  filter(response_tidy_corG, trait == "trait1", variable != "var") %>% 
    mutate(type = "c1_select") %>% select(-trait)
)

variance_tidy1 <- bind_rows(
  variance_tidy,
  filter(response_tidy_corG, variable == "var") %>% mutate(type = "c1_select")
)



## Plot the genetic correlation in the TP
correlations_tidy1 %>% 
  filter(type == "tp_base", variable == "cor", selection == "mean") %>% 
  ggplot(aes(x = value, fill = gencor)) + 
  geom_density(alpha = 0.5) + 
  facet_grid( ~ arch) +
  theme_acs()



## We have a problem with the close linkage and loose linkage groups.
## Problem may have been solved



## Combine and summarize the response results
response_tomodel <- bind_rows(response_tidy_part1, response_tidy_part2)
  
  
## Summarize over iterations
response_summary <- response_tomodel %>% 
  filter(!is.na(value)) %>%
  group_by(trait1_h2, trait2_h2, gencor, arch, trait, intensity, selection, variable) %>%
  summarize_at(vars(value), funs(mean, sd, n())) %>%
  ungroup() %>%
  mutate(stat = qt(p = 1 - (alpha / 2), df = n - 1) * (sd / sqrt(n) ),
         lower = mean - stat, upper = mean + stat)
    

### Plotting 

# Choose a selection intensity to investigate
i_sp <- 0.10


# Response to selection relative to that based on the mean
vrbl <- "stand_response_index"

## Response to selection of the index - relative to base population
g_stand_resp <- response_summary %>%
  filter(variable == vrbl) %>%
  ggplot(aes(x = intensity, y = mean, shape = selection, fill = selection)) +
  geom_point(size = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  ylab("Standardized response (relative to base population)") +
  xlab("Proportion of selected individuals") +
  facet_grid(trait1_h2 + trait2_h2 ~ arch + gencor, 
             labeller = labeller(arch = label_value, .default = label_both)) +
  theme_acs()

# Save
ggsave(filename = "gencor_sim_stand_response.jpg", plot = g_stand_resp, path = fig_dir,
       height = 12, width = 12, dpi = 1000)



## Take a subset
g_stand_resp_sub <- response_summary %>%
  filter(variable == vrbl) %>%
  # filter(trait1_h2 == 1, gencor != 0) %>%
  filter(trait1_h2 == 0.6, trait2_h2 != 1, selection != "musp") %>%
  ggplot(aes(x = intensity, y = mean, shape = selection, fill = selection)) +
  geom_point(size = 0.5) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  ylab("Standardized response (relative to base population)") +
  xlab("Proportion of selected individuals") +
  facet_grid(trait1_h2 + trait2_h2 ~ arch + gencor, 
             labeller = labeller(arch = label_value, .default = label_both)) +
  theme_acs()

# Save
ggsave(filename = "gencor_sim_stand_response_subset.jpg", plot = g_stand_resp_sub, path = fig_dir,
       height = 3.5, width = 10, dpi = 1000)



## Highlight one selection intensity
g_stand_resp_isp <- response_summary %>%
  filter(variable == "stand_response_index", intensity == i_sp) %>%
  ggplot(aes(x = arch, y = mean, fill = selection)) +
  geom_col(position = "dodge") +
  geom_errorbar(aes(ymin = lower, ymax = upper), position = position_dodge(0.9), width = 0.5) +
  ylab("Standardized response (relative to base population)") +
  xlab("Genetic correlation architecture") +
  facet_grid(trait1_h2 + trait2_h2 ~ gencor, 
             labeller = labeller(arch = label_value, .default = label_both)) +
  theme_acs() +
  labs(subtitle = paste("i =", i_sp))

# Save
ggsave(filename = "gencor_sim_stand_response_isp.jpg", plot = g_stand_resp_isp, path = fig_dir,
       height = 12, width = 6, dpi = 1000)


## Highlight one selection intensity
g_stand_resp_isp_sub <- response_summary %>%
  filter(variable == "stand_response_index", intensity == i_sp) %>%
  # filter(trait1_h2 == 1) %>%
  filter(trait1_h2 == 0.6, trait2_h2 != 1, selection != "musp") %>%
  ggplot(aes(x = arch, y = mean, fill = selection)) +
  geom_col(position = "dodge") +
  geom_errorbar(aes(ymin = lower, ymax = upper), position = position_dodge(0.9), width = 0.5) +
  ylab("Standardized response (relative to base population)") +
  xlab("Genetic correlation architecture") +
  facet_grid(trait1_h2 + trait2_h2 ~ gencor, 
             labeller = labeller(arch = label_value, .default = label_both)) +
  theme_acs() +
  labs(subtitle = paste("i =", i_sp))

# Save
ggsave(filename = "gencor_sim_stand_response_isp_subset.jpg", plot = g_stand_resp_isp_sub, path = fig_dir,
       height = 4, width = 6, dpi = 1000)


## Model
fit_resp <- response_tomodel %>% 
  filter(variable == "stand_response_index", intensity %in% i_sp) %>%
  mutate_if(is.character, as.factor) %>%
  mutate(intensity = as.factor(intensity)) %>%
  lm(value ~ trait1_h2 + trait2_h2 + gencor + arch + selection + gencor:arch + selection:arch, data = .)

anova(fit_resp)
plot(effects::allEffects(fit_resp))

## Notes:
## 1. There is a strong interaction between architecture and intended genetic correlation, where pleiotropic architecture
## is highly influenced by the genetic correlation
## 2. There is moderate interaction of architecture and selection strategy, suggesting that accounting for genetic correlation
## is better when the architecture is not pleiotropic

## Fit a model without random and corG selection
## Model
fit_resp <- response_tomodel %>% 
  filter(variable == "stand_response_index", intensity %in% i_sp, selection %in% c("mean", "muspC")) %>%
  mutate_if(is.character, as.factor) %>%
  mutate(intensity = as.factor(intensity)) %>%
  lm(value ~ trait1_h2 + trait2_h2 + gencor + arch + selection + gencor:arch + selection:arch, data = .)

anova(fit_resp)
plot(effects::allEffects(fit_resp))







# Response to selection relative to that based on the mean
vrbl <- "relative_response_index"

g_relative_resp <- response_summary %>%
  filter(variable == vrbl) %>%
  ggplot(aes(x = intensity, y = mean, shape = selection, fill = selection)) +
  geom_point(size = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  ylab("Standardized response (relative to mean selection)") +
  xlab("Proportion of selected individuals") +
  facet_grid(trait1_h2 + trait2_h2 ~ arch + gencor, 
             labeller = labeller(arch = label_value, .default = label_both)) +
  theme_acs()

# Save
ggsave(filename = "gencor_sim_relative_response.jpg", plot = g_relative_resp, path = fig_dir,
       height = 12, width = 12, dpi = 1000)

# Response relative to selection on the mean - subset
g_relative_resp_sub <- response_summary %>%
  filter(variable == vrbl) %>%
  # filter(trait1_h2 == 1, gencor != 0) %>%
  filter(trait1_h2 == 0.6, trait2_h2 != 1, selection != "musp") %>%
  ggplot(aes(x = intensity, y = mean, shape = selection, fill = selection)) +
  geom_hline(yintercept = 0, size = 0.25) + 
  geom_point(size = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  ylab("Standardized response (relative to mean selection)") +
  xlab("Proportion of selected individuals") +
  facet_grid(trait1_h2 + trait2_h2 ~ arch + gencor, 
             labeller = labeller(arch = label_value, .default = label_both)) +
  theme_acs() +
  theme(axis.text.x = element_text(size = 5))

# Save
ggsave(filename = "gencor_sim_relative_response_sub.jpg", plot = g_relative_resp_sub, path = fig_dir,
       height = 3.5, width = 10, dpi = 1000)



## Highlight one selection intensity
g_relative_resp_isp <- response_summary %>%
  filter(variable == "relative_response_index", intensity == i_sp, selection != "mean") %>%
  ggplot(aes(x = arch, y = mean, fill = selection)) +
  geom_col(position = "dodge") +
  geom_errorbar(aes(ymin = lower, ymax = upper), position = position_dodge(0.9), width = 0.5) +
  ylab("Standardized response (relative to base population)") +
  xlab("Genetic correlation architecture") +
  facet_grid(trait1_h2 + trait2_h2 ~ gencor, 
             labeller = labeller(arch = label_value, .default = label_both)) +
  theme_acs() +
  labs(subtitle = paste("i =", i_sp))

# Save
ggsave(filename = "gencor_sim_relative_response_isp.jpg", plot = g_relative_resp_isp, path = fig_dir,
       height = 12, width = 6, dpi = 1000)


## Highlight one selection intensity
g_relative_resp_isp_sub <- response_summary %>%
  filter(variable == vrbl, intensity == i_sp, selection != "mean") %>%
  # filter(trait1_h2 == 1) %>%
  filter(trait1_h2 == 0.6, trait2_h2 != 1, selection != "musp") %>%
  ggplot(aes(x = arch, y = mean, fill = selection)) +
  geom_col(position = "dodge") +
  geom_errorbar(aes(ymin = lower, ymax = upper), position = position_dodge(0.9), width = 0.5) +
  ylab("Standardized response (relative to base population)") +
  xlab("Genetic correlation architecture") +
  facet_grid(trait1_h2 + trait2_h2 ~ gencor, 
             labeller = labeller(arch = label_value, .default = label_both)) +
  theme_acs() +
  labs(subtitle = paste("i =", i_sp))

# Save
ggsave(filename = "gencor_sim_relative_response_isp_subset.jpg", plot = g_relative_resp_isp_sub, path = fig_dir,
       height = 5, width = 6, dpi = 1000)


## Model
fit_resp <- response_tomodel %>% 
  filter(variable == "relative_response_index", intensity %in% i_sp, selection != "mean") %>%
  mutate_if(is.character, as.factor) %>%
  mutate(intensity = as.factor(intensity)) %>%
  lm(value ~ trait1_h2 + trait2_h2 + gencor + arch + selection + gencor:arch + selection:arch, data = .)

anova(fit_resp)
plot(effects::allEffects(fit_resp))

## Notes
## 1. The advantage of selecting on musp or muspC is greater when heritability is lower


## Model
fit_resp <- response_tomodel %>% 
  filter(variable == "relative_response_index", intensity %in% i_sp, selection %in% c("muspC")) %>%
  mutate_if(is.character, as.factor) %>%
  mutate(intensity = as.factor(intensity)) %>%
  lm(value ~ trait1_h2 + trait2_h2 + gencor + arch + gencor:arch, data = .)

anova(fit_resp)
plot(effects::allEffects(fit_resp))

## Notes
## 1. architecture x selection is not significant
## 2. musp and muspC outperform mu, especially if the architecture is not pleiotropy
## 





### Examine the genetic correlations
correlations_tidy2 <- correlations_tidy1 %>% 
  filter(is.na(intensity) | intensity == i_sp)

correlations_summ <- correlations_tidy2 %>%
  group_by(trait1_h2, trait2_h2, gencor, arch, selection, variable, type) %>%   
  summarize_at(vars(value), funs(mean(., na.rm = T), sd(., na.rm = T), n())) %>%
  ungroup() %>%
  mutate(stat = qt(p = 1 - (alpha / 2), df = n - 1) * (sd / sqrt(n) ),
         lower = mean - stat, upper = mean + stat,
         type = factor(type, c("tp_base", "tp_select", "c1_all", "c1_select")))


# variance_tidy2 <- variance_tidy1 %>% 
#   filter(is.na(intensity) | intensity == i_sp)
# 
# variance_summ <- variance_tidy2 %>%
#   group_by(trait1_h2, trait2_h2, gencor, arch, selection, variable, type, trait) %>%   
#   summarize_at(vars(value), funs(mean(., na.rm = T), sd(., na.rm = T), n())) %>%
#   ungroup() %>%
#   mutate(stat = qt(p = 1 - (alpha / 2), df = n - 1) * (sd / sqrt(n) ),
#          lower = mean - stat, upper = mean + stat,
#          type = factor(type, c("tp_base", "tp_select", "c1_all", "c1_select")))
# 



## Plot the progression of the correlations
g_delta_corG <- correlations_summ %>% 
  filter(variable == "cor") %>% 
  # filter(trait1_h2 == 0.6, trait2_h2 != 1) %>% 
  ggplot(aes(x = type, y = mean, color = selection, group = selection, ymin = lower, ymax = upper)) +
  geom_point() + 
  geom_line() + 
  geom_errorbar(width = 0.2) + 
  facet_grid(trait1_h2 + trait2_h2 ~ arch + gencor) + 
  theme_acs()

# Save
ggsave(filename = "gencor_delta_corG_selection.jpg", plot = g_delta_corG, path = fig_dir, width = 10, height = 10, dpi = 1000)


## Plot the progression of the correlations - sample
g_delta_corG <- correlations_summ %>% 
  filter(variable == "cor") %>% 
  filter(trait2_h2 == 0.3) %>%
  ggplot(aes(x = type, y = mean, color = selection, group = selection, ymin = lower, ymax = upper)) +
  geom_point() + 
  geom_line() + 
  geom_errorbar(width = 0.2) + 
  facet_grid(trait1_h2 + trait2_h2 ~ arch + gencor, labeller = labeller(trait1_h2 = label_both, trait2_h2 = label_both)) + 
  theme_acs()

# Save
ggsave(filename = "gencor_delta_corG_selection_sample.jpg", plot = g_delta_corG, path = fig_dir, width = 10, height = 4, dpi = 1000)



## Plot the progression of genetic variance


















### Genetic correlation recurrent selection simulation

# Create a vector of colors to use
selection_replace <- c(mean = "Family\nmean", muspC = "Correlated\nresponse", rand = "Random")
selection_color <- set_names(umn_palette(2, 5)[3:5], selection_replace)
  
  

# Load the simulation results
load(file.path(result_dir, "popvar_gencor_recurrent_selection_simulation_results.RData"))


## Unnest
sim_selection_tidy <- popvar_gencor_selection_simulation_out %>%
  unnest(results) %>%
  mutate(sd = sqrt(var)) %>%
  select(-var) %>%
  gather(variable, value, mean, sd, cor) %>% 
  filter(!(variable == "cor" & trait == "trait2")) %>%
  mutate_at(vars(trait1_h2, trait2_h2, gencor, selection, arch, population), as.factor) %>%
  mutate(arch = factor(str_replace_all(arch, arch_replace), level = arch_replace),
         selection = factor(str_replace_all(selection, selection_replace), level = selection_replace))

## Add the base population variables for the response to selection
sim_selection_response <- sim_selection_tidy %>%
  filter(variable != "cor") %>%
  left_join(., spread(subset(sim_selection_tidy, cycle == "0" & variable != "cor", -c(cycle, population)), variable, value)) %>%
  group_by(trait1_h2, trait2_h2, gencor, selection, arch, iter, trait, variable) %>%
  mutate(response = (value - value[1]) / sd) %>%
  ungroup()

# Create an index
sim_selection_response_index <- sim_selection_response %>% 
  filter(variable == "mean") %>% 
  group_by(trait1_h2, trait2_h2, gencor, selection, arch, iter, cycle, population) %>% 
  summarize(response = mean(response)) %>%
  ungroup() %>%
  mutate(trait = "index", variable = "mean")


## Calculate the covariance between traits
sim_selection_covariance <- sim_selection_tidy %>%
  filter(variable == "cor") %>% 
  left_join(., sim_selection_tidy %>% filter(variable == "sd") %>% spread(trait, value), 
            by = c("trait1_h2", "trait2_h2", "gencor", "selection", "arch", "iter", "cycle", "population")) %>% 
  mutate(variable = "cov", response = (value * (trait1 * trait2))) %>% 
  select(trait1_h2:trait, variable, response)
  


## Combine
sim_selection_summ <- bind_rows(sim_selection_response, sim_selection_response_index, 
                                rename(filter(sim_selection_tidy, variable == "cor"), response = value),
                                sim_selection_covariance) %>%
  group_by(trait1_h2, trait2_h2, gencor, arch, selection, cycle, population, trait, variable) %>%
  summarize_at(vars(response), funs(mean(., na.rm = T), sd(., na.rm = T), n())) %>%
  ungroup() %>%
  mutate(stat = qt(p = 1 - (alpha / 2), df = n - 1) * (sd / sqrt(n) ),
         lower = mean - stat, upper = mean + stat)
  



## Fit a model to look at response at the parents at two cycles
fit_list <- sim_selection_response_index %>% 
  filter(cycle %in% c(2, 11), population == "parents") %>%
  split(.$cycle) %>%
  map(~lm(response ~ trait2_h2 + gencor + selection + arch + selection:trait2_h2 + gencor:arch + selection:arch, data = .))

map(fit_list, anova)

effs_list <- map(fit_list, ~effects::allEffects(.))
eff_df_list <- map(effs_list, ~map(., as.data.frame))

# Replot
plot_list <- pmap(list(eff_df_list, names(eff_df_list)), function(eff_df, cycle) {

  cycle1 <- as.numeric(cycle) - 1
  
  g1 <- eff_df$`trait2_h2:selection` %>%
    ggplot(aes(x = trait2_h2, y = fit, ymin = lower, ymax = upper, color = selection)) + 
    geom_point(position = position_dodge(0.75)) + 
    geom_errorbar(position = position_dodge(0.75), width = 0.5) +
    scale_color_manual(name = "Selection\nmethod", values = selection_color) +
    labs(x = "Trait 2 Heritability", y = "Standardized genotypic value (index)") +
    theme_acs() +
    labs(subtitle = paste("Cycle", cycle1)) + 
    theme(legend.position = c(0.80, 0.25))
  
  g2 <- eff_df$`selection:arch` %>%
    mutate(arch = factor(arch, levels = arch_replace)) %>%
    ggplot(aes(x = arch, y = fit, ymin = lower, ymax = upper, color = selection)) + 
    geom_point(position = position_dodge(0.75)) + 
    geom_errorbar(position = position_dodge(0.75), width = 0.5) +
    scale_color_manual(name = "Selection\nmethod", values = selection_color, guide = FALSE) +
    labs(x = "Genetic architecture", y = "Standardized genotypic value (index)") +
    theme_acs()
  
  g3 <- eff_df$`gencor:arch` %>%
    mutate(arch = factor(arch, levels = arch_replace)) %>%
    ggplot(aes(x = gencor, y = fit, ymin = lower, ymax = upper, color = arch)) + 
    geom_point(position = position_dodge(0.75)) + 
    geom_errorbar(position = position_dodge(0.75), width = 0.5) +
    scale_color_discrete(name = "Genetic\narchitecture") +
    labs(x = "Genetic correlation", y = "Standardized genotypic value (index)") +
    theme_acs() +
    theme(legend.position = c(0.25, 0.75), axis.title.y = element_blank())
  
  if (cycle1 != 1) {
    plot_grid(g1, g2, g3, nrow = 1, align = "h")
  } else {
    plot_grid(
      g1 + theme(legend.position = "none"), 
      g2, 
      g3 + theme(legend.position = "none"), 
      nrow = 1, align = "h")
    
  }
    
  
})


# Combine and Save
g_fit_combine <- plot_grid(plotlist = plot_list, ncol = 1)
ggsave(filename = "gencor_index_model.jpg", plot = g_fit_combine, path = fig_dir, width = 8, height = 6, dpi = 1000)


## Notes:
## 1. Selection on the correlated response is usually better, but not significantly bettern than selection on the family
## mean. This advantange is less when the architecture is defined by loose linkage
## 2. Seleciton using the correlated response is advantageous when the trait being indirectly selected is low in heritability.


# Fit a model for trait1 and trait2 independently
resp_fit <- sim_selection_response %>%
  filter(variable == "mean") %>%
  filter(cycle == 2, population == "parents") %>%
  group_by(trait) %>%
  do(fit = lm(response ~ trait2_h2 + gencor + selection + arch + selection:trait2_h2 + gencor:arch, data = .))

resp_fit$fit %>% map(anova)

effs_list <- resp_fit$fit %>% map(~effects::allEffects(.)) %>% setNames(., resp_fit$trait)
plot(effs_list)
eff_df <- map(effs, as.data.frame)








## Plot
# Parents selected from first cycle
sim_selection_summ %>% 
  filter(cycle == 2, population == "parents", variable == "mean") %>% 
  ggplot(aes(x = arch, y = mean, color = selection, shape = trait, ymin = lower, ymax = upper)) +
  geom_point(position = position_dodge(0.75)) +
  geom_errorbar(position = position_dodge(0.75), width = 0.5) +
  facet_grid(trait1_h2 + trait2_h2 ~ gencor) +
  theme_acs()



# Index response
g_index_response <- sim_selection_summ %>% 
  filter(trait == "index", population != "parents") %>% 
  ggplot(aes(x = cycle, y = mean, color = selection, ymin = lower, ymax = upper, fill = selection)) + 
  geom_point(size = 1) +
  geom_line() +
  geom_ribbon(alpha = 0.2) +
  facet_grid(trait1_h2 + trait2_h2 ~ arch + gencor, 
             labeller = labeller(gencor = function(x) paste("Correlation:", x),
                                 trait1_h2 = function(x) paste("Trait 1 h2:", x),
                                 trait2_h2 = function(x) paste("Trait 2 h2:", x))) +
  scale_color_manual(name = "Selection\nmethod", values = selection_color) +
  scale_fill_manual(name = "Selection\nmethod", values = selection_color) +
  scale_x_continuous(breaks = seq(0, 10, 2)) +
  ylab("Standardized genotypic mean (index)") +
  xlab("Cycle") +
  theme_acs()

ggsave(filename = "gencor_index_response.jpg", plot = g_index_response, path = fig_dir,
       height = 4, width = 8, dpi = 1000)




## Individual trait response
g_trait_response <- sim_selection_summ %>% 
  filter(trait != "index", variable == "mean", population != "parents") %>% 
  # Nudge the trait2 data upwards
  mutate(mean = ifelse(trait == "trait2", mean + 2, mean), lower = ifelse(trait == "trait2", lower + 2, lower),
         upper = ifelse(trait == "trait2", upper + 2, upper)) %>%
  ggplot(aes(x = cycle, y = mean, color = selection, ymin = lower, ymax = upper, fill = selection, shape = trait, lty = trait)) + 
  geom_point() +
  geom_line() +
  geom_ribbon(alpha = 0.25) +
  scale_color_manual(name = "Selection\nmethod", values = selection_color) +
  scale_fill_manual(name = "Selection\nmethod", values = selection_color) +
  scale_linetype_discrete(name = "Trait") +
  scale_shape_discrete(name = "Trait") +
  scale_x_continuous(breaks = seq(0, 10, 2)) +
  scale_y_continuous(sec.axis = sec_axis(~ . - 2, name = "Standardized genotypic mean (Trait 2)")) +
  labs(y = "Standardized genotypic mean (Trait 1)", x = "Cycle") +
  facet_grid(trait1_h2 + trait2_h2 ~ arch + gencor, 
             labeller = labeller(gencor = function(x) paste("Correlation:", x),
                                 trait1_h2 = function(x) paste("Trait 1 h2:", x),
                                 trait2_h2 = function(x) paste("Trait 2 h2:", x))) +
  theme_acs() +
  theme(strip.placement = "outside")

ggsave(filename = "gencor_trait_response.jpg", plot = g_trait_response, path = fig_dir,
       height = 4, width = 8, dpi = 1000)


# Value to shift the y axis
y_shift <- 1

## Standard deviations of each trait
g_trait_genvar <- sim_selection_summ %>% 
  filter(variable == "sd", population != "parents") %>% 
  # Nudge the trait2 data upwards
  mutate(mean = ifelse(trait == "trait2", mean + y_shift, mean), lower = ifelse(trait == "trait2", lower + y_shift, lower),
         upper = ifelse(trait == "trait2", upper + y_shift, upper)) %>%
  ggplot(aes(x = cycle, y = mean, color = selection, ymin = lower, ymax = upper, fill = selection, shape = trait, lty = trait)) + 
  geom_point() +
  geom_line() +
  geom_ribbon(alpha = 0.25) +
  scale_color_manual(name = "Selection\nmethod", values = selection_color) +
  scale_fill_manual(name = "Selection\nmethod", values = selection_color) +
  scale_linetype_discrete(name = "Trait") +
  scale_shape_discrete(name = "Trait") +
  scale_x_continuous(breaks = seq(0, 10, 2)) +
  scale_y_continuous(sec.axis = sec_axis(~ . - y_shift, name = "Standardized genotypic standard deviation (Trait 2)")) +
  labs(y = "Standardized genotypic standard deviation (Trait 1)", x = "Cycle") +
  facet_grid(trait1_h2 + trait2_h2 ~ arch + gencor, 
             labeller = labeller(gencor = function(x) paste("Correlation:", x),
                                 trait1_h2 = function(x) paste("Trait 1 h2:", x),
                                 trait2_h2 = function(x) paste("Trait 2 h2:", x))) +
  theme_acs() +
  theme(strip.placement = "outside")

ggsave(filename = "gencor_trait_genvar.jpg", plot = g_trait_genvar, path = fig_dir,
       height = 4, width = 8, dpi = 1000)








# Correlation
g_correlation <- sim_selection_summ %>% 
  filter(variable == "cor", population != "parents") %>% 
  # Create a grouping factor
  unite(group, arch, selection, gencor, trait1_h2, trait2_h2, sep = "_", remove = FALSE) %>%
  ggplot(aes(x = cycle, y = mean, color = selection, ymin = lower, ymax = upper, fill = selection, shape = population, group = group)) + 
  geom_point(size = 1) +
  geom_line() +
  geom_ribbon(alpha = 0.25) +
  ylab("Genetic correlation") +
  scale_x_continuous(breaks = seq(0, 10, 2)) +
  scale_color_manual(name = "Selection\nmethod", values = selection_color) +
  scale_fill_manual(name = "Selection\nmethod", values = selection_color) +
  facet_grid(trait1_h2 + trait2_h2 ~ arch + gencor) +
  theme_acs()

ggsave(filename = "gencor_correlation.jpg", plot = g_correlation, path = fig_dir,
       height = 4, width = 8, dpi = 1000)



g_correlation_bulmer <- sim_selection_summ %>% 
  filter(variable == "cor") %>%
  # filter(selection == "Correlated response") %>%
  # Create a grouping factor
  unite(group, arch, selection, gencor, trait1_h2, trait2_h2, sep = "_", remove = FALSE) %>%
  mutate(cycle = as.factor(cycle)) %>%
  ggplot(aes(x = cycle, y = mean, color = selection, ymin = lower, ymax = upper, fill = selection, shape = population, group = group)) + 
  geom_point(size = 1) +
  geom_line() +
  # geom_ribbon(alpha = 0.25) +
  facet_grid(trait1_h2 + trait2_h2 + selection ~ arch + gencor, 
             labeller = labeller(gencor = function(x) paste("Correlation:", x),
                                 trait1_h2 = function(x) paste("Trait 1 h2:", x),
                                 trait2_h2 = function(x) paste("Trait 2 h2:", x))) +
  scale_color_manual(name = "Selection\nmethod", values = selection_color) +
  scale_fill_manual(name = "Selection\nmethod", values = selection_color) +
  scale_shape_discrete(name = "Population") +
  ylab("Genetic correlation") +
  xlab("Cycle") +
  theme_acs()

ggsave(filename = "gencor_correlation_bulmer.jpg", plot = g_correlation_bulmer, path = fig_dir,
       height = 10, width = 8, dpi = 1000)




# Covariance
g_covariance <- sim_selection_summ %>% 
  filter(variable == "cov", population != "parents") %>% 
  filter(population != "tp") %>%
  # Create a grouping factor
  unite(group, arch, selection, gencor, trait1_h2, trait2_h2, sep = "_", remove = FALSE) %>%
  ggplot(aes(x = cycle, y = mean, color = selection, ymin = lower, ymax = upper, fill = selection, group = group)) + 
  geom_point() +
  geom_line() +
  geom_ribbon(alpha = 0.25) +
  scale_x_continuous(breaks = seq(0, 10, 2)) +
  scale_color_manual(name = "Selection\nmethod", values = selection_color) +
  scale_fill_manual(name = "Selection\nmethod", values = selection_color) +
  ylab("Genetic covariance") +
  xlab("Cycle") +
  facet_grid(trait1_h2 + trait2_h2 ~ arch + gencor) +
  theme_acs() +
  theme(legend.position = c(0.93, 0.15))

ggsave(filename = "gencor_covariance.jpg", plot = g_covariance, path = fig_dir,
       height = 5, width = 7, dpi = 1000)


