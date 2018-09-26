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
fit <- lm(accuracy ~ trait1_h2 + trait2_h2 + nQTL + tp_size + gencor + arch + model + gencor:arch + trait1_h2:trait2_h2 + model:arch + model:nQTL, 
          data = sim_summary_tidy, subset = parameter == "corG")
anova(fit)

# Effect plot
plot(effects::allEffects(fit))

## Notes
## 1. Heritability is significant - duh. There is a synergistic effect between the trait heritabilities. You need
## both to be high to get high accuracy.
## 2. 100 QTL is more accuracy than 30 QTL - look back at Zhong2007 for clues here
## 3. TP size is significant - duh
## 4. Accuracy is much lower under pleiotropic genetic architecture. This makes sense because PopVar will
## assign non-zero marker effects to most markers, which will assumed to be segregating. The true QTL
## will not segregate, of course. 
## 5. No difference between RR-BLUP and BayesC on average, but BayesC does improve the predictions
## under pleiotropic genetic architecture

## Model just for pleiotropy
fit <- lm(accuracy ~ trait1_h2 + trait2_h2 + nQTL + tp_size + gencor + model + trait1_h2:trait2_h2 + model:nQTL, 
          data = sim_summary_tidy, subset = parameter == "corG" & arch == "Pleiotropy")
anova(fit)

# Effect plot
plot(effects::allEffects(fit))


## Notes:
## 1. On-average, BayesC does better than RRBLUP when the number of QTL is small and pleiotropy is the architecture


## TP size breaks for the scale
tp_size_breaks <- parse_number(unique(sim_out_summ$tp_size))


## Plot using points
g_accuracy <- sim_out_summ %>% 
  filter(parameter == "corG", variable == "accuracy") %>%
  mutate(trait2_h2 = factor(trait2_h2, levels = rev(levels(trait2_h2)))) %>%
  ggplot(aes(x = tp_size, y = mean, ymin = lower, ymax = upper, color = arch)) + 
  geom_point(size = 1) +
  geom_line() + 
  geom_errorbar(width = 0.25) +
  scale_color_discrete(name = NULL) +
  facet_grid(trait2_h2 + trait1_h2 ~ nQTL + gencor + model, labeller = label_both) +
  theme_acs() +
  theme(legend.position = c(0.85, 0.75), legend.key.height = unit(0.5, "lines"))



## Plot using points
g_accuracy1 <- sim_out_summ %>% 
  filter(parameter == "corG", variable == "accuracy", gencor == 0.5, nQTL == 100) %>%
  mutate(trait2_h2 = factor(trait2_h2, levels = rev(levels(trait2_h2))),
         tp_size = parse_number(tp_size)) %>%
  ggplot(aes(x = tp_size, y = mean, ymin = lower, ymax = upper, color = arch, lty = model)) + 
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
ggsave(filename = "gencor_simulation_accuracy.jpg", plot = g_accuracy1, path = fig_dir,
       height = 5, width = 5, dpi = 1000)


## Compare RRBLUP with BayesC at two TP sizes
sim_out_summ %>% 
  filter(parameter == "corG", variable == "accuracy", gencor == 0.5, trait1_h2 == 1, trait2_h2 == 1, tp_size %in% c(150, 600)) %>%
  ggplot(aes(x = tp_size, y = mean, ymin = lower, ymax = upper, fill = model)) +
  geom_col(position = "dodge") +
  geom_errorbar(position = position_dodge(0.9), width = 0.5) +
  facet_grid(nQTL ~ arch) +
  theme_acs()




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
  theme_acs()


# bias
g_pred_bias_corG <- pred_results_summ %>%
  filter(parameter == "corG", variable == "bias") %>%
  ggplot(aes(x = pLinkage, y = dLinkage, fill = mean)) +
  geom_tile() +
  scale_fill_gradient2(low = "red", high = "blue", name = "Bias") +
  facet_grid(~ gencor) +
  theme_acs()

# Combine
g_pred_corG <- plot_grid(g_pred_acc_corG, g_pred_bias_corG, ncol = 1, align = "hv")
ggsave(filename = "gencor_arch_space_pred.jpg", plot = g_pred_corG, path = fig_dir, width = 6, height = 5, dpi = 1000)




## Fit a model for correlation
fit <- lm(value ~ gencor + dLinkage + pLinkage + gencor:pLinkage, data = correlation_tidy, subset = variable == "tp_gencor" & dLinkage != 0)
anova(fit)
plot(effects::allEffects(fit))

## Notes
## 1.


## Fit a model for accuracy
fit <- lm(accuracy ~ gencor + dLinkage + pLinkage + gencor:pLinkage, data = predictions_tidy, subset = parameter == "corG" & dLinkage != 0)
anova(fit)
plot(effects::allEffects(fit))

# Treat pLinkage as numeric
fit <- lm(accuracy ~ gencor + dLinkage + pLinkage + gencor:pLinkage, data = mutate(predictions_tidy, pLinkage = parse_number(pLinkage)), 
          subset = parameter == "corG" & dLinkage != 0)
anova(fit)
plot(effects::allEffects(fit))

## Note
## 1. Upward trend in prediction accuracy with decreasing pleiotropy, as expected.


## Summarize this trend to visualize differently
predictions_linakge_summ <- predictions_tidy %>% 
  filter(dLinkage != 0, parameter == "corG") %>%
  mutate(dLinkage = parse_number(dLinkage)) %>%
  group_by(gencor, pLinkage, dLinkage) %>%
  summarize_at(vars(accuracy), funs(mean, sd))

## Plot lines for each dLinkage
g_pred_dlink <- predictions_linakge_summ %>% 
  ggplot(aes(x = pLinkage, y = mean, group = dLinkage, color = dLinkage)) + 
  geom_point(size = 1) + 
  geom_smooth(method = "lm", se = FALSE, lwd = 0.75) +
  facet_grid(~ gencor) +
  theme_acs()

ggsave(filename = "gencor_accuracy_dL_pL.jpg", plot = g_pred_dlink, path = fig_dir, width = 6, height = 3, dpi = 1000)

## Plot just dLinkage == 5
g_pred_dLink5 <- predictions_linakge_summ %>% 
  filter(dLinkage == 5) %>%
  ggplot(aes(x = pLinkage, y = mean, group = 1)) + 
  geom_point(size = 1) + 
  geom_smooth(method = "lm", se = T, lwd = 0.75) +
  facet_grid(~ gencor) +
  ylab("Prediction accuracy") +
  labs(subtitle = "dLinkage = 5 cM") +
  theme_acs()

ggsave(filename = "gencor_accuracy_dL_pL_example.jpg", plot = g_pred_dLink5, path = fig_dir, width = 6, height = 3, dpi = 1000)





















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













