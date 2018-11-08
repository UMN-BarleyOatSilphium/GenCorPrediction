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
library(broom)
library(modelr)

# Significance level
alpha <- 0.05

# Create a replacement vector
arch_replace <- c(loose_link = "Loose Linkage", close_link = "Tight Linkage", pleio = "Pleiotropy")
# Create a vector to replace the parameters in graphs
param_replace <- c("mu" = "Family Mean", "varG" = "Genetic Variance", "corG" = "Genetic Correlation")




# Load the simulation results
load(file.path(result_dir, "popvar_gencor_simulation_prediction_results.RData"))
load(file.path(result_dir, "popvar_gencor_simulation_prediction_missing_results.RData"))

popvar_prediction_simulation_out <- bind_rows(popvar_prediction_simulation_out, missing_out)


## Intended number of combinations
n_expected <- popvar_prediction_simulation_out %>% summarize_at(vars(trait1_h2:arch, iter, model), n_distinct) %>% prod


## Are there any missing combinations?
n_expected - nrow(popvar_prediction_simulation_out)

(missing <- popvar_prediction_simulation_out %>%
  distinct(trait1_h2, trait2_h2, nQTL, tp_size, gencor, arch, model, iter) %>%
  mutate_all(as.factor) %>% 
  mutate(obs = T) %>% 
  complete(trait1_h2, trait2_h2, nQTL, tp_size, gencor, arch, model, iter, fill = list(obs = F)) %>% 
  filter(!obs) %>%
  mutate_at(vars(trait1_h2, trait2_h2, nQTL, tp_size, gencor, iter), parse_number) %>%
  mutate_at(vars(arch, model), parse_character))

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



## What are the genetic correlations in the tp for each tp size and intended gencor?
## Summarize the mean and sd



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


## Summarize all of the results
sim_out_summ <- sim_summary_tidy %>%
  gather(variable, value, accuracy, bias) %>% 
  group_by(trait1_h2, trait2_h2, nQTL, tp_size, gencor, arch, model, trait, parameter, variable) %>% 
  summarize_at(vars(value), funs(mean(., na.rm = T), sd(., na.rm = T), n())) %>% 
  ungroup() %>%
  mutate(stat = qt(p = 1 - (alpha / 2), df = n - 1) * (sd / sqrt(n) ),
         lower = mean - stat, upper = mean + stat,
         round_mean = round(mean, 2))




### Fit a model for all traits and parameters
models <- sim_summary_tidy %>% 
  group_by(trait, parameter) %>% 
  do(fit = lm(accuracy ~ trait1_h2 + trait2_h2 + nQTL + tp_size + gencor + arch + model, data = .)) %>%
  ungroup()

## Variance explained
models %>% 
  mutate(anova = map(fit, ~tidy(anova(.)) %>% mutate(per_exp = sumsq / sum(sumsq)) %>% arrange(desc(per_exp)))) %>% 
  unnest(anova) %>%
  select(trait, parameter, term, per_exp) %>%
  spread(term, per_exp)

## We might expect an interaction between:
## 1. nQTL and model
## 2. arch and model
## 3. The three-way interaction above
models2 <- sim_summary_tidy %>% 
  group_by(trait, parameter) %>% 
  do(fit = lm(accuracy ~ trait1_h2 + trait2_h2 + nQTL + tp_size + gencor + arch + model + nQTL:model + model:arch + nQTL:model:arch + trait1_h2:trait2_h2, data = .)) %>%
  ungroup()

## Variance explained
models2 %>% 
  mutate(anova = map(fit, ~tidy(anova(.)) %>% mutate(per_exp = sumsq / sum(sumsq)) %>% arrange(desc(per_exp)))) %>% 
  unnest(anova) %>%
  select(trait, parameter, term, per_exp) %>%
  spread(term, per_exp)


# Effect plot
models2_effs <- models2 %>%
  mutate(effs = map(fit, ~effects::allEffects(.)))

# Plot
plot(subset(models2_effs, parameter == "corG", effs, drop = T)[[1]])
plot(subset(models2_effs, parameter == "covG", effs, drop = T)[[1]])

plot(subset(models2_effs, parameter == "varG" & trait == "trait1", effs, drop = T)[[1]])
plot(subset(models2_effs, parameter == "varG" & trait == "trait2", effs, drop = T)[[1]])

# Note:
# 1. Predictions of the covariance are much more accurate under pleiotropy when using BayesC



## Highlight examples
corG_effs <- lapply(X = subset(models2_effs, parameter == "corG", effs, drop = T)[[1]], as.data.frame)

# Heritability
g_h2_summary <- corG_effs$`trait1_h2:trait2_h2` %>%
  mutate(trait2_h2 = factor(trait2_h2, levels = rev(unique(trait2_h2)))) %>%
  ggplot(aes(x = trait1_h2, y = fit, ymin = lower, ymax = upper, color = trait2_h2, group = trait2_h2)) + 
  geom_pointrange(size = 0.5) + 
  geom_line() +
  # geom_errorbar(width = 0.2) +
  scale_color_manual(name = "Trait 2\nHeritability", values = umn_palette(2)[3:5], guide = guide_legend(title.position = "left")) +
  xlab("Trait 1 Heritability") +
  ylab("Prediction accuracy") +
  scale_y_continuous(breaks = pretty) +
  theme_classic() +
  theme(legend.position = c(0.25, 0.90), legend.key.height = unit(0.75, "line"))

# ## Summarize using a fixed TP, architecture, gencor, nQTL, and model
# g_h2_summary <- sim_out_summ %>%
#   filter(variable == "accuracy", parameter == "corG", tp_size == 300) %>%
#   filter_at(vars(nQTL, tp_size, gencor, arch, model), all_vars(. == .[1])) %>%
#   mutate(trait2_h2 = factor(trait2_h2, levels = rev(unique(trait2_h2)))) %>%
#   {ggplot(data = ., aes(x = trait1_h2, y = mean, ymin = lower, ymax = upper, color = trait2_h2, group = trait2_h2, lty = nQTL)) + 
#   geom_pointrange() + 
#   geom_line() +
#   # geom_errorbar(width = 0.2) +
#   scale_color_manual(name = "Trait 2\nHeritability", values = umn_palette(2)[3:5], guide = guide_legend(title.position = "left")) +
#   xlab("Trait 1 Heritability") +
#   ylab("Prediction accuracy") +
#   scale_y_continuous(breaks = pretty) +
#   labs(caption = )
#   theme_classic() +
#   theme(legend.position = c(0.25, 0.90), legend.key.height = unit(0.75, "line"))}


# TP size
g_tp_summary <- corG_effs$tp_size %>%
  ggplot(aes(x = tp_size, y = fit, ymin = lower, ymax = upper, group = 1)) + 
  geom_pointrange(size = 0.5) + 
  geom_line() +
  # geom_errorbar(width = 0.2) +
  xlab("Training population size") +
  ylab("Prediction accuracy") +
  scale_y_continuous(breaks = pretty) +
  theme_classic() +
  theme(legend.position = c(0.25, 0.80), legend.key.height = unit(0.75, "line"))

# Model, nQTL, and architecture
g_model_summary <- corG_effs$`nQTL:arch:model` %>%
  mutate(arch = factor(arch, levels = arch_replace)) %>%
  ggplot(aes(x = arch, y = fit, ymin = lower, ymax = upper, color = model, shape = nQTL)) + 
  geom_pointrange(position = position_dodge(0.5), size = 0.5) + 
  # geom_errorbar(position = position_dodge(0.5), width = 0.5) +
  scale_color_manual(name = "Model", values = umn_palette(3)[3:4], guide = guide_legend(title.position = "left")) +
  scale_shape_discrete(guide = guide_legend(title.position = "left")) +
  xlab("Genetic architecture") +
  ylab("Prediction accuracy") +
  scale_y_continuous(breaks = pretty) +
  theme_classic() +
  theme(legend.position = c(0.50, 0.90), legend.key.height = unit(0.75, "line"), legend.spacing.y = unit(0.5, "line"),
        legend.direction = "horizonal", legend.box = "horizontal")


ylimit <- c(0.30, 0.70)
# Put everything on the same y scale
g_combine_summary <- plot_grid(
  g_h2_summary + ylim(ylimit), 
  g_tp_summary + ylim(ylimit) + theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()), 
  g_model_summary + ylim(ylimit) + theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()),
  nrow = 1, rel_widths = c(1, 0.90, 0.90))
  

ggsave(filename = "gencor_accuracy_summary1.jpg", plot = g_combine_summary, path = fig_dir,
       height = 3, width = 10, dpi = 1000)


##






## TP size breaks for the scale
tp_size_breaks <- parse_number(unique(sim_out_summ$tp_size))


## Plot using points
g_accuracy_list <- sim_out_summ %>% 
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
g_accuracy2_list <- sim_out_summ %>% 
  filter(variable == "accuracy", parameter == "corG") %>%
  mutate(trait2_h2 = factor(trait2_h2, levels = rev(levels(trait2_h2))),
         tp_size = parse_number(tp_size)) %>%
  mutate(group = paste0("Cor", str_remove_all(gencor, "\\."), "_QTL", nQTL)) %>%
  split(.$group) %>%
  map(~{
    df <- .
    ggplot(data = df, aes(x = tp_size, y = mean, ymin = lower, ymax = upper, color = arch, lty = model)) +
      geom_point(size = 1) +
      geom_line() + 
      geom_errorbar(width = 25) +
      scale_color_discrete(name = NULL) +
      scale_linetype_discrete(name = NULL) + 
      scale_x_continuous(breaks = tp_size_breaks) + 
      facet_grid(trait2_h2 ~ trait1_h2, labeller = label_both) +
      labs(caption = bquote(Base~Cor[G]*":"~.(unique(parse_number(df$gencor)))~", nQTL:"*~.(unique(parse_number(df$nQTL))))) +
      ylim(c(0.1, 0.90)) +
      theme_acs() +
      theme(legend.position = c(0.85, 0.75), legend.key.height = unit(0.5, "lines"))
  })

# Save
for (i in seq_along(g_accuracy2_list)) {
  filename <- paste0("gencor_simulation_accuracy", names(g_accuracy2_list)[i], ".jpg")
  ggsave(filename = filename, plot = g_accuracy2_list[[i]], path = fig_dir, height = 5, width = 5, dpi = 1000)
}





## Compare prediction accuracies for the family mean, variance, and correlation
# Get the models for the mean and variance, then extract the same effects as above
models2_effs_df <- models2_effs %>% 
  bind_cols(., as_data_frame(transpose(.$effs))) %>%
  filter(parameter != "musp") %>%
  mutate_at(vars(tp_size:`nQTL:arch:model`), ~map(., as.data.frame)) %>%
  mutate(parameter = str_replace_all(parameter, param_replace) %>% factor(., levels = param_replace),
         trait = str_to_title(trait))



# Heritability
g_h2_summary1 <- models2_effs_df %>% 
  unnest(`trait1_h2:trait2_h2`) %>%
  mutate(trait2_h2 = factor(trait2_h2, levels = rev(unique(trait2_h2)))) %>%
  ggplot(aes(x = trait1_h2, y = fit, ymin = lower, ymax = upper, color = trait2_h2, group = trait2_h2)) + 
  geom_pointrange(size = 0.5) + 
  geom_line() +
  # geom_errorbar(width = 0.2) +
  scale_color_manual(name = "Trait 2\nHeritability", values = umn_palette(2)[3:5]) +
  xlab("Trait 1 Heritability") +
  ylab("Prediction accuracy") +
  scale_y_continuous(breaks = pretty) +
  facet_grid(trait ~ parameter) +
  theme_presentation2() +
  theme(legend.position = c(0.85, 0.25), legend.key.height = unit(0.75, "line"))


ggsave(filename = "gencor_all_param_accuracy_heritability.jpg", plot = g_h2_summary1, path = fig_dir,
       height = 5, width = 7, dpi = 1000)
       

# TP size
g_tpsize_summary1 <- models2_effs_df %>% 
  unnest(tp_size) %>%
  ggplot(aes(x = tp_size, y = fit, ymin = lower, ymax = upper, group = 1)) + 
  geom_pointrange(size = 0.5) + 
  geom_line() +
  # geom_errorbar(width = 0.2) +
  xlab("Training population size") +
  ylab("Prediction accuracy") +
  scale_y_continuous(breaks = pretty) +
  facet_grid(trait ~ parameter) +
  theme_presentation2() +
  theme(legend.position = c(0.85, 0.25), legend.key.height = unit(0.75, "line"))

ggsave(filename = "gencor_all_param_accuracy_tpsize.jpg", plot = g_tpsize_summary1, path = fig_dir,
       height = 5, width = 7, dpi = 1000)



# Model and number of QTL
g_model_summary1 <- models2_effs_df %>% 
  unnest(`nQTL:arch:model`) %>%
  mutate(arch = str_replace_all(arch, pattern = " ", "\n"),
         arch = factor(arch, levels = str_replace_all(arch_replace, pattern = " ", "\n"))) %>%
  ggplot(aes(x = arch, y = fit, ymin = lower, ymax = upper, color = model, shape = nQTL)) + 
  geom_pointrange(position = position_dodge(0.5), size = 0.5) + 
  # geom_errorbar(position = position_dodge(0.5), width = 0.5) +
  scale_color_manual(name = "Model", values = umn_palette(3)[3:4], guide = guide_legend(title.position = "left")) +
  scale_shape_discrete(guide = guide_legend(title.position = "left")) +
  xlab("Genetic architecture") +
  ylab("Prediction accuracy") +
  scale_y_continuous(breaks = pretty) +
  facet_grid(trait ~ parameter) +
  theme_presentation2() +
  theme(legend.position = c(0.85, 0.25), legend.key.height = unit(0.75, "line"), legend.spacing.y = unit(0.5, "line"),
        axis.text.x = element_text(angle = 35, hjust = 1))

ggsave(filename = "gencor_all_param_accuracy_modelarch.jpg", plot = g_model_summary1, path = fig_dir,
       height = 5.5, width = 8, dpi = 1000)



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
# load(file.path(result_dir, "popvar_gencor_space_simulation_results_original.RData"))


# Mutate the architecture space combinations
sim_out1 <- popvar_corG_space_simulation_out %>% 
  mutate(probcor = map(probcor, ~`names<-`(as.data.frame(.), c("dLinkage", "pLinkage")) %>% tail(., 1))) %>% # The tail is used to remove the probabilities of pleiotropy))
  unnest(probcor) %>%
  mutate(dLinkage = ifelse(pLinkage == 0, 0, dLinkage),
         dLinkageFactor = ifelse(dLinkage == 0, "[0, 0]", paste0("(", round(dLinkage) - 5, ", ", dLinkage, "]"))) %>%
  bind_cols(., as_data_frame(transpose(.$results))) %>%
  select(-results)
  
## Are there any missing combinations?
sim_out1 %>%
  # Remove pleiotrpic 
  filter(pLinkage != 0) %>%
  distinct(trait1_h2, trait2_h2, gencor, dLinkage, pLinkage, iter) %>%
  mutate_all(as.factor) %>% 
  mutate(obs = T) %>% 
  complete(trait1_h2, trait2_h2, gencor, dLinkage, pLinkage, iter, fill = list(obs = F)) %>% 
  filter(!obs)

## Good!

# Tidy the results and extract the probability of linkage and the degree of linkage
sim_results_tidy <- sim_out1 %>%
  mutate_at(vars(trait1_h2:gencor, dLinkage, pLinkage), as.factor) %>%
  mutate(dLinkageFactor = factor(dLinkageFactor, levels = c("[0, 0]", head(unique(dLinkageFactor), -1))))
  

## Extract each dataset
correlation_tidy <- unnest(sim_results_tidy, other)
predictions_tidy <- unnest(sim_results_tidy, summary)


## Extract the training population genetic correlation
base_cor_summ <- correlation_tidy %>%
  group_by(trait1_h2, trait2_h2, nQTL, tp_size, gencor, dLinkage, dLinkageFactor, pLinkage, variable) %>%
  summarize_at(vars(value), funs(mean(., na.rm = T), sd(., na.rm = T), n())) %>%
  ungroup()

base_cor_summ1 <- bind_rows(filter(base_cor_summ, pLinkage != 0),
                            base_cor_summ %>% filter(pLinkage == 0) %>% select(-dLinkageFactor) %>% 
                              left_join(., distinct(ungroup(base_cor_summ), gencor, dLinkageFactor)) )




## Plot
g_base_cor <- base_cor_summ1 %>%
  filter(variable == "tp_gencor") %>%
  ggplot(aes(x = pLinkage, y = dLinkageFactor, fill = mean)) +
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
  group_by(trait1_h2, trait2_h2, nQTL, tp_size, gencor, dLinkage, dLinkageFactor, pLinkage, parameter, variable) %>%
  summarize_at(vars(value), funs(mean(., na.rm = T), sd(., na.rm = T), n())) %>%
  mutate(se = sd / sqrt(n), stat = se * qt(p = 1 - (alpha / 2), df = n - 1), lower = mean - stat, upper = mean + stat) %>%
  ungroup()


## Plot results for genetic correlation
g_pred_acc_corG <- pred_results_summ %>%
  filter(parameter == "corG", variable == "accuracy") %>%
  ggplot(aes(x = pLinkage, y = dLinkageFactor, fill = mean)) +
  geom_tile() +
  scale_fill_gradient(limits = c(0.40, 0.72), low = "white", high = "green", name = "Prediction\naccuracy") +
  facet_grid(~ gencor) +
  ylab("Maximum distance between QTL (cM)") +
  xlab("Proportion of linked QTL pairs") +
  theme_acs()


# bias
g_pred_bias_corG <- pred_results_summ1 %>%
  filter(parameter == "corG", variable == "bias") %>%
  ggplot(aes(x = pLinkage, y = dLinkageFactor, fill = mean)) +
  geom_tile() +
  scale_fill_gradient2(low = "red", high = "blue", name = "Bias") +
  facet_grid(~ gencor) +
  ylab("Maximum distance between QTL (cM)") +
  xlab("Proportion of linked QTL pairs") +
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




# For each gencor and dLinkage, fit a model regressing accuracy on pLinkage
fit_list <- predictions_tidy_tomodel %>% 
  group_by(gencor, dLinkage) %>% 
  do(fit = lm(accuracy ~ pLinkage, data = .))



# Treat pLinkage as numeric
fit <- lm(accuracy ~ gencor + dLinkage + pLinkage + gencor:pLinkage + dLinkage:pLinkage, data = predictions_tidy_tomodel,
          subset = dLinkage != 50)
anova(fit)
plot(effects::allEffects(fit))

## Note
## 1. Upward trend in prediction accuracy with decreasing pleiotropy, as expected.
## UPDATE: Negative trend in prediction accuracy with decreasing pleiotropy - not expected.


predictions_tidy %>%
  mutate(pLinkage = parse_number(pLinkage)) %>%
  ggplot(aes(x = pLinkage, y = accuracy, color = gencor)) + 
  geom_smooth(method = "lm") + 
  facet_wrap(~ dLinkage, ncol = 5) +
  theme_acs()

ggsave(filename = "gencor_plinkage_accuracy.jpg", plot = g_pLinkage_accuracy, path = fig_dir, width = 3, height = 3, dpi = 1000)


# Plot the relationship between degree of linkage and accuracy
# (assuming that 100% of the architecture is due to that degree of linkage)
g_pred_linkage1 <- pred_results_summ %>% 
  filter(pLinkage %in% c(0, 1), parameter == "corG", variable == "accuracy") %>%
  ggplot(aes(x = dLinkageFactor, y = mean, color = gencor, group = gencor)) + 
  geom_point() + 
  geom_smooth(method = "lm", se = FALSE) + 
  ylab("Accuracy") +
  xlab("Interval between QTL for each trait (cM)") + 
  scale_color_discrete(name = "Genetic\ncorrelation") +
  theme_acs() +
  theme(legend.position = c(0.80, 0.80))

ggsave(filename = "prediction_accuracy_space_linkage1.jpg", plot = g_pred_linkage1, path = fig_dir, width = 4.5, height = 4, dpi = 1000)

#
















### Genetic correlation recurrent selection simulation

# Create a vector of colors to use
selection_replace <- c(mean = "Family\nmean", muspC = "Correlated\nresponse", rand = "Random")
selection_color <- set_names(umn_palette(2, 5)[3:5], selection_replace)
  
  

# Load the simulation results
files <- list.files(result_dir, pattern = "recurrent", full.names = TRUE)
selection_simulation_out <- vector("list", length(files))

for (i in seq_along(files)) {
  load(file = files[i])
  selection_simulation_out[[i]] <- popvar_gencor_selection_simulation_out %>%
    mutate(iter = iter + ((i - 1) * max(iter)))
}

popvar_gencor_selection_simulation_out <- bind_rows(selection_simulation_out)


## Unnest
sim_selection_tidy <- popvar_gencor_selection_simulation_out %>%
  unnest(results) %>%
  mutate(sd = sqrt(var)) %>%
  select(-var) %>%
  gather(variable, value, mean, cor, sd) %>% 
  filter(!(variable == "cor" & trait == "trait2")) %>%
  mutate_at(vars(trait1_h2, trait2_h2, gencor, selection, arch, population), as.factor) %>%
  mutate(arch = factor(str_replace_all(arch, arch_replace), level = arch_replace),
         selection = factor(str_replace_all(selection, selection_replace), level = selection_replace))

## Separate out the haplotype frequency results
sim_allele_freq_tidy <- sim_selection_tidy %>% 
  select(trait1_h2:population, contains("freq")) %>% 
  distinct() %>% 
  mutate(unfav_freq1 = unfav_freq,
         unfav_freq2 = 1 - (pos_freq + neg_freq + unfav_freq),
         unfav_freq = unfav_freq1 + unfav_freq2) %>% 
  gather(variable, response, pos_freq, neg_freq, unfav_freq1, unfav_freq2, unfav_freq)


## Separate out the LD of haplotypes
sim_allele_LD_tidy <- sim_selection_tidy %>%
  select(trait1_h2:population, hap_LD) %>% 
  distinct() %>%
  gather(variable, response, hap_LD)


## Add the base population variables for the response to selection
sim_selection_response <- sim_selection_tidy %>%
  select(-contains("freq"), -hap_LD) %>%
  filter(variable != "cor") %>%
  left_join(., spread(subset(sim_selection_tidy, cycle == "0" & variable != "cor", -c(cycle, population)), variable, value)) %>%
  mutate(response = ifelse(variable == "mean", (value - mean) / sd, 1 + ((value - sd) / sd )))

# Create an index
sim_selection_response_index <- sim_selection_response %>% 
  filter(variable == "mean") %>% 
  group_by(trait1_h2, trait2_h2, gencor, selection, arch, iter, cycle, population) %>% 
  summarize(response = mean(response, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(trait = "index", variable = "mean")


## Calculate the covariance between traits
sim_selection_covariance <- sim_selection_tidy %>%
  select(-contains("freq"), -hap_LD) %>%
  filter(variable == "cor") %>% 
  left_join(., sim_selection_tidy %>% filter(variable == "sd") %>% spread(trait, value), 
            by = c("trait1_h2", "trait2_h2", "gencor", "selection", "arch", "iter", "cycle", "population")) %>% 
  mutate(variable = "cov", response = (value * (trait1 * trait2))) %>% 
  select(trait1_h2:trait, variable, response)
  


## Combine
sim_selection_summ <- bind_rows(sim_selection_response, sim_selection_response_index, 
                                rename(filter(sim_selection_tidy, variable == "cor"), response = value),
                                sim_selection_covariance, sim_allele_freq_tidy, sim_allele_LD_tidy) %>%
  filter(!is.na(response)) %>%
  group_by(trait1_h2, trait2_h2, gencor, arch, selection, cycle, population, trait, variable) %>%
  summarize_at(vars(response), funs(mean(., na.rm = T), sd(., na.rm = T), n())) %>%
  ungroup() %>%
  mutate(stat = qt(p = 1 - (alpha / 2), df = n - 1) * (sd / sqrt(n) ),
         lower = mean - stat, upper = mean + stat)
  





## Plot


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




## Highlight 3 cycles
g_index_list <- sim_selection_summ %>% 
  filter(trait == "index", population == "parents") %>% 
  filter(cycle %in% c(2, 6, max(cycle)), trait2_h2 == 0.3) %>% 
  mutate(cycle = cycle - 1, cycle = as.factor(cycle)) %>%
  ggplot(aes(x = arch, y = mean, ymin = lower, ymax = upper, color = selection, shape = cycle)) + 
  geom_point(position = position_dodge(0.9)) + 
  # geom_line() +
  geom_errorbar(position = position_dodge(0.9), width = 0.5) +
  scale_shape_discrete(name = "Cycle", guide = guide_legend(title.position = "left")) +
  scale_color_manual(name = "Selection\nmethod", values = selection_color, guide = guide_legend(title.position = "left")) +
  facet_grid(~ gencor, labeller = labeller(gencor = function(x) str_c("Genetic Correlation: ", x))) +
  ylab("Standardized genotypic mean (index)") +
  xlab("Genetic architecture") +
  theme_acs() +
  theme(legend.position = c(0.25, 0.80), legend.spacing.y = unit(0.2, "line"), legend.box = "horizontal",
        legend.background = element_rect(fill = alpha("white", 0)))

ggsave(filename = "gencor_index_two_cycles.jpg", plot = g_index_list, path = fig_dir, width = 5, height = 3, dpi = 1000)




# Value to shift the y axis
y_shift <- 4
## Individual trait response
g_trait_response <- sim_selection_summ %>% 
  filter(trait != "index", variable == "mean", population != "parents") %>% 
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
  scale_x_continuous(breaks = pretty) +
  # scale_y_continuous(sec.axis = sec_axis(~ . - y_shift), name = "Standardized genotypic mean (Trait 2)")) +
  labs(y = "Standardized genotypic mean", x = "Cycle") +
  facet_grid(trait1_h2 + trait2_h2 ~ arch + gencor, 
             labeller = labeller(gencor = function(x) paste("Correlation:", x),
                                 trait1_h2 = function(x) paste("Trait 1 h2:", x),
                                 trait2_h2 = function(x) paste("Trait 2 h2:", x))) +
  theme_acs() +
  theme(strip.placement = "outside")

ggsave(filename = "gencor_trait_response.jpg", plot = g_trait_response, path = fig_dir,
       height = 4, width = 10, dpi = 1000)


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

# ggsave(filename = "gencor_covariance_tp.jpg", plot = g_covariance, path = fig_dir,
#        height = 5, width = 7, dpi = 1000)
# 



haplo_replace <- c("pos_freq" = "Favorable\ncoupling", "neg_freq" = "Unfavorable\ncoupling", "unfav_freq" = "Repulsion")

## Change in haplotype frequencies
g_haplotype <- sim_selection_summ %>% 
  filter(str_detect(variable, "freq"), !str_detect(variable, "freq[1-2]")) %>% # Use for haplotype_frequency
  # filter(str_detect(variable, "freq"), variable != "unfav_freq") %>% # Use for haplotype_frequency_separate
  filter(!(str_detect(variable, "unfav_freq") & arch == "Pleiotropy"), population != "parents") %>%
  mutate(variable = str_replace_all(variable, haplo_replace)) %>%
  ggplot(aes(x = cycle, y = mean, ymin = lower, ymax = upper, color = selection, fill = selection, lty = variable)) +
  geom_point(size = 0.5) +
  geom_line(lwd = 0.5) +
  geom_ribbon(alpha = 0.25) +
  scale_x_continuous(breaks = seq(0, 10, 2)) +
  scale_color_manual(name = "Selection\nmethod", values = selection_color) +
  scale_fill_manual(name = "Selection\nmethod", values = selection_color) +
  scale_linetype_discrete(name = "Haplotype\nphase") +
  ylab("Haplotype frequency") +
  xlab("Cycle") +
  facet_grid(trait1_h2 + trait2_h2 ~ arch + gencor) +
  theme_acs()

# ggsave(filename = "haplotype_frequency_separate.jpg", plot = g_haplotype, path = fig_dir, width = 8, height = 4, dpi = 1000)
ggsave(filename = "haplotype_frequency.jpg", plot = g_haplotype, path = fig_dir, width = 8, height = 4, dpi = 1000)




## Plot LD
g_LD <- sim_selection_summ %>% 
  filter(variable == "hap_LD", population != "parents", arch != "Pleiotropy") %>%
  ggplot(aes(x = cycle, y = mean, ymin = lower, ymax = upper, color = selection, fill = selection)) +
  geom_point(size = 0.5) +
  geom_line(lwd = 0.5) +
  geom_ribbon(alpha = 0.25) +
  # scale_x_continuous(breaks = seq(0, 10, 2)) +
  scale_color_manual(name = "Selection\nmethod", values = selection_color) +
  scale_fill_manual(name = "Selection\nmethod", values = selection_color) +
  ylab("Mean haplotype linkage disequilibrium") +
  xlab("Cycle") +
  facet_grid(trait1_h2 + trait2_h2 ~ arch + gencor) +
  theme_acs()

ggsave(filename = "haplotype_LD.jpg", plot = g_LD, path = fig_dir, width = 6, height = 4, dpi = 1000)



















