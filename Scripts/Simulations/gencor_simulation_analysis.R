## GenCorPrediction - genetic correlation simulation analysis
## 
## Author: Jeff Neyhart
## Last modified: March 27, 2019
## 

# Load the source script
repo_dir <- getwd()
source(file.path(repo_dir, "source.R"))

# Other libraries
library(cowplot)
library(broom)
library(modelr)
library(effects)

# Significance level
alpha <- 0.05

# Create a replacement vector
arch_replace <- c(loose_link = "Loose Linkage", close_link = "Tight Linkage", pleio = "Pleiotropy")
# Create a vector to replace the parameters in graphs
param_replace <- c("mu" = "Family Mean", "varG" = "Genetic Variance", "corG" = "Genetic Correlation")
# another vector, but in expression format
param_replace_exp <- c("mu" = "hat(mu)", "varG" = "hat(sigma)[G]^2", "corG" = "hat(r)[G]")


# Create a vector of colors to use
selection_replace <- c(mean = "FM", muspC = "CPM", corG = "Genetic\ncorrelation", rand = "Random")
selection_color <- set_names(c(umn_palette(2, 5)[3:5], "grey75"), selection_replace)










#### Prediction accuracy simulation

# Load the simulation results
files <- list.files(result_dir, pattern = "simulation_prediction", full.names = TRUE)
prediction_simulation_out <- vector("list", length(files))

for (i in seq_along(files)) {
  load(file = files[i])
  prediction_simulation_out[[i]] <- bind_rows(popvar_prediction_simulation_out)
}

popvar_prediction_simulation_out <- bind_rows(prediction_simulation_out)


## Determine missing combinations
(missing_cases <- popvar_prediction_simulation_out %>% 
    select(-input, -results) %>% 
    mutate_all(as.factor) %>% 
    anti_join(x = complete_(., names(.)), y = .))



## Tidy the results

pred_sim_tidy <- popvar_prediction_simulation_out %>%
  bind_cols(., as_data_frame(transpose(popvar_prediction_simulation_out$results))) %>%
  select(-results) %>%
  mutate_at(vars(-input, -summary, -metadata), as.factor) %>%
  mutate(arch = factor(str_replace_all(arch, arch_replace), levels = arch_replace),
         tp_meta = map(metadata, "tp_summ"),
         mar_meta = map(metadata, "mar_eff_meta"),
         pred_exp_corG = map(metadata, "pred_exp"))
         

sim_summary_tidy <- pred_sim_tidy %>% 
  unnest(summary)

sim_tp_meta_tidy <- pred_sim_tidy %>% 
  unnest(tp_meta)











## What are the genetic correlations in the tp for each tp size and intended gencor?
## Summarize the mean and sd

# Plot
g_base_corG_list <- sim_tp_meta_tidy %>% 
  filter(trait == "trait1") %>%
  group_by(tp_size, arch, gencor, nQTL) %>%
  mutate(n = n()) %>%
  ungroup() %>%
  split(.$tp_size) %>%
  map(~{
    df <- .
    
    ggplot(data = df, aes(x = cor, fill = gencor)) + 
      geom_density(alpha = 0.5) + 
      scale_fill_brewer(name = expression("Target"~r[G]), palette = "Set2", guide = guide_legend(title.position = "left", label.position = "right")) +
      facet_grid(nQTL ~ arch, labeller = labeller(nQTL = label_both), switch = "y") +
      xlim(c(-1, 1)) +
      ylab("Density") +
      xlab(expression("Realized base"~r[G])) + 
      labs(subtitle = bquote(N[TP]~"="~.(unique(parse_number(df$tp_size)))~", n ="~.(unique(df$n)))) + 
      # theme_acs() +
      theme_genetics(base_size = 10) +
      theme(legend.position = "top", legend.direction = "horizontal", strip.placement = "outside")
    
  })

## Combine
g_base_corG_combine <- plot_grid(plotlist = map(g_base_corG_list, ~. + theme(legend.position = "none")), 
                                                ncol = 2, labels = LETTERS[seq_along(g_base_corG_list)])

## Add legend
g_base_corG_combine1 <- plot_grid(g_base_corG_combine, get_legend(g_base_corG_list[[1]]), ncol = 1, rel_heights = c(1, 0.05))


ggsave(filename = "gencor_base_corG.jpg", plot = g_base_corG_combine1, path = fig_dir, height = 6, width = 8, dpi = 1000)
# ggsave(filename = "gencor_base_corG_presentation.jpg", plot = g_base_corG, path = fig_dir, height = 6.5, width = 8, dpi = 1000)



## For each nTP, compare the distributions of base correlations
base_corG_dist_test <- sim_tp_meta_tidy %>% 
  filter(trait == "trait1") %>%
  group_by(tp_size) %>%
  do({
    df <- .
    
    ## Nest the correlation by nQTL, gencor, and arch
    df_nest <- df %>%
      select(nQTL, gencor, arch, cor) %>%
      nest(cor) %>%
      mutate(data = map(data, "cor"))
      
    
    ## Create pairs of comparisons
    ## Then perform ks.test for each pair and correct for multiple testing
    comp_pars_test <- crossing(dist1 = df_nest, dist2 = df_nest) %>% 
      filter(!(gencor == gencor1 & arch == arch1 & nQTL == nQTL1)) %>%
      mutate(test = map2(.x = data, .y = data1, .f = ~ks.test(x = .x, y = .y)),
             statistic = map_dbl(test, "statistic"),
             p_value = map_dbl(test, "p.value"),
             p_value_adj = p.adjust(p = p_value, method = "bonf"))
    
    ## Return the tests and results
    comp_pars_test %>%
      select(-data, -data1, -test)
    
  })
    
## Notes
## 1. There were no identical distributions when comparing different target genetic correlations
## 






## Fit a model to compare prediction accuracy
models2 <- sim_summary_tidy %>% 
  mutate(zscore = ztrans(accuracy)) %>%
  group_by(trait, parameter) %>%
  do(fit = {
    df <- .
    
    ## Fit a full model for correlation
    fit <- lm(zscore ~ trait1_h2 + trait2_h2 + nQTL + tp_size + gencor + arch + model + nQTL:model + gencor:arch +
                nQTL:arch + model:arch + nQTL:model:arch, data = df)
    
    # Backwards regression
    step(fit, direction = "backward")
    
  }) %>% ungroup()


models2$fit %>% map(anova)


## Effect of nQTL




## Sort terms based on variance explained for correlations
var_exp1 <- models2 %>% 
  mutate(anova = map(fit, ~tidy(anova(.)) %>% mutate(prop_exp = sumsq / sum(sumsq)) %>% arrange(desc(prop_exp)))) %>% 
  unnest(anova) %>%
  select(trait, parameter, term, prop_exp) %>%
  unite(trait_param, c("trait", "parameter"), sep = "_") %>%
  mutate(per_exp = paste0(formatC(prop_exp * 100, digits = 3), "%"))




models2_all_effs <- models2 %>%
  mutate(effs = map(fit, allEffects))

# Plot effects
plot(subset(models2_all_effs, parameter == "corG", effs, drop = T)[[1]])
plot(subset(models2_all_effs, parameter == "covG", effs, drop = T)[[1]])

plot(subset(models2_all_effs, parameter == "varG" & trait == "trait1", effs, drop = T)[[1]])
plot(subset(models2_all_effs, parameter == "varG" & trait == "trait2", effs, drop = T)[[1]])

plot(subset(models2_all_effs, parameter == "mu" & trait == "trait1", effs, drop = T)[[1]])
plot(subset(models2_all_effs, parameter == "mu" & trait == "trait2", effs, drop = T)[[1]])



## Model bias
# Subset
bias_models <- sim_summary_tidy %>%
  # filter(parameter == "corG", !is.na(bias)) %>%
  filter(abs(bias) < 10) %>%
  group_by(trait, parameter) %>%
  do(model = {
    
    df <- .
    
    bias_model <- lm(bias ~ trait1_h2 + trait2_h2 + nQTL + tp_size + gencor + arch + model + nQTL:model + gencor:arch +
                       nQTL:arch + model:arch + nQTL:model:arch, data = df)
    
    step(bias_model, direction = "backward")
    
  }) %>% ungroup()

## Terms
bias_models %>% mutate(terms = map_chr(model, ~as.character(terms(.))[3]))
    

## Variance bias
bias_models %>% 
  filter(parameter == "varG") %>% 
  mutate(effects = map(model, ~Effect(focal.predictors = c("arch", "model"), .) %>% as.data.frame())) %>% 
  unnest(effects) %>% 
  mutate(arch = factor(arch, levels = arch_replace)) %>%
  ggplot(aes(x = arch, y = fit)) +
  geom_col() +
  facet_grid(model ~ trait)

## Covariance bias
plot(Effect(focal.predictors = c("model", "arch"), subset(bias_models, parameter == "covG", model, drop = T)[[1]]))
    

## Effects of model and architecture
bias_effect <- Effect(focal.predictors = c("model", "arch"), subset(bias_models, parameter == "corG", model, drop = T)[[1]]) %>% as.data.frame()

## Combine bias for correlation, covariance, and variance
all_bias <- bias_models %>% 
  filter(parameter %in% c("corG", "covG", "varG")) %>%
  mutate(effects = map(model, ~Effect(focal.predictors = c("arch", "model"), .) %>% as.data.frame())) %>% 
  unnest(effects) %>%
  mutate(arch = factor(arch, levels = arch_replace),
         parameter = str_replace_all(string = parameter, pattern = c("corG" = "Correlation", "covG" = "Covariance", "varG" = "Variance")),
         trait = str_to_title(trait), 
         group = ifelse(parameter == "Variance", paste(trait, parameter), parameter))

## function for naming
lablr <- function(x, )

g_all_bias <- all_bias %>% 
  filter(parameter != "varG") %>%
  ggplot(aes(x = model, y = fit, fill = arch, ymin = lower, ymax = upper)) +
  geom_col(position = position_dodge(0.9)) +
  geom_errorbar(position = position_dodge(0.9), width = 0.5) +
  scale_fill_brewer(palette = "Set1", name = NULL) +
  scale_y_continuous(breaks = pretty) +
  ylab("Bias") +
  xlab("Model") +
  facet_grid(~ group) +
  theme_genetics() +
  theme(legend.position = "bottom", legend.text = element_text(size = 6), legend.title = element_text(size = 8),
        axis.text.x = element_text(angle = 45, hjust = 1))
  

ggsave(filename = "model_arch_bias.jpg", plot = g_all_bias, path = fig_dir, width = 10, height = 10, units = "cm", dpi = 1000)












## Summarize over iterations
pred_sim_summary <- sim_summary_tidy %>% 
  gather(variable, value, accuracy, bias) %>% 
  ## Filter bias
  filter(!(variable == "bias" & abs(value) > 10)) %>%
  group_by(trait1_h2, trait2_h2, nQTL, tp_size, gencor, arch, model, trait, parameter, variable) %>%
  summarize_at(vars(value), funs(mean(., na.rm = T), sd(., na.rm = T), n())) %>%
  mutate(se = sd / sqrt(n), stat = se * qt(p = 1 - (alpha / 2), df = n - 1), lower = mean - stat, upper = mean + stat) %>%
  ungroup() %>%
  rename(fit = mean)



## Predictions of corG
# Color for number of QTL
color_qtl <- neyhart_palette("umn2")[3:4]

## Rename BayesC to BayesCpi
model_rename <- function(x) ifelse(x == "BayesC", expression("BayesC"*pi), x)


## Plot the effects of a subset. Break down by... Heritability?
g_pred_corG_list <- pred_sim_summary %>%
  filter(parameter == "corG", variable == "accuracy") %>%
  # filter(trait1_h2 != 1, trait2_h2 != 1, gencor == 0.5) %>%
  mutate(group = paste(nQTL, model, sep = "_"),
         herit = paste0("h[1]^2==~", trait1_h2, "~'/'~h[2]^2==~", trait2_h2)) %>%
  split(.$gencor) %>%
  map(~{
    ggplot(data = ., aes(x = tp_size, y = fit, color = nQTL, lty = model, shape = model, group = group)) +
      geom_ribbon(aes(ymin = lower, ymax = upper, fill = nQTL), alpha = 0.15, color = 0) +
      # geom_point() +
      geom_line(size = 0.5) +
      # geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.5) +
      scale_linetype_discrete(name = "Model", labels = model_rename) +
      scale_shape_discrete(name = "Model", labels = model_rename) +
      scale_color_manual(values = color_qtl) +
      scale_fill_manual(values = color_qtl) +
      scale_y_continuous(name = expression(r[G]~" prediction accuracy"), breaks = pretty, limits = c(0.15, 0.85)) +
      scale_x_discrete(name = "Training population size") +
      facet_grid(arch ~ herit, labeller = labeller(herit = label_parsed), switch = "y") +
      theme_genetics(base_size = 10) +
      theme(legend.position = "bottom") +
      labs(subtitle = bquote(r[G(0)]==.(unique(parse_number(.$gencor)))))
  })

## Combine
g_pred_corG_combine <- plot_grid(plotlist = map(g_pred_corG_list, ~. + theme(legend.position = "none")) , ncol = 1, labels = LETTERS[seq_along(g_pred_corG_list)])
g_pred_corG_combine1 <- plot_grid(g_pred_corG_combine, get_legend(g_pred_corG_list[[1]]), ncol = 1, rel_heights = c(1, 0.05))

ggsave(filename = "gencor_accuracy_full.jpg", plot = g_pred_corG_combine1, path = fig_dir, width = 10, height = 15, dpi = 1000)
  



## Subset for manuscript figure
g_pred_corG_paper <- pred_sim_summary %>%
  filter(parameter == "corG", variable == "accuracy") %>%
  filter(trait1_h2 != 1, trait2_h2 != 1, gencor == 0.5) %>%
  mutate(group = paste(nQTL, model, sep = "_"),
         herit = paste0("h[1]^2==~", trait1_h2, "~'/'~h[2]^2==~", trait2_h2)) %>%
  ggplot(data = ., aes(x = tp_size, y = fit, color = nQTL, lty = model, shape = model, group = group)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = nQTL), alpha = 0.15, color = 0) +
  geom_point(size = 0.75) +
  geom_line() +
  # geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.5) +
  scale_linetype_discrete(name = "Model", labels = model_rename) +
  scale_shape_discrete(name = "Model", labels = model_rename) +
  scale_color_manual(values = color_qtl) + 
  scale_fill_manual(values = color_qtl) +
  scale_y_continuous(name = expression(r[G]~" prediction accuracy"), breaks = pretty) +
  scale_x_discrete(name = "Training population size") +
  facet_grid(arch ~ herit, labeller = labeller(herit = label_parsed), switch = "y") +
  labs(subtitle = bquote(r[G(0)]=="0.5")) +
  theme_genetics(base_size = 8) +
  theme(legend.position = "bottom")

# Save
ggsave(filename = "gencor_accuracy_paper.jpg", plot = g_pred_corG_paper, path = fig_dir, width = 10, height = 10, units = "cm", dpi = 1000)





param_colors <- neyhart_palette("fall")[c(1,3,5)]

## Highlight the difference between correlation and mean/variance
pred_comparison_all_list <- pred_sim_summary %>%
  filter(parameter %in% c("mu", "corG", "varG"), variable == "accuracy") %>%
  mutate(group = paste(trait, parameter, sep = "_"),
         herit = paste0("h[1]^2==~", trait1_h2, "~'/'~h[2]^2==~", trait2_h2),
         parameter = factor(str_replace_all(parameter, param_replace), levels = param_replace)) %>%
  split(list(.$gencor, .$model, .$nQTL)) %>%
  map(~ggplot(data = ., aes(x = tp_size, y = fit, color = parameter, lty = trait, group = group)) +
        geom_point() +
        geom_line() +
        # geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.5) +
        scale_linetype_discrete(name = "Trait", guide = guide_legend(nrow = 2), labels = str_to_title) +
        scale_color_manual(name = "Parameter", values = param_colors, guide = guide_legend(nrow = 2)) +
        scale_y_continuous(name = "Prediction accuracy", breaks = pretty) +
        scale_x_discrete(name = "Training population size") +
        facet_grid(arch ~ herit, labeller = labeller(herit = label_parsed), switch = "y") +
        labs(subtitle = bquote(r[G(0)]==.(unique(parse_number(.$gencor)))*", Model:"~.(unique(as.character(.$model)))*", nQTL ="~.(unique(parse_number(.$nQTL))))) +
        theme_presentation2(base_size = 6) +
        theme(legend.position = "bottom") )


g_pred_comparison_all <- plot_grid(plotlist = map(pred_comparison_all_list, ~. + theme(legend.position = "none")), nrow = 4)
g_pred_comparison_all1 <- plot_grid(g_pred_comparison_all, get_legend(pred_comparison_all_list[[1]]), ncol = 1, rel_heights = c(1, 0.05))
# Save
ggsave(filename = "parameter_accuracy_compare_all.jpg", plot = g_pred_comparison_all1, width = 15, height = 20, path = fig_dir, dpi = 500)




## Subset for manuscript
g_pred_comparison_paper <- pred_sim_summary %>%
  filter(parameter %in% c("mu", "corG", "varG"), variable == "accuracy") %>%
  filter(gencor == 0.5, model == "RRBLUP", nQTL == 100) %>%
  filter_at(vars(contains("h2")), all_vars(. != 1)) %>%
  mutate(group = paste(trait, parameter, sep = "_"),
         herit = paste0("h[1]^2==~", trait1_h2, "~'/'~h[2]^2==~", trait2_h2),
         parameter = factor(str_replace_all(parameter, param_replace), levels = param_replace)) %>%
  ggplot(data = ., aes(x = tp_size, y = fit, color = parameter, lty = trait, group = group)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = parameter), alpha = 0.15, color = 0) +
  geom_point(size = 0.5) +
  geom_line() +
  scale_linetype_discrete(name = NULL, guide = guide_legend(nrow = 1), labels = str_to_title) +
  scale_color_manual(name = "Parameter", values = param_colors, guide = guide_legend(nrow = 1),
                     labels = list(bquote(hat(mu)), bquote(hat(sigma)[G]^2), bquote(hat(r)[G]))) +
  scale_fill_manual(name = "Parameter", values = param_colors, guide = guide_legend(nrow = 1),
                    labels = list(bquote(hat(mu)), bquote(hat(sigma)[G]^2), bquote(hat(r)[G]))) +
  scale_y_continuous(name = "Prediction accuracy", breaks = pretty) +
  scale_x_discrete(name = "Training population size") +
  facet_grid(arch ~ herit, labeller = labeller(herit = label_parsed), switch = "y") +
  labs(subtitle = bquote(r[G(0)]=="0.5, Model: RRBLUP, nQTL = 100")) +
  theme_genetics(base_size = 8) +
  theme(legend.position = "bottom")

# Save
ggsave(filename = "parameter_accuracy_compare_all_paper.jpg", plot = g_pred_comparison_paper, width = 10, height = 10, 
       units = "cm",path = fig_dir, dpi = 1000)



## Combine the correlation predictions and comparison of predictions into one figure
g_accuracy_and_compare <- plot_grid(g_pred_corG_paper, g_pred_comparison_paper, labels = LETTERS[1:2], ncol = 1)

ggsave(filename = "meta_accuracy_figure_paper.jpg", plot = g_accuracy_and_compare, width = 10, height = 20, 
       units = "cm",path = fig_dir, dpi = 1000)



## For each heritability, tp size, etc, calculate the difference between family mean - variance - correlation, then
## find the average
pred_sim_compare <- pred_sim_summary %>%
  filter(parameter %in% c("mu", "corG", "varG"), variable == "accuracy")

pred_sim_compare1 <- filter(pred_sim_compare, parameter == "corG") %>% 
  select(-trait) %>% crossing(., data_frame(trait = c("trait1", "trait2"))) %>%
  bind_rows(filter(pred_sim_compare, parameter != "corG")) %>%
  select(trait1_h2:parameter, trait, fit) %>%
  spread(parameter, fit) %>%
  mutate(mu_minus_varG = mu - varG, 
         mu_minus_corG = mu - corG,
         varG_minus_corG = varG - corG)

pred_sim_compare1 %>%
  select(corG:varG_minus_corG) %>%
  gather(difference, value) %>%
  group_by(difference) %>%
  summarize_at(vars(value), funs(mean, min, max))

# difference       mean     min   max
# 1 corG            0.477  0.183  0.805
# 2 mu              0.873  0.638  1.000
# 3 mu_minus_corG   0.395  0.146  0.744
# 4 mu_minus_varG   0.215  0.0335 0.416
# 5 varG            0.657  0.338  0.955
# 6 varG_minus_corG 0.180 -0.0943 0.500





## Range in prediction accuracies
pred_sim_summary %>% filter(parameter == "corG", variable == "accuracy") %>% pull(fit) %>% range()

# 0.1443822 0.8133484

# For each condition, calculate the average change in accuracy when going from the lowest TP size to highest TP
pred_sim_summary %>% filter(parameter == "corG", variable == "accuracy") %>% filter(tp_size %in% range(parse_number(tp_size))) %>% group_by(trait1_h2, trait2_h2, nQTL, gencor, arch, model) %>% arrange(trait1_h2, trait2_h2, nQTL, gencor, arch, model) %>% do(acc_change = {max(.$fit) / min(.$fit)}) %>% unnest() %>% pull(acc_change) %>% mean()

# 1.505811


## For each condition, calculate the difference between accuracy between the genetic architectures
pred_sim_summary %>% 
  filter(parameter == "corG", variable == "accuracy") %>% 
  group_by(trait1_h2, trait2_h2, nQTL, tp_size, gencor, model) %>%
  do({
    df <- .
    crossing(select(df, arch, fit), select(df, arch, fit)) %>% 
      filter(arch != arch1) %>% 
      mutate(diff = fit - fit1, per_diff = (fit - fit1) / fit1)
  }) %>%
  group_by(arch, arch1) %>%
  summarize(mean_diff = mean(diff), mean_per_diff = mean(per_diff), sd = sd(diff), n = n()) %>%
  ungroup() %>%
  mutate(se = sd / sqrt(n), stat = se * qt(p = alpha / 2, df = n - 1, lower.tail = FALSE),
         lower = mean_diff - stat, upper = mean_diff + stat) %>%
  filter(sign(mean_diff) == 1)

# arch          arch1         mean_diff mean_per_diff     sd     n      se    stat  lower  upper
# 1 Loose Linkage Pleiotropy       0.0173        0.0619 0.0615   432 0.00296 0.00582 0.0115 0.0231
# 2 Tight Linkage Loose Linkage    0.0707        0.168  0.0349   432 0.00168 0.00330 0.0674 0.0740
# 3 Tight Linkage Pleiotropy       0.0880        0.234  0.0523   432 0.00252 0.00495 0.0830 0.0929


## Do same for nQTL
pred_sim_summary %>% 
  filter(parameter == "corG", variable == "accuracy") %>% 
  group_by(trait1_h2, trait2_h2, arch, tp_size, gencor, model) %>%
  do({
    df <- .
    crossing(select(df, nQTL, fit), select(df, nQTL, fit)) %>% filter(nQTL != nQTL1) %>% mutate(diff = fit - fit1, per_diff = (fit - fit1) / fit1)
  }) %>%
  group_by(arch, nQTL, nQTL1) %>%
  summarize(mean_diff = mean(diff), mean_per_diff = mean(per_diff), sd = sd(diff), n = n()) %>%
  ungroup() %>%
  mutate(se = sd / sqrt(n), stat = se * qt(p = alpha / 2, df = n - 1, lower.tail = FALSE),
         lower = mean_diff - stat, upper = mean_diff + stat) %>%
  filter(sign(mean_diff) == 1)

# arch          nQTL  nQTL1 mean_diff mean_per_diff     sd     n      se    stat  lower  upper
# 1 Loose Linkage 100   30       0.0530        0.134  0.0351   216 0.00239 0.00471 0.0482 0.0577
# 2 Tight Linkage 100   30       0.0323        0.0721 0.0335   216 0.00228 0.00449 0.0278 0.0368
# 3 Pleiotropy    30    100      0.0358        0.0915 0.0330   216 0.00224 0.00442 0.0314 0.0403

## Do same for model
pred_sim_summary %>% 
  filter(parameter == "corG", variable == "accuracy") %>% 
  group_by(trait1_h2, trait2_h2, arch, tp_size, gencor, nQTL) %>%
  do({
    df <- .
    crossing(select(df, model, fit), select(df, model, fit)) %>% filter(model != model1) %>% mutate(diff = fit - fit1, per_diff = (fit - fit1) / fit1)
  }) %>%
  group_by(arch, nQTL, model, model1) %>%
  summarize(mean_diff = mean(diff), mean_per_diff = mean(per_diff), sd = sd(diff), n = n()) %>%
  ungroup() %>%
  mutate(se = sd / sqrt(n), stat = se * qt(p = alpha / 2, df = n - 1, lower.tail = FALSE),
         lower = mean_diff - stat, upper = mean_diff + stat) %>%
  filter(sign(mean_diff) == 1)

# arch          nQTL  model  model1 mean_diff mean_per_diff     sd     n      se    stat    lower   upper
# 1 Loose Linkage 30    BayesC RRBLUP   0.0144        0.0341  0.0358   108 0.00344 0.00682  0.00759  0.0212 
# 2 Loose Linkage 100   RRBLUP BayesC   0.00460       0.0149  0.0306   108 0.00295 0.00585 -0.00125  0.0104 
# 3 Tight Linkage 30    BayesC RRBLUP   0.0150        0.0305  0.0381   108 0.00367 0.00727  0.00776  0.0223 
# 4 Tight Linkage 100   RRBLUP BayesC   0.00178       0.00768 0.0346   108 0.00333 0.00661 -0.00482  0.00839
# 5 Pleiotropy    30    BayesC RRBLUP   0.0197        0.0460  0.0353   108 0.00339 0.00673  0.0130   0.0265 
# 6 Pleiotropy    100   BayesC RRBLUP   0.00636       0.0181  0.0308   108 0.00297 0.00588  0.000481 0.0122




pred_sim_summary %>%
  filter(variable == "accuracy", parameter %in% c("mu", "varG", "corG")) %>%
  select(-sd:-upper) %>%
  spread(parameter, fit) %>%
  group_by(trait1_h2, trait2_h2, nQTL, tp_size, gencor, arch, model) %>%
  mutate(corG = corG[1]) %>%
  ungroup() %>%
  ## Calculate difference
  mutate(mu_varG = mu - varG, mu_corG = mu - corG, varG_corG = varG - corG) %>%
  summarize_at(vars(mu_varG, mu_corG, varG_corG), mean)
  


## Create a supplemental table
supp_table_2 <- pred_sim_summary %>%
  filter(variable == "accuracy", parameter %in% c("mu", "varG", "corG")) %>%
  mutate_at(vars(fit, lower, upper), formatC, digits = 2, format = "f") %>%
  mutate(annotation = paste0(fit, " (", lower, ", ", upper, ")")) %>%
  select(trait1_h2:parameter, annotation) %>%
  mutate(trait = str_to_title(trait)) %>%
  unite(group, trait, parameter, sep = " ") %>% 
  spread(group, annotation) %>% 
  rename_all(str_to_title) %>% 
  select(Trait1_h2:Model, CorG = `Trait1 Corg`, `Trait1 Mu`, `Trait1 Varg`, `Trait2 Mu`, `Trait2 Varg`)

write_csv(x = supp_table_2, path = file.path(fig_dir, "compare_parameter_accuracy.csv"))







### Genetic correlation recurrent selection simulation


## Load the simulation results
# First find the files. If more than one file is present, combine the files
files <- list.files(result_dir, pattern = "recurrent", full.names = TRUE)

if (length(files) > 1) {
  selection_simulation_out <- vector("list", length(files))

  for (i in seq_along(files)) {
    load(file = files[i])
    selection_simulation_out[[i]] <- popvar_gencor_selection_simulation_out
  }
  
  ## Combine and save
  popvar_gencor_selection_simulation_out <- bind_rows(selection_simulation_out) %>%
    group_by(trait1_h2, trait2_h2, gencor, selection, arch) %>%
    mutate(iter = seq(n())) %>%
    ungroup()

} else {
  # Otherwise read in the data
  load(files)
  popvar_gencor_selection_simulation_out <- selection_simulation_out
  
}


## Determine missing combinations
(missing_cases <- popvar_gencor_selection_simulation_out %>% 
  select(-input, -results) %>% 
  mutate_all(as.factor) %>% 
  anti_join(x = complete_(., names(.)), y = .))



## Unnest
sim_selection_tidy <- popvar_gencor_selection_simulation_out %>%
  unnest(results) %>%
  mutate(sd = sqrt(var)) %>%
  select(-var) %>%
  gather(variable, value, mean, cor, sd) %>% 
  filter(!(variable == "cor" & trait == "trait2")) %>%
  mutate(arch = factor(str_replace_all(arch, arch_replace), level = arch_replace),
         selection = factor(str_replace_all(selection, selection_replace), level = selection_replace))

## Separate out the haplotype frequency results
sim_allele_freq_tidy <- sim_selection_tidy %>% 
  filter(trait == "trait1", variable == "mean") %>%
  select(trait1_h2:population, favorable, antagonistic1, antagonistic2, unfavorable) %>% 
  mutate(antagonistic = antagonistic1 + antagonistic2) %>% ## Combined frequency of both antagonistic haplotypes
  select(-antagonistic1, -antagonistic2) %>%
  gather(variable, response, favorable, unfavorable, contains("antagonistic")) %>%
  # Calculate the change in frequency from the start
  left_join(., filter(., cycle == 0) %>% select(trait1_h2:iter, variable, base_response = response)) %>%
  mutate(response = response - base_response) %>%
  select(-base_response)


## Separate QTL fixation proportion
sim_qtl_fixed_tidy <- sim_selection_tidy %>%
  filter(!(trait == "trait2" & arch == "Pleiotropy")) %>%
  select(trait1_h2:population, prop_fixed) %>%
  distinct() %>%
  gather(variable, response, prop_fixed)
  



## Add the base population variables for the response to selection
sim_selection_response <- sim_selection_tidy %>%
  select(trait1_h2:population, trait, variable, value) %>%
  filter(variable %in% c("mean", "sd"), !is.na(value)) %>%
  left_join(., spread(select(filter(., cycle == 0), -population, -cycle), variable, value)) %>%
  mutate(response = ifelse(variable == "mean", (value - mean) / sd, (value^2) / (sd^2))) ## For variance, calculate the proportion remaining
  

# Create an index
sim_selection_response_index <- sim_selection_response %>% 
  filter(variable == "mean") %>% 
  group_by(trait1_h2, trait2_h2, gencor, selection, arch, iter, cycle, population) %>% 
  summarize(response = mean(response, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(trait = "index", variable = "mean")


## Create a table of the covariance and correlation
sim_selection_cor_cov <- sim_selection_tidy %>%
  select(trait1_h2:population, variable, response = value) %>%
  filter(variable == "cor") %>%
  spread(variable, response) %>%
  ## Add the sd
  left_join(., filter(sim_selection_response, variable == "sd") %>% select(-variable, -mean:-response) %>% spread(trait, value)) %>%
  mutate(cov = cor * (trait1 * trait2)) %>%
  group_by(trait1_h2, trait2_h2, gencor, selection, arch, iter) %>%
  # mutate(cov = scale(cov)) %>%
  ungroup() %>%
  select(-trait1, -trait2) %>%
  gather(variable, response, cor, cov)
  


## Combine
sim_selection_summ <- bind_rows(sim_selection_response, sim_selection_response_index, sim_selection_cor_cov, sim_allele_freq_tidy, sim_qtl_fixed_tidy) %>%
  filter(!is.na(response)) %>%
  filter(gencor != 0) %>% # Don't look at gencor == 0
  group_by(trait1_h2, trait2_h2, gencor, arch, selection, cycle, population, trait, variable) %>%
  summarize_at(vars(response), funs(mean(., na.rm = T), sd(., na.rm = T), n())) %>%
  ungroup() %>%
  mutate(stat = qt(p = 1 - (alpha / 2), df = n - 1) * (sd / sqrt(n) ),
         lower = mean - stat, upper = mean + stat) %>%
  ## add annotation for heritability and genetic correlation
  mutate(herit = paste0("h[1]^2:~", trait1_h2, "~'/'~h[2]^2:~", trait2_h2),
         cor = paste0("r[G(0)]:~", gencor))
         





## Plot


# Index response - separate by correlation
# Subset data
data_toplot <- sim_selection_summ %>% 
  filter(trait == "index", population != "parents")

g_index_response <- data_toplot %>% 
  split(.$gencor) %>%
  map(~{
    ggplot(data = ., aes(x = cycle, y = mean, color = selection, ymin = lower, ymax = upper, fill = selection)) + 
      # geom_point(size = 0.5) +
      geom_line(lwd = 0.5) +
      geom_ribbon(alpha = 0.2, lwd = 0, color = 0) +
      # geom_errorbar(width = 0.5) +
      facet_grid(herit ~ arch, labeller = labeller(herit = label_parsed), switch = "y") +
      # facet_grid(arch ~ herit, labeller = labeller(herit = label_parsed), switch = "y") +
      scale_color_manual(name = "Selection\nmethod", values = selection_color) +
      scale_fill_manual(name = "Selection\nmethod", values = selection_color) +
      scale_x_continuous(breaks = seq(0, 10, 2)) +
      scale_y_continuous(breaks = pretty, limits = c(0, max(data_toplot$upper))) +
      ylab("Response (index)") +
      xlab("Cycle") +
      theme_genetics(base_size = 10) +
      theme(legend.position = "bottom") + 
      labs(subtitle = bquote(r[G(0)]==.(unique(parse_number(.$gencor)))))
  })

## Combine
g_index_combine <- plot_grid(plotlist = map(rev(g_index_response), ~. + theme(legend.position = "none")), ncol = 1, labels = LETTERS[1:3])
# g_index_combine <- plot_grid(plotlist = map(g_index_response, ~. + theme(legend.position = "none")), nrow = 1, labels = LETTERS[1:3])
g_index_combine1 <- plot_grid(g_index_combine, get_legend(g_index_response[[1]]), ncol = 1, rel_heights = c(1, 0.06))

ggsave(filename = "gencor_index_response.jpg", plot = g_index_combine1, path = fig_dir, height = 8, width = 4, dpi = 1000)


## Marginal gain of musp to mu
sim_selection_summ %>% 
  filter(trait == "index", population != "parents") %>%
  group_by(trait1_h2, trait2_h2, gencor, arch, cycle) %>% 
  do({
    df <- .
    crossing(selection = select(df, selection, mean), selection1 = selection) %>% 
      filter(selection != selection1) %>% 
      mutate(diff = mean - mean1, per_diff = (mean - mean1) / mean1)
  }) %>%
  group_by(trait1_h2, trait2_h2, gencor, arch) %>%
  ## Only look at mu and musp
  filter(selection == "CPM" & selection1 != "Random", trait2_h2 == 0.3) %>%
  filter(diff == max(diff)) %>%
  ungroup() %>% 
  arrange(trait2_h2, desc(diff))

# trait1_h2 trait2_h2 gencor arch          cycle selection  mean selection1 mean1  diff per_diff
# 1       0.6       0.3    0.5 Pleiotropy       10 CPM        4.12 FM          3.80 0.319   0.0839
# 2       0.6       0.3    0.5 Loose Linkage    10 CPM        2.93 FM          2.79 0.145   0.0522
# 3       0.6       0.3   -0.5 Tight Linkage    10 CPM        2.39 FM          2.26 0.137   0.0608
# 4       0.6       0.3    0.5 Tight Linkage    10 CPM        2.95 FM          2.82 0.133   0.0471
# 5       0.6       0.3   -0.5 Pleiotropy       10 CPM        1.98 FM          1.87 0.113   0.0603
# 6       0.6       0.3   -0.5 Loose Linkage     0 CPM        0    FM          0    0     NaN  
 
## Look at an example where the difference was zero
sim_selection_summ %>% 
  filter(trait == "index", population != "parents") %>%
  group_by(trait1_h2, trait2_h2, gencor, arch, cycle) %>% 
  do({
    df <- .
    crossing(selection = select(df, selection, mean), selection1 = selection) %>% 
      filter(selection != selection1) %>% 
      mutate(diff = mean - mean1, per_diff = (mean - mean1) / mean1)
  }) %>%
  group_by(trait1_h2, trait2_h2, gencor, arch) %>%
  ## Only look at mu and musp
  filter(selection == "CPM" & selection1 != "Random", arch == "Loose Linkage", gencor == -0.5, trait2_h2 == 0.3)




# Value to shift the y axis
y_shift_resp <- 2

data_toplot <- sim_selection_summ %>% 
  filter(trait != "index", population != "parents", variable == "mean") %>% 
  mutate_at(vars(mean, lower, upper), funs(ifelse(trait == "trait2", . + y_shift_resp, .))) 

# Trait response - separate by correlation
g_trait_response <- data_toplot %>%
  split(.$gencor) %>% 
  map(~{
    ggplot(data = ., aes(x = cycle, y = mean, color = selection, ymin = lower, ymax = upper, fill = selection, lty = trait, shape = trait)) + 
      # geom_point(size = 1) +
      geom_line(lwd = 0.5) +
      geom_ribbon(alpha = 0.2, lwd = 0, color = 0) +
      # geom_errorbar(width = 0.5) +
      facet_grid(herit ~ arch, labeller = labeller(herit = label_parsed), switch = "y") +
      # facet_grid(arch ~ herit, labeller = labeller(herit = label_parsed), switch = "y") +
      scale_color_manual(name = "Selection method", values = selection_color, guide = guide_legend(title.position = "top")) +
      scale_fill_manual(name = "Selection method", values = selection_color, guide = guide_legend(title.position = "top")) +
      scale_linetype_discrete(name = "Trait", labels = function(x) str_to_title(str_replace(x, "([a-z])([0-9])", "\\1 \\2")),
                              guide = guide_legend(title.position = "top")) +
      scale_shape_discrete(name = "Trait", labels = function(x) str_to_title(str_replace(x, "([a-z])([0-9])", "\\1 \\2")),
                              guide = guide_legend(title.position = "top")) +
      scale_x_continuous(breaks = seq(0, 10, 2), name = "Cycle") +
      scale_y_continuous(breaks = pretty, limits = c(0, max(data_toplot$upper)), name = "Response") +
      # scale_y_continuous(breaks = pretty, name = "Standardized genotypic mean (Trait 1)", limits = y_limit,
      #                    sec.axis = sec_axis(name = "Standardized genotypic mean (Trait 2)", trans = ~ . + y_shift)) +
      theme_genetics(base_size = 10) +
      theme(legend.position = "bottom", legend.text = element_text(size = 8)) +
      labs(subtitle = bquote(r[G(0)]==.(unique(parse_number(.$gencor)))))
    
  })

## Combine
g_trait_combine <- plot_grid(plotlist = map(rev(g_trait_response), ~. + theme(legend.position = "none")), ncol = 1, labels = LETTERS[1:3])
# g_trait_combine <- plot_grid(plotlist = map(g_trait_response, ~. + theme(legend.position = "none")), nrow = 1, labels = LETTERS[1:3])
g_trait_combine1 <- plot_grid(g_trait_combine, get_legend(g_trait_response[[1]]), ncol = 1, rel_heights = c(1, 0.075))

ggsave(filename = "gencor_trait_response.jpg", plot = g_trait_combine1, path = fig_dir, height = 8, width = 4, dpi = 1000)


## Look at traits individually
sim_selection_summ %>% 
  filter(trait != "index", variable == "mean", population != "parents") %>%
  group_by(trait1_h2, trait2_h2, gencor, arch, cycle, trait) %>% 
  do({
    df <- .
    crossing(selection = select(df, selection, mean), selection1 = selection) %>% 
      filter(selection != selection1) %>% 
      mutate(diff = mean - mean1, per_diff = (mean - mean1) / mean1)
  }) %>%
  group_by(trait1_h2, trait2_h2, gencor, arch, trait) %>%
  ## Only look at mu and musp
  filter(selection == "CPM" & selection1 != "Random") %>%
  filter(diff == max(diff)) %>%
  ungroup() %>% 
  arrange(trait, trait2_h2, desc(diff)) %>%
  as.data.frame()









## Combine index and trait response plots

## Alternate plot
g_index_response_alt <- sim_selection_summ %>% 
  filter(trait == "index", population != "parents") %>%
  filter(trait2_h2 == 0.3) %>%
  # filter(gencor == -0.5) %>%
  ggplot(data = ., aes(x = cycle, y = mean, color = selection, ymin = lower, ymax = upper, fill = selection)) + 
  # geom_point(size = 0.5) +
  geom_line(lwd = 0.5) +
  geom_ribbon(alpha = 0.2, lwd = 0, color = 0) +
  # geom_errorbar(width = 0.5) +
  facet_grid(cor ~ arch, labeller = labeller(herit = label_parsed, cor = label_parsed), switch = "y") +
  # facet_grid(arch ~ herit + cor, labeller = labeller(herit = label_parsed, cor = label_parsed), switch = "y") +
  scale_color_manual(name = "Selection\nmethod", values = selection_color) +
  scale_fill_manual(name = "Selection\nmethod", values = selection_color) +
  scale_x_continuous(breaks = seq(0, 10, 2)) +
  scale_y_continuous(breaks = pretty) +
  ylab("Response (index)") +
  xlab("Cycle") +
  theme_genetics(base_size = 10) +
  theme(legend.position = "bottom", legend.text = element_text(size = 8)) #, strip.text = element_text(size = 6))


## Extra plot for presentations
ggsave(filename = "gencor_index_alt_presenation.jpg", plot = g_index_response_alt, path = fig_dir, height = 10, width = 10, units = "cm", dpi = 1000)



# Alternate trait response plot
g_trait_response_alt <- sim_selection_summ %>% 
  filter(trait != "index", population != "parents", variable == "mean") %>% 
  mutate_at(vars(mean, lower, upper), funs(ifelse(trait == "trait2", . + y_shift_resp, .))) %>%
  filter(trait2_h2 == 0.3) %>%
  ggplot(data = ., aes(x = cycle, y = mean, color = selection, ymin = lower, ymax = upper, fill = selection, lty = trait)) + 
  # geom_point(size = 0.5) +
  geom_line(lwd = 0.5) +
  geom_ribbon(alpha = 0.2, lwd = 0, color = 0) +
  # geom_errorbar(width = 0.5) +
  facet_grid(cor ~ arch, labeller = labeller(herit = label_parsed, cor = label_parsed), switch = "y") +
  # facet_grid(arch ~ herit + cor, labeller = labeller(herit = label_parsed, cor = label_parsed), switch = "y") +
  scale_color_manual(name = "Selection method", values = selection_color, guide = guide_legend(title.position = "left")) +
  scale_fill_manual(name = "Selection method", values = selection_color, , guide = guide_legend(title.position = "left")) +
  scale_linetype_discrete(name = "Trait", labels = function(x) str_to_title(str_replace(x, "([a-z])([0-9])", "\\1 \\2")),
                          guide = guide_legend(title.position = "left")) +
  scale_x_continuous(breaks = seq(0, 10, 2), name = "Cycle") +
  scale_y_continuous(breaks = pretty, name = "Response (per trait)") +
  # scale_y_continuous(breaks = pretty, name = "Standardized genotypic mean (Trait 1)", limits = y_limit,
  #                    sec.axis = sec_axis(name = "Standardized genotypic mean (Trait 2)", trans = ~ . + y_shift)) +
  theme_genetics(base_size = 10) +
  theme(legend.position = "bottom", legend.text = element_text(size = 8)) #, strip.text = element_text(size = 6))


## Combine index with trait
g_index_trait_alt_combine <- plot_grid(plotlist = map(list(g_index_response_alt, g_trait_response_alt), ~. + theme(legend.position = "none")), 
                                   ncol = 2, labels = LETTERS[1:2])
g_index_trait_alt_combine1 <- plot_grid(g_index_trait_alt_combine, get_legend(g_trait_response_alt), ncol = 1, rel_heights = c(1, 0.06))

# ggsave(filename = "gencor_index_trait_alt_response_paper.jpg", plot = g_index_trait_alt_combine1, path = fig_dir, 
#        height = 20, width = 10, units = "cm", dpi = 1000)

## Save the alternate where h2_2 == 0.6 is removed
ggsave(filename = "gencor_index_trait_alt2_response_paper.jpg", plot = g_index_trait_alt_combine1, path = fig_dir, 
       height = 10, width = 20, units = "cm", dpi = 1000)













# Value to shift the y axis
y_shift <- 1

## Standard deviations of each trait
g_trait_genvar_list <- sim_selection_summ %>% 
  filter(variable == "sd", population != "parents") %>% 
  # Nudge the trait2 data upwards
  mutate_at(vars(mean, lower, upper), funs(ifelse(trait == "trait2", . + y_shift, .))) %>%
  split(.$gencor) %>%
  map(~{
    ggplot(data = ., aes(x = cycle, y = mean, color = selection, lty = trait, ymin = lower, ymax = upper, fill = selection)) + 
      # geom_point(size = 0.5) +
      geom_line(lwd = 0.5) +
      geom_ribbon(alpha = 0.2, lwd = 0, color = 0) +
      # geom_errorbar(width = 0.5) +
      facet_grid(herit ~ arch, labeller = labeller(herit = label_parsed), switch = "y") +
      # facet_grid(arch ~ herit, labeller = labeller(herit = label_parsed), switch = "y") +
      scale_color_manual(name = "Selection method", values = selection_color, guide = guide_legend(title.position = "top")) +
      scale_fill_manual(name = "Selection method", values = selection_color, guide = guide_legend(title.position = "top")) +
      scale_linetype_discrete(name = "Trait", labels = function(x) str_to_title(str_replace(x, "([a-z])([0-9])", "\\1 \\2")),
                              guide = guide_legend(title.position = "top")) +
      scale_x_continuous(breaks = seq(0, 10, 2)) +
      scale_y_continuous(breaks = pretty) +
      ylab("Standardized genetic variance") +
      xlab("Cycle") +
      theme_genetics() +
      theme(legend.position = "bottom", legend.text = element_text(size = 8)) + 
      labs(subtitle = bquote(r[G(0)]==.(unique(parse_number(.$gencor)))))
  })
  
## combine
g_trait_genvar_combine <- plot_grid(plotlist = map(g_trait_genvar_list, ~. + theme(legend.position = "none")), ncol = 1,
                                    labels = LETTERS[1:3])
g_trait_genvar_combine1 <- plot_grid(g_trait_genvar_combine, get_legend(g_trait_genvar_list[[1]]), ncol = 1, rel_heights = c(1, 0.07))

ggsave(filename = "gencor_trait_genvar.jpg", plot = g_trait_genvar_combine1, path = fig_dir, height = 8, width = 4, dpi = 1000)










### Genetic correlation across cycles
# y axis limit
y_limit <- sim_selection_summ %>%
  filter(variable == "cor", population != "parents") %>% 
  select(lower, upper) %>% unlist() %>% range()
  

g_gencor_list <- sim_selection_summ %>%
  filter(variable == "cor", population != "parents") %>% 
  split(.$gencor) %>%
  map(~{
    ggplot(data = ., aes(x = cycle, y = mean, color = selection, ymin = lower, ymax = upper, fill = selection)) + 
      # geom_point(size = 0.5) +
      geom_line(lwd = 0.5) +
      geom_ribbon(alpha = 0.2, lwd = 0, color = 0) +
      # geom_errorbar(width = 0.5) +
      facet_grid(herit ~ arch, labeller = labeller(herit = label_parsed), switch = "y") +
      # facet_grid(arch ~ herit, labeller = labeller(herit = label_parsed), switch = "y") +
      scale_color_manual(name = "Selection method", values = selection_color) +
      scale_fill_manual(name = "Selection method", values = selection_color) +
      scale_x_continuous(breaks = seq(0, 10, 2)) +
      scale_y_continuous(breaks = pretty, limits = y_limit) +
      ylab("Genetic correlation") +
      xlab("Cycle") +
      theme_genetics() +
      theme(legend.position = "bottom", legend.text = element_text(size = 8)) + 
      labs(subtitle = bquote(r[G(0)]==.(unique(parse_number(.$gencor)))))
  })


g_gencor_combine <- plot_grid(plotlist = map(g_gencor_list, ~. + theme(legend.position = "none")), ncol = 1, labels = LETTERS[1:3])
g_gencor_combine1 <- plot_grid(g_gencor_combine, get_legend(g_gencor_list[[1]]), ncol = 1, rel_heights = c(1, 0.05))

ggsave(filename = "gencor_correlation.jpg", plot = g_gencor_combine1, path = fig_dir, height = 8, width = 4, dpi = 1000)



### Genetic covariance across cycles
data_toplot <- sim_selection_summ %>%
  filter(variable == "cov", population != "parents") #%>% filter(cycle != 0)

# y axis limit
y_limit <- data_toplot %>% 
  select(lower, upper) %>% unlist() %>% range()


g_gencov_list <- data_toplot %>%
  split(.$gencor) %>%
  map(~{
    ggplot(data = ., aes(x = cycle, y = mean, color = selection, ymin = lower, ymax = upper, fill = selection)) + 
      # geom_point(size = 0.5) +
      geom_line(lwd = 0.5) +
      geom_ribbon(alpha = 0.2, lwd = 0, color = 0) +
      # geom_errorbar(width = 0.5) +
      facet_grid(herit ~ arch, labeller = labeller(herit = label_parsed), switch = "y") +
      # facet_grid(arch ~ herit, labeller = labeller(herit = label_parsed), switch = "y") +
      scale_color_manual(name = "Selection method", values = selection_color) +
      scale_fill_manual(name = "Selection method", values = selection_color) +
      scale_x_continuous(breaks = seq(0, 10, 2)) +
      scale_y_continuous(breaks = pretty, limits = y_limit) +
      ylab("Genetic covariance") +
      xlab("Cycle") +
      theme_genetics() +
      theme(legend.position = "bottom", legend.text = element_text(size = 8))  +
      labs(subtitle = bquote(r[G(0)]==.(unique(parse_number(.$gencor)))))
  })


g_gencov_combine <- plot_grid(plotlist = map(g_gencov_list, ~. + theme(legend.position = "none")), ncol = 1, labels = LETTERS[1:3])
g_gencov_combine1 <- plot_grid(g_gencov_combine, get_legend(g_gencov_list[[1]]), ncol = 1, rel_heights = c(1, 0.05))

ggsave(filename = "gencor_covariance.jpg", plot = g_gencov_combine1, path = fig_dir, height = 8, width = 4, dpi = 1000)





## Manuscript plots
y_shift <- 1

g_genvar_alt <- sim_selection_summ %>% 
  filter(population != "parents", variable == "sd") %>%
  filter(trait2_h2 == 0.3) %>%
  mutate_at(vars(mean, lower, upper), funs(ifelse(trait == "trait2", . + y_shift, .))) %>%
  ggplot(data = ., aes(x = cycle, y = mean, color = selection, lty = trait, ymin = lower, ymax = upper, fill = selection)) + 
  # geom_point(size = 0.5) +
  geom_line(lwd = 0.5) +
  geom_ribbon(alpha = 0.2, lwd = 0, color = 0) +
  # geom_errorbar(width = 0.5) +
  facet_grid(cor ~ arch, labeller = labeller(herit = label_parsed, cor = label_parsed), switch = "y") +
  # facet_grid(arch ~ herit + cor, labeller = labeller(herit = label_parsed, cor = label_parsed), switch = "y") +
  scale_linetype_discrete(name = "Trait", labels = function(x) str_to_title(str_replace(x, "([a-z])([0-9])", "\\1 \\2")),
                          guide = guide_legend(title.position = "top")) +
  scale_color_manual(name = "Selection method", values = selection_color, guide = guide_legend(title.position = "top")) +
  scale_fill_manual(name = "Selection method", values = selection_color, guide = guide_legend(title.position = "top")) +
  scale_x_continuous(breaks = seq(0, 10, 2)) +
  scale_y_continuous(breaks = pretty) +
  ylab("Standardized genetic variance") +
  xlab("Cycle") +
  theme_genetics() +
  theme(legend.position = "bottom", legend.text = element_text(size = 8) )# , strip.text = element_text(size = 6))

g_gencor_alt <- sim_selection_summ %>% 
  filter(population != "parents", variable == "cor") %>% 
  filter(trait2_h2 == 0.3) %>%
  ggplot(data = ., aes(x = cycle, y = mean, color = selection, ymin = lower, ymax = upper, fill = selection)) + 
  # geom_point(size = 0.5) +
  geom_line(lwd = 0.5) +
  geom_ribbon(alpha = 0.2, lwd = 0, color = 0) +
  # geom_errorbar(width = 0.5) +
  facet_grid(cor ~ arch, labeller = labeller(herit = label_parsed, cor = label_parsed), switch = "y") +
  # facet_grid(arch ~ herit + cor, labeller = labeller(herit = label_parsed, cor = label_parsed), switch = "y") +
  scale_color_manual(name = "Selection method", values = selection_color, guide = guide_legend(title.position = "top")) +
  scale_fill_manual(name = "Selection method", values = selection_color, guide = guide_legend(title.position = "top")) +
  scale_x_continuous(breaks = seq(0, 10, 2)) +
  scale_y_continuous(breaks = pretty) +
  ylab("Genetic correlation") +
  xlab("Cycle") +
  theme_genetics() +
  theme(legend.position = "bottom", legend.text = element_text(size = 8) )# , strip.text = element_text(size = 6))

g_gencov_alt <- sim_selection_summ %>% 
  filter(population != "parents", variable == "cov") %>% 
  filter(trait2_h2 == 0.3) %>%
  # filter(cycle != 0) %>%
  ggplot(data = ., aes(x = cycle, y = mean, color = selection, ymin = lower, ymax = upper, fill = selection)) + 
  # geom_point(size = 0.5) +
  geom_line(lwd = 0.5) +
  geom_ribbon(alpha = 0.2, lwd = 0, color = 0) +
  # geom_errorbar(width = 0.5) +
  facet_grid(cor ~ arch, labeller = labeller(herit = label_parsed, cor = label_parsed), switch = "y") +
  # facet_grid(arch ~ herit + cor, labeller = labeller(herit = label_parsed, cor = label_parsed), switch = "y") +
  scale_color_manual(name = "Selection method", values = selection_color, guide = guide_legend(title.position = "top")) +
  scale_fill_manual(name = "Selection method", values = selection_color, guide = guide_legend(title.position = "top")) +
  scale_x_continuous(breaks = seq(0, 10, 2)) +
  scale_y_continuous(breaks = pretty) +
  ylab("Standardized genetic covariance") +
  xlab("Cycle") +
  theme_genetics() +
  theme(legend.position = "bottom", legend.text = element_text(size = 8) )# , strip.text = element_text(size = 6))

## Combine
g_variability_alt <- plot_grid(plotlist = map(list(g_genvar_alt, g_gencor_alt, g_gencov_alt), ~. + theme(legend.position = "none")),
                                 ncol = 1, labels = LETTERS[1:3])
g_variability_alt1 <- plot_grid(g_variability_alt, get_legend(g_genvar_alt), ncol = 1, rel_heights = c(1, 0.05))

ggsave(filename = "gencor_cor_cov_paper.jpg", plot = g_variability_alt1, path = fig_dir, height = 25, width = 10, units = "cm", dpi = 1000)



##
g_variability_alt <- plot_grid(plotlist = map(list(g_genvar_alt, g_gencor_alt), ~. + theme(legend.position = "none")), ncol = 1, 
                               labels = LETTERS[1:2], align = "hv", axis = "tblr")

g_variability_alt1 <- plot_grid(g_variability_alt, get_legend(g_genvar_alt), ncol = 1, rel_heights = c(1, 0.07))

ggsave(filename = "gencor_cor_cov_paper2.jpg", plot = g_variability_alt1, path = fig_dir, height = 17, width = 10, units = "cm", dpi = 1000)




## Types of haplotypes
haplotype_types <- c("favorable", "antagonistic", "antagonistic1", "antagonistic2", "unfavorable")


## Change in haplotype frequencies
g_haplotype_list <- sim_selection_summ %>% 
  filter(population != "parents", variable %in% haplotype_types) %>% # Use for haplotype_frequency
  mutate(variable = factor(str_to_title(variable), levels = str_to_title(haplotype_types)))  %>%
  split(.$gencor) %>%
  map(~{
    ggplot(data = ., aes(x = cycle, y = mean, color = selection, ymin = lower, ymax = upper, fill = selection, lty = variable)) + 
      # geom_point(size = 0.5) +
      geom_line(lwd = 0.5) +
      geom_ribbon(alpha = 0.2, lwd = 0, color = 0) +
      # geom_errorbar(width = 0.5) +
      facet_grid(herit ~ arch, labeller = labeller(herit = label_parsed), switch = "y") +
      # facet_grid(arch ~ herit, labeller = labeller(herit = label_parsed), switch = "y") +
      scale_color_manual(name = "Selection method", values = selection_color) +
      scale_fill_manual(name = "Selection method", values = selection_color) +
      scale_x_continuous(breaks = seq(0, 10, 2)) +
      scale_y_continuous(breaks = pretty) +
      scale_linetype_discrete(name = "Haplotype") +
      ylab("Haplotype frequency (deviation from base population)") +
      xlab("Cycle") +
      theme_genetics() +
      theme(legend.position = "bottom", legend.text = element_text(size = 8))  +
      labs(subtitle = bquote(r[G(0)]==.(unique(parse_number(.$gencor)))))
  })

g_haplotype_combine <- plot_grid(plotlist = map(g_haplotype_list, ~. + theme(legend.position = "none")), ncol = 1, labels = LETTERS[1:3])
g_haplotype_combine1 <- plot_grid(g_haplotype_combine, get_legend(g_haplotype_list[[1]]), ncol = 1, rel_heights = c(1, 0.05))

ggsave(filename = "gencor_haplotype_freq.jpg", plot = g_haplotype_combine1, path = fig_dir, height = 8, width = 4, dpi = 1000)


## Alternative
g_haplotype_alt <- sim_selection_summ %>% 
  filter(population != "parents", variable %in% haplotype_types) %>% # Use for haplotype_frequency
  filter(trait2_h2 == 0.3) %>%
  mutate(variable = factor(str_to_title(variable), levels = str_to_title(haplotype_types)))  %>%
  ggplot(data = ., aes(x = cycle, y = mean, color = selection, ymin = lower, ymax = upper, fill = selection, lty = variable)) + 
  # geom_point(size = 0.5) +
  geom_line(lwd = 0.5) +
  geom_ribbon(alpha = 0.2, lwd = 0, color = 0) +
  # geom_errorbar(width = 0.5) +
  # facet_grid(herit ~ arch, labeller = labeller(herit = label_parsed), switch = "y") +
  facet_grid(cor ~ arch, labeller = labeller(herit = label_parsed, cor = label_parsed), switch = "y") +
  scale_color_manual(name = "Selection\nmethod", values = selection_color) +
  scale_fill_manual(name = "Selection\nmethod", values = selection_color) +
  # scale_color_manual(name = "Selection\nmethod", values = selection_color, guide = guide_legend(nrow = 2, title.position = "top")) +
  # scale_fill_manual(name = "Selection\nmethod", values = selection_color, guide = guide_legend(nrow = 2, title.position = "top")) +
  scale_x_continuous(breaks = seq(0, 10, 2)) +
  scale_y_continuous(breaks = pretty) +
  scale_linetype_discrete(name = "Haplotype") +
  ylab("Change in haplotype frequency") +
  xlab("Cycle") +
  theme_genetics() +
  theme(legend.position = "bottom", legend.direction = "horizontal", legend.box = "vertical", legend.key.width = unit(0.8, "line"),
        legend.text = element_text(size = 7), legend.box.just = "center")

# ## Create the base plot
# g_haplotype_alt_use <- g_haplotype_alt +
#   scale_linetype_discrete(name = "Haplotype\nphase", guide = FALSE)
# 
# ## Create a plot from which to extract a legned
# g_haplotype_alt_legend <- g_haplotype_alt +
#   scale_color_manual(name = "Selection\nmethod", values = selection_color, guide = FALSE) +
#   scale_fill_manual(name = "Selection\nmethod", values = selection_color, guide = FALSE) +
#   theme(legend.position = "right", legend.direction = "vertical")
# 
# g_haplotype_alt_combine <- plot_grid(g_haplotype_alt_use, get_legend(g_haplotype_alt_legend), nrow = 1, rel_widths = c(1, 0.05))
#   

ggsave(filename = "gencor_haplotype_freq_alt.jpg", plot = g_haplotype_alt, path = fig_dir, height = 15, width = 10, units = "cm", dpi = 1000)
# ggsave(filename = "gencor_haplotype_freq_alt.jpg", plot = g_haplotype_alt_combine, path = fig_dir, height = 15, width = 10, units = "cm", dpi = 1000)




## Proportion of fixed QTL
g_fixed_qtl <- sim_selection_summ %>%
  filter(variable == "prop_fixed", population != "parents") %>% # Use for haplotype_frequency
  split(.$gencor) %>%
  map(~{
    ggplot(data = ., aes(x = cycle, y = mean, color = selection, ymin = lower, ymax = upper, fill = selection, lty = trait)) + 
      # geom_point(size = 0.5) +
      geom_line(lwd = 0.5) +
      geom_ribbon(alpha = 0.2, lwd = 0, color = 0) +
      # geom_errorbar(width = 0.5) +
      facet_grid(herit ~ arch, labeller = labeller(herit = label_parsed), switch = "y") +
      # facet_grid(arch ~ herit, labeller = labeller(herit = label_parsed), switch = "y") +
      scale_color_manual(name = "Selection method", values = selection_color) +
      scale_fill_manual(name = "Selection method", values = selection_color) +
      scale_x_continuous(breaks = seq(0, 10, 2)) +
      scale_y_continuous(breaks = pretty) +
      scale_linetype_discrete(name = "Trait") +
      ylab("Haplotype frequency (deviation from base population)") +
      xlab("Cycle") +
      theme_genetics() +
      theme(legend.position = "bottom", legend.text = element_text(size = 8))  +
      labs(subtitle = bquote(r[G(0)]==.(unique(parse_number(.$gencor)))))
  })


g_fixed_qtl_combine <- plot_grid(plotlist = map(g_fixed_qtl, ~. + theme(legend.position = "none")), ncol = 1, labels = LETTERS[1:3])
g_fixed_qtl_combine1 <- plot_grid(g_fixed_qtl_combine, get_legend(g_fixed_qtl[[1]]), ncol = 1, rel_heights = c(1, 0.05))

ggsave(filename = "gencor_fixed_qtl.jpg", plot = g_fixed_qtl_combine1, path = fig_dir, height = 8, width = 4, dpi = 1000)



g_fixed_qtl_alt <- sim_selection_summ %>%
  filter(variable == "prop_fixed", population != "parents") %>% # Use for haplotype_frequency
  filter(trait2_h2 == 0.3) %>%
  ggplot(data = ., aes(x = cycle, y = mean, color = selection, ymin = lower, ymax = upper, fill = selection, lty = trait)) + 
  # geom_point(size = 0.5) +
  geom_line(lwd = 0.5) +
  geom_ribbon(alpha = 0.2, lwd = 0, color = 0) +
  # geom_errorbar(width = 0.5) +
  facet_grid(cor ~ arch, labeller = labeller(cor = label_parsed), switch = "y") +
  scale_color_manual(name = "Selection\nmethod", values = selection_color, guide = FALSE) +
  scale_fill_manual(name = "Selection\nmethod", values = selection_color, guide = FALSE) +
  scale_x_continuous(breaks = seq(0, 10, 2)) +
  scale_y_continuous(breaks = pretty) +
  scale_linetype_discrete(name = NULL, labels = str_to_title) +
  ylab("Proportion of fixed QTL") +
  xlab("Cycle") +
  theme_genetics() +
  theme(legend.position = c(0.20, 0.13), legend.text = element_text(size = 8), legend.box = "vertical", legend.box.just = "center",
        legend.margin = margin())

g_fixed_qtl_alt1 <- g_fixed_qtl_alt + 
  scale_linetype_discrete(guide = FALSE) +
  scale_color_manual(name = "Selection\nmethod", values = selection_color) +
  scale_fill_manual(name = "Selection\nmethod", values = selection_color) +
  theme(legend.position = "bottom", legend.text = element_text(size = 8), legend.box = "vertical")

g_fixed_qtl_alt_combine <- plot_grid(g_fixed_qtl_alt, get_legend(g_fixed_qtl_alt1), ncol = 1, rel_heights = c(1, 0.09))

  

ggsave(filename = "gencor_fixed_qtl_alt.jpg", plot = g_fixed_qtl_alt_combine, path = fig_dir, height = 12, width = 10, units = "cm", dpi = 1000)


### Combine haplotype and qtl fixation graphs

plot_list <- list(g_haplotype_alt + scale_fill_manual(guide = FALSE, values = selection_color) + scale_color_manual(guide = FALSE, values = selection_color), 
                  g_fixed_qtl_alt)

##
g_frequency_combine <- plot_grid(plotlist = plot_list, ncol = 1, labels = LETTERS[1:2], rel_heights = c(1, 0.85))
g_frequency_combine1 <- plot_grid(g_frequency_combine, get_legend(g_fixed_qtl_alt1), ncol = 1, rel_heights = c(1, 0.07), axis = "tblr")

ggsave(filename = "gencor_haplo_qtl_freq.jpg", plot = g_frequency_combine1, path = fig_dir, height = 20, width = 10, units = "cm", dpi = 1000)



