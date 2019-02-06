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
library(effects)

# Significance level
alpha <- 0.05

# Create a replacement vector
arch_replace <- c(loose_link = "Loose Linkage", close_link = "Tight Linkage", pleio = "Pleiotropy")
# Create a vector to replace the parameters in graphs
param_replace <- c("mu" = "Family Mean", "varG" = "Genetic Variance", "corG" = "Genetic Correlation")

# Create a vector of colors to use
selection_replace <- c(mean = "Family\nmean", muspC = "Correlated\nresponse", corG = "Genetic\ncorrelation", rand = "Random")
selection_color <- set_names(c(umn_palette(2, 5)[3:5], "grey75"), selection_replace)










#### Prediction accuracy simulation

# Load the simulation results
files <- list.files(result_dir, pattern = "simulation_prediction", full.names = TRUE)
prediction_simulation_out <- vector("list", length(files))

for (i in seq_along(files)) {
  load(file = files[i])
  prediction_simulation_out[[i]] <- popvar_prediction_simulation_out
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




## For each condition, plot the distribution of pi
sim_mar_meta_tidy <- bind_rows(
  pred_sim_tidy %>% filter(model == "BayesC") %>% unnest(mar_meta),
  pred_sim_tidy %>% filter(model == "RRBLUP") %>% mutate(param = "pi", trait1 = 1, trait2 = 1)
) %>% select(-summary:-pred_exp_corG, -mar_meta) %>%
  gather(trait, pi, trait1, trait2)


sim_mar_meta_tidy %>%
  filter(tp_size == 600, gencor == 0.5) %>%
  ggplot(aes(x = arch, y = pi, fill = model)) + 
  geom_boxplot(position = position_dodge(0.9), alpha = 0.5) +
  facet_grid(trait1_h2 + trait2_h2 ~ nQTL + trait, labeller = labeller(nQTL = label_both), switch = "y")



## Get the expected and predicted corG for each simulation
sim_corG_tidy <- pred_sim_tidy %>% 
  select(trait1_h2:iter, pred_exp_corG) %>%
  unnest(pred_exp_corG)

## Plot a density of the expected corG for different conditions
data_toplot <- sim_corG_tidy %>%
  filter(!is.na(expectation)) %>%
  filter(tp_size == 600, gencor == 0.5)

ggplot(data_toplot, aes(x = expectation)) +
  geom_density(data = data_toplot, aes(x = prediction, fill = model), alpha = 0.3) +
  geom_density(fill = "grey85", alpha = 0.3) +
  facet_grid(trait1_h2 + trait2_h2 ~ nQTL + arch)







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

## Sort terms based on variance explained for correlations
var_exp1 <- models2 %>% 
  mutate(anova = map(fit, ~tidy(anova(.)) %>% mutate(prop_exp = sumsq / sum(sumsq)) %>% arrange(desc(prop_exp)))) %>% 
  unnest(anova) %>%
  select(trait, parameter, term, prop_exp) %>%
  unite(trait_param, c("trait", "parameter"), sep = "_") %>%
  mutate(per_exp = paste0(formatC(prop_exp * 100, digits = 3), "%"))


## What is the effect of correlation?




models2_all_effs <- models2 %>%
  mutate(effs = map(fit, allEffects))

# Plot effects
plot(subset(models2_effs, parameter == "corG", effs, drop = T)[[1]])
# plot(subset(models2_effs, parameter == "covG", effs, drop = T)[[1]])

plot(subset(models2_effs, parameter == "varG" & trait == "trait1", effs, drop = T)[[1]])
plot(subset(models2_effs, parameter == "varG" & trait == "trait2", effs, drop = T)[[1]])

plot(subset(models2_effs, parameter == "mu" & trait == "trait1", effs, drop = T)[[1]])
plot(subset(models2_effs, parameter == "mu" & trait == "trait2", effs, drop = T)[[1]])



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
    

## Variance bias
bias_models %>% 
  filter(parameter == "varG") %>% 
  mutate(effects = map(model, ~Effect(focal.predictors = c("arch", "model"), .) %>% as.data.frame())) %>% 
  unnest(effects) %>% 
  mutate(arch = factor(arch, levels = arch_replace)) %>%
  ggplot(aes(x = arch, y = fit)) +
  geom_col() +
  facet_grid(model ~ trait)
    


corG_bias_model <- lm(bias ~ trait1_h2 + trait2_h2 + nQTL + tp_size + gencor + arch + model + nQTL:model + gencor:arch +
                        nQTL:arch + model:arch + nQTL:model:arch, data = bias_tidy)

corG_bias_model1 <- step(corG_bias_model, direction = "backward")

## Effects of model and architecture
bias_effect <- Effect(focal.predictors = c("model", "arch"), subset(bias_models, parameter == "corG", model, drop = T)[[1]]) %>% as.data.frame()


### Plot
g_bias_example <- bias_effect %>%
  mutate(arch = factor(arch, levels = arch_replace)) %>%
  ggplot(aes(x = model, y = fit, fill = arch, ymin = lower, ymax = upper)) +
  geom_col(position = position_dodge(0.9)) +
  geom_errorbar(position = position_dodge(0.9), width = 0.5) +
  scale_fill_brewer(palette = "Set1", name = NULL) +
  scale_y_continuous(breaks = pretty) +
  ylab("Bias") +
  xlab("Model") +
  theme_genetics() +
  theme(legend.position = "bottom", legend.text = element_text(size = 6), legend.title = element_text(size = 8))

ggsave(filename = "model_arch_bias.jpg", plot = g_bias_example, path = fig_dir, width = 10, height = 10, units = "cm", dpi = 1000)







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


## Plot the effects of a subset. Break down by... Heritability?
g_pred_corG_list <- pred_sim_summary %>%
  filter(parameter == "corG", variable == "accuracy") %>%
  # filter(trait1_h2 != 1, trait2_h2 != 1, gencor == 0.5) %>%
  mutate(group = paste(nQTL, model, sep = "_"),
         herit = paste0("h[1]^2==~", trait1_h2, "~'/'~h[2]^2==~", trait2_h2)) %>%
  split(.$gencor) %>%
  map(~{
    ggplot(data = ., aes(x = tp_size, y = fit, color = nQTL, lty = model, shape = model, group = group)) +
      geom_point() +
      geom_line() +
      # geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.5) +
      scale_linetype_discrete(name = "Model") +
      scale_shape_discrete(name = "Model") +
      scale_color_manual(values = color_qtl) +
      scale_y_continuous(name = "Prediction accuracy", breaks = pretty, limits = c(0.15, 0.85)) +
      scale_x_discrete(name = "Training population size") +
      facet_grid(arch ~ herit, labeller = labeller(herit = label_parsed), switch = "y") +
      theme_genetics(base_size = 10) +
      theme(legend.position = "bottom", strip.placement = "outside") +
      labs(subtitle = bquote(r[G(0)]==.(unique(parse_number(.$gencor)))))
  })

## Combine
g_pred_corG_combine <- plot_grid(plotlist = map(g_pred_corG_list, ~. + theme(legend.position = "none")) , ncol = 1, labels = LETTERS[seq_along(g_pred_corG_list)])
g_pred_corG_combine1 <- plot_grid(g_pred_corG_combine, get_legend(g_pred_corG_list[[1]]), ncol = 1, rel_heights = c(1, 0.05))

ggsave(filename = "gencor_accuracy_full.jpg", plot = g_pred_corG_combine1, path = fig_dir, width = 10, height = 12, dpi = 1000)
  

## Subset for manuscript figure
g_pred_corG_paper <- pred_sim_summary %>%
  filter(parameter == "corG", variable == "accuracy") %>%
  filter(trait1_h2 != 1, trait2_h2 != 1, gencor == 0.5) %>%
  mutate(group = paste(nQTL, model, sep = "_"),
         herit = paste0("h[1]^2==~", trait1_h2, "~'/'~h[2]^2==~", trait2_h2)) %>%
  ggplot(data = ., aes(x = tp_size, y = fit, color = nQTL, lty = model, shape = model, group = group)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = nQTL), alpha = 0.15, color = 0) +
  # geom_point(size = 1) +
  geom_line(size = 0.5) +
  # geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.5) +
  scale_linetype_discrete(name = "Model") +
  scale_shape_discrete(name = "Model") +
  scale_color_manual(values = color_qtl) + 
  scale_fill_manual(values = color_qtl) +
  scale_y_continuous(name = "Prediction accuracy", breaks = pretty) +
  scale_x_discrete(name = "Training population size") +
  facet_grid(arch ~ herit, labeller = labeller(herit = label_parsed), switch = "y") +
  labs(subtitle = bquote(r[G(0)]=="0.5")) +
  theme_genetics(base_size = 8) +
  theme(legend.position = "bottom")

# Save
ggsave(filename = "gencor_accuracy_paper.jpg", path = fig_dir, width = 10, height = 10, units = "cm", dpi = 1000)


## Highlight the difference between correlation and mean/variance
g_pred_other_paper <- pred_sim_summary %>%
  filter(parameter != "musp", variable == "accuracy") %>%
  filter(gencor == 0.5, model == "RRBLUP", nQTL == 100) %>%
  mutate(group = paste(trait, parameter, sep = "_"),
         herit = paste0("h[1]^2==~", trait1_h2, "~'/'~h[2]^2==~", trait2_h2)) %>%
  ggplot(data = ., aes(x = tp_size, y = fit, color = parameter, lty = trait, group = group)) +
  geom_point() +
  geom_line() +
  # geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.5) +
  scale_linetype_discrete(name = "Trait") +
  scale_color_discrete(name = "Parameter") +
  scale_y_continuous(name = "Prediction accuracy", breaks = pretty) +
  scale_x_discrete(name = "Training population size") +
  facet_grid(arch ~ herit, labeller = labeller(herit = label_parsed)) +
  theme_presentation2(base_size = 10) +
  theme(legend.position = "bottom")



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




## Bias























### 
### Simulation results for one cycle of selection
### 

# Load the results
load(file.path(result_dir, "popvar_gencor_selection_simulation_results.RData"))

## Tidy the results
cycle1_selection_tidy <- popvar_gencor_cycle1_selection_simulation_out %>%
  select(-results, -input) %>% 
  bind_cols(., as_data_frame(transpose(popvar_gencor_cycle1_selection_simulation_out$results)))


## Model the response results
cycle1_selection_response <- cycle1_selection_tidy %>%
  unnest(response) %>%
  mutate_at(vars(trait1_h2, trait2_h2, gencor, intensity, trait), as.factor) %>%
  mutate(arch = factor(str_replace_all(arch, arch_replace), levels = arch_replace),
         selection = factor(str_replace_all(selection, selection_replace)),
         nPop = as.factor(parse_number(nPop)))

## Calculate the response of the index by averaging the response of both traits
cycle1_selection_response_index <- cycle1_selection_response %>% 
  group_by(trait1_h2, trait2_h2, gencor, arch, iter, selection, nPop, intensity) %>% 
  summarize(response = mean(response)) %>% 
  ungroup() %>%
  mutate(trait = "index")

## Calculate a mean and 95% confidence interval
cycle1_selection_response_summary <- cycle1_selection_response %>%
  bind_rows(., cycle1_selection_response_index) %>%
  select(trait1_h2:trait, response, stand_sd, cor) %>%
  gather(parameter, estimate, response, stand_sd, cor) %>%
  group_by(trait1_h2, trait2_h2, gencor, arch, selection, nPop, intensity, trait, parameter) %>%
  summarize_at(vars(estimate), funs(mean(., na.rm = T), sd(., na.rm = T), n())) %>%
  ungroup() %>%
  mutate(stat = qt(p = 1 - (alpha / 2), df = n - 1) * (sd / sqrt(n) ),
         lower = mean - stat, upper = mean + stat)



## Plot
# Response to selection
## Number of selection candidates
nCandidates <- 1200



# Value to add to trait 2 to fit on the same graph
y_nudge <- 2

g_cycle1_response <- cycle1_selection_response_summary %>%
  filter(parameter == "response", trait != "index") %>%
  mutate(intensity = parse_number(intensity)) %>%
  mutate_at(vars(mean, lower, upper), funs(ifelse(trait == "trait2", . + y_nudge, .))) %>%
  ggplot(aes(x = intensity, y = mean, shape = nPop, lty = trait)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = selection), alpha = 0.2) +
  geom_point(aes(color = selection), size = 1) +
  geom_line(aes(color = selection), lwd = 0.25) +
  scale_color_manual(values = selection_color, name = "Selection method", guide = guide_legend(nrow = 2, title.position = "top")) +
  scale_fill_manual(values = selection_color, name = "Selection method", , guide = guide_legend(nrow = 2, title.position = "top")) +
  scale_shape_discrete(name = "# Pops. / Pop. Size", labels = function(x) paste(x, "/", nCandidates / parse_number(x)),
                       guide = guide_legend(nrow = 2, title.position = "top")) +
  scale_linetype_discrete(name = "Trait", guide = guide_legend(direction = "vertical"), labels = str_to_title) +
  scale_y_continuous(breaks = pretty, name = "Response (Trait 1)", sec.axis = sec_axis(name = "Response (Trait 2)", trans = ~ . - y_nudge)) +
  scale_x_continuous(breaks = pretty, name = "Selection intensity") +
  facet_grid(trait1_h2 + trait2_h2 ~ gencor + arch) +
  theme_presentation2(base_size = 8) +
  theme(legend.position = "bottom")

ggsave(filename = "cycle1_trait_response.jpg", plot = g_cycle1_response, path = fig_dir, height = 10, width = 10, dpi = 1000)


## Just look at the correlated response and family mean
g_cycle1_response_subset <- cycle1_selection_response_summary %>%
  filter(parameter == "response", trait != "index", str_detect(selection, "Correlated|Family")) %>%
  mutate(intensity = parse_number(intensity)) %>%
  mutate_at(vars(mean, lower, upper), funs(ifelse(trait == "trait2", . + y_nudge, .))) %>%
  ggplot(aes(x = intensity, y = mean, shape = nPop, lty = trait)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = selection), alpha = 0.2) +
  geom_point(aes(color = selection), size = 1) +
  geom_line(aes(color = selection), lwd = 0.25) +
  scale_color_manual(values = selection_color, name = "Selection method", guide = guide_legend(nrow = 2, title.position = "top")) +
  scale_fill_manual(values = selection_color, name = "Selection method", , guide = guide_legend(nrow = 2, title.position = "top")) +
  scale_shape_discrete(name = "# Pops. / Pop. Size", labels = function(x) paste(x, "/", nCandidates / parse_number(x)),
                       guide = guide_legend(nrow = 2, title.position = "top")) +
  scale_linetype_discrete(name = "Trait", guide = guide_legend(direction = "vertical"), labels = str_to_title) +
  scale_y_continuous(breaks = pretty, name = "Response (Trait 1)", sec.axis = sec_axis(name = "Response (Trait 2)", trans = ~ . - y_nudge)) +
  scale_x_continuous(breaks = pretty, name = "Selection intensity") +
  facet_grid(trait1_h2 + trait2_h2 ~ gencor + arch) +
  theme_presentation2() +
  theme(legend.position = "bottom")

ggsave(filename = "cycle1_trait_response_subset.jpg", plot = g_cycle1_response_subset, path = fig_dir, height = 10, width = 10, dpi = 1000)





# Response to selection - index
g_cycle1_response_index <- cycle1_selection_response_summary %>%
  filter(parameter == "response", trait == "index") %>%
  mutate(intensity = parse_number(intensity)) %>%
  ggplot(aes(x = intensity, y = mean, shape = nPop)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = selection), alpha = 0.2) +
  # geom_point(aes(color = selection), size = 2) +
  # geom_line(aes(color = selection)) +
  geom_point(aes(color = selection), size = 0.5) +
  geom_line(aes(color = selection), lwd = 0.25) +
  scale_color_manual(values = selection_color, name = "Selection method", guide = guide_legend(nrow = 2, title.position = "top")) +
  scale_fill_manual(values = selection_color, name = "Selection method", , guide = guide_legend(nrow = 2, title.position = "top")) +
  scale_shape_discrete(name = "# Pops. / Pop. Size", labels = function(x) paste(x, "/", nCandidates / parse_number(x)),
                       guide = guide_legend(nrow = 2, title.position = "top")) +
  scale_y_continuous(breaks = pretty, name = "Response (Trait 1)") +
  scale_x_continuous(breaks = pretty, name = "Selection intensity") +
  facet_grid(trait1_h2 + trait2_h2 ~ gencor + arch, scales = "free_y") +
  theme_presentation2(base_size = 8) +
  theme(legend.position = "bottom")

ggsave(filename = "cycle1_index_response.jpg", plot = g_cycle1_response_index, path = fig_dir, height = 8, width = 10, dpi = 1000)


# Response to selection - index
# Only musp / family mean
g_cycle1_response_index_subset <- cycle1_selection_response_summary %>%
  filter(parameter == "response", trait == "index", str_detect(selection, "Correlated|Family")) %>%
  mutate(intensity = parse_number(intensity)) %>%
  ggplot(aes(x = intensity, y = mean, shape = nPop)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = selection), alpha = 0.2) +
  # geom_point(aes(color = selection), size = 2) +
  # geom_line(aes(color = selection)) +
  geom_point(aes(color = selection), size = 0.5) +
  geom_line(aes(color = selection), lwd = 0.25) +
  scale_color_manual(values = selection_color, name = "Selection method", guide = guide_legend(nrow = 2, title.position = "top")) +
  scale_fill_manual(values = selection_color, name = "Selection method", , guide = guide_legend(nrow = 2, title.position = "top")) +
  scale_shape_discrete(name = "# Pops. / Pop. Size", labels = function(x) paste(x, "/", nCandidates / parse_number(x)),
                       guide = guide_legend(nrow = 2, title.position = "top")) +
  scale_y_continuous(breaks = pretty, name = "Response (Trait 1)") +
  scale_x_continuous(breaks = pretty, name = "Selection intensity") +
  facet_grid(trait1_h2 + trait2_h2 ~ gencor + arch, scales = "free_y") +
  theme_presentation2(base_size = 8) +
  theme(legend.position = "bottom")

ggsave(filename = "cycle1_index_response_subset.jpg", plot = g_cycle1_response_index_subset, path = fig_dir, height = 8, width = 10, dpi = 1000)






# Genetic standard deviation
g_cycle1_sd <- cycle1_selection_response_summary %>%
  filter(parameter == "stand_sd", trait != "index") %>%
  mutate(intensity = parse_number(intensity)) %>%
  ggplot(aes(x = intensity, y = mean, shape = nPop)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = selection), alpha = 0.2) +
  # geom_point(aes(color = selection), size = 2) +
  # geom_line(aes(color = selection)) +
  geom_point(aes(color = selection), size = 0.5) +
  geom_line(aes(color = selection), lwd = 0.25) +
  scale_color_manual(values = selection_color, name = "Selection method", guide = guide_legend(nrow = 2, title.position = "top")) +
  scale_fill_manual(values = selection_color, name = "Selection method", , guide = guide_legend(nrow = 2, title.position = "top")) +
  scale_shape_discrete(name = "# Pops. / Pop. Size", labels = function(x) paste(x, "/", nCandidates / parse_number(x)),
                       guide = guide_legend(nrow = 2, title.position = "top")) +
  scale_y_continuous(breaks = pretty, name = "Response (Trait 1)") +
  scale_x_continuous(breaks = pretty, name = "Selection intensity") +
  facet_grid(trait + trait1_h2 + trait2_h2 ~ gencor + arch, scales = "free_y") +
  theme_presentation2(base_size = 8) +
  theme(legend.position = "bottom")


ggsave(filename = "cycle1_trait_sd.jpg", plot = g_cycle1_sd, path = fig_dir, height = 5, width = 10, dpi = 1000)


# Genetic correlation
g_cycle1_correlation <- cycle1_selection_response_summary %>%
  filter(parameter == "cor", trait == "trait1") %>%
  mutate(intensity = parse_number(intensity)) %>%
  ggplot(aes(x = intensity, y = mean, shape = nPop)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = selection), alpha = 0.2) +
  # geom_point(aes(color = selection), size = 2) +
  # geom_line(aes(color = selection)) +
  geom_point(aes(color = selection), size = 0.5) +
  geom_line(aes(color = selection), lwd = 0.25) +
  scale_color_manual(values = selection_color, name = "Selection method", guide = guide_legend(nrow = 2, title.position = "top")) +
  scale_fill_manual(values = selection_color, name = "Selection method", , guide = guide_legend(nrow = 2, title.position = "top")) +
  scale_shape_discrete(name = "# Pops. / Pop. Size", labels = function(x) paste(x, "/", nCandidates / parse_number(x)),
                       guide = guide_legend(nrow = 2, title.position = "top")) +
  scale_y_continuous(breaks = pretty, name = "Response (Trait 1)") +
  scale_x_continuous(breaks = pretty, name = "Selection intensity") +
  facet_grid(gencor ~  arch, scales = "free_y") +
  theme_presentation2(base_size = 8) +
  theme(legend.position = "bottom")


ggsave(filename = "cycle1_genetic_correlation.jpg", plot = g_cycle1_correlation, path = fig_dir, height = 5, width = 5, dpi = 1000)



### For each architecture and correlation, determine the intensity, nPop, and selection method that resulted in the highest gain
cycle1_selection_response_summary %>%
  filter(parameter == "response") %>% 
  group_by(trait, gencor, arch) %>% 
  top_n(n = 1, wt = mean) %>% 
  arrange(trait, arch)




















### Genetic correlation recurrent selection simulation


# Load the simulation results
files <- list.files(result_dir, pattern = "recurrent", full.names = TRUE)
selection_simulation_out <- vector("list", length(files))

for (i in seq_along(files)) {
  load(file = files[i])
  selection_simulation_out[[i]] <- popvar_gencor_selection_simulation_out
}

popvar_gencor_selection_simulation_out <- bind_rows(selection_simulation_out)


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
  mutate_at(vars(trait1_h2, trait2_h2, gencor, selection, arch, population), as.factor) %>%
  mutate(arch = factor(str_replace_all(arch, arch_replace), level = arch_replace),
         selection = factor(str_replace_all(selection, selection_replace), level = selection_replace))

## Separate out the haplotype frequency results
sim_allele_freq_tidy <- sim_selection_tidy %>% 
  filter(trait == "trait1", variable == "mean") %>%
  select(trait1_h2:population, contains("freq")) %>% 
  mutate(unfav_freq1 = unfav_freq,
         unfav_freq2 = 1 - (pos_freq + neg_freq + unfav_freq), 
         unfav_freq = unfav_freq + unfav_freq2) %>% ## Combined frequency of both antagonistic haplotypes
  gather(variable, response, pos_freq, neg_freq, unfav_freq1, unfav_freq2, unfav_freq) %>%
  # Calculate the change in frequency from the start
  left_join(., filter(., cycle == 0) %>% select(trait1_h2:iter, variable, base_response = response)) %>%
  mutate(response = response - base_response)


## Separate out the LD of haplotypes
sim_allele_LD_tidy <- sim_selection_tidy %>%
  select(trait1_h2:population, hap_LD) %>% 
  distinct() %>%
  gather(variable, response, hap_LD)


## Add the base population variables for the response to selection
sim_selection_response <- sim_selection_tidy %>%
  select(-contains("freq"), -hap_LD) %>%
  filter(variable == "mean") %>%
  left_join(., spread(subset(sim_selection_tidy, cycle == "0" & variable != "cor", -c(cycle, population)), variable, value)) %>%
  mutate(response = ifelse(variable == "mean", (value - mean) / sd, 1 + ((value - sd) / sd )))

## Genetic variance
sim_selection_genvar <- sim_selection_tidy %>%
  select(-contains("freq"), -hap_LD) %>%
  filter(variable == "sd") %>%
  left_join(., spread(subset(sim_selection_tidy, cycle == "0" & variable == "sd", -c(cycle, population)), variable, value)) %>% # Re-organize to show the value in the base population
  mutate(value = value^2, sd = sd^2, # Square the genetic standard deviation
         response = value / sd) # Calculate the proportion of variance remaining

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
  select(trait1_h2:trait, variable, response) %>%
  left_join(., select(rename(filter(., population == "tp"), base_cov = response), -population, -cycle), 
            by = c("trait1_h2", "trait2_h2", "gencor", "selection", "arch", "iter", "trait", "variable")) %>%
  mutate(response = scale(response))
  


## Combine
sim_selection_summ <- bind_rows(sim_selection_response, sim_selection_response_index, sim_selection_genvar,
                                rename(filter(sim_selection_tidy, variable == "cor"), response = value),
                                sim_selection_covariance, sim_allele_freq_tidy, sim_allele_LD_tidy) %>%
  filter(!is.na(response)) %>%
  filter(gencor != 0) %>% droplevels() %>% # Don't look at gencor == 0
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
summ1 <- sim_selection_summ %>% 
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
  filter(selection == "Correlated\nresponse" & selection1 != "Random") %>%
  filter(diff == max(diff)) %>%
  ungroup() %>% arrange(desc(diff))

filter(summ1, diff %in% c(min(diff), max(diff)))

# trait1_h2 trait2_h2 gencor arch          cycle selection               mean selection1        mean1    diff per_diff
# 1 0.6       0.3       0.5    Pleiotropy       10 "Correlated\nresponse"  4.16 "Family\nmean"  3.83 0.331    0.0864
# 2 0.6       0.6       -0.5   Tight Linkage     8 "Correlated\nresponse"  2.36 "Family\nmean"  2.33 0.0358   0.0154








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




## Combine index and trait response plots

## Alternate plot
g_index_response_alt <- sim_selection_summ %>% 
  filter(trait == "index", population != "parents") %>%
  filter(trait2_h2 == 0.3) %>%
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
  filter(variable == "cov", population != "parents") %>% filter(cycle != 0)

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
  filter(cycle != 0) %>%
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

# ggsave(filename = "gencor_cor_cov_paper.jpg", plot = g_variability_alt1, path = fig_dir, height = 25, width = 10, units = "cm", dpi = 1000)
ggsave(filename = "gencor_cor_cov_paper2.jpg", plot = g_variability_alt1, path = fig_dir, height = 25, width = 10, units = "cm", dpi = 1000)





haplo_replace <- c("pos_freq" = "Favorable\nCoupling", "neg_freq" = "Unfavorable\nCoupling", "unfav_freq" = "Repulsion")

## Change in haplotype frequencies
g_haplotype_list <- sim_selection_summ %>% 
  filter(str_detect(variable, "freq"), !str_detect(variable, "freq[1-2]")) %>% # Use for haplotype_frequency
  filter(!(str_detect(variable, "unfav_freq") & arch == "Pleiotropy"), population != "parents") %>%
  mutate(variable = str_replace_all(variable, haplo_replace))  %>%
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
      scale_linetype_discrete(name = "Haplotype\nphase") +
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
  filter(str_detect(variable, "freq"), !str_detect(variable, "freq[1-2]")) %>% # Use for haplotype_frequency
  filter(!(str_detect(variable, "unfav_freq") & arch == "Pleiotropy"), population != "parents", trait2_h2 != 0.3) %>%
  mutate(variable = factor(str_replace_all(variable, haplo_replace), levels = haplo_replace)) %>%
  ggplot(data = ., aes(x = cycle, y = mean, color = selection, ymin = lower, ymax = upper, fill = selection, lty = variable)) + 
  # geom_point(size = 0.5) +
  geom_line(lwd = 0.5) +
  geom_ribbon(alpha = 0.2, lwd = 0, color = 0) +
  # geom_errorbar(width = 0.5) +
  # facet_grid(herit ~ arch, labeller = labeller(herit = label_parsed), switch = "y") +
  facet_grid(cor ~ arch, labeller = labeller(herit = label_parsed, cor = label_parsed), switch = "y") +
  scale_color_manual(name = "Selection\nmethod", values = selection_color) +
  scale_fill_manual(name = "Selection\nmethod", values = selection_color) +
  scale_x_continuous(breaks = seq(0, 10, 2)) +
  scale_y_continuous(breaks = pretty) +
  scale_linetype_discrete(name = "Haplotype\nphase") +
  ylab("Change in haplotype frequency") +
  xlab("Cycle") +
  theme_genetics() +
  theme(legend.position = "bottom", legend.direction = "horizontal", legend.box = "vertical", legend.key.width = unit(0.8, "line"),
        legend.spacing.y = unit(0.25, "line"), legend.text = element_text(size = 8), legend.box.just = "left")

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








## Change in haplotype LD
g_haplotype_LD_list <- sim_selection_summ %>% 
  filter(variable == "hap_LD", population != "parents", arch != "Pleiotropy") %>%
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
      scale_y_continuous(breaks = pretty) +
      ylab("Mean haplotype linkage disequilibrium") +
      xlab("Cycle") +
      theme_genetics() +
      theme(legend.position = "bottom", legend.text = element_text(size = 8))  +
      labs(subtitle = bquote(r[G(0)]==.(unique(parse_number(.$gencor)))))
  })

g_haplotype_combine <- plot_grid(plotlist = map(g_haplotype_list, ~. + theme(legend.position = "none")), ncol = 1, labels = LETTERS[1:3])
g_haplotype_combine1 <- plot_grid(g_haplotype_combine, get_legend(g_haplotype_list[[1]]), ncol = 1, rel_heights = c(1, 0.05))

ggsave(filename = "gencor_haplotype_freq.jpg", plot = g_haplotype_combine1, path = fig_dir, height = 8, width = 4, dpi = 1000)








#######
####### Appendix
####### 
####### 

# 
# 
# ### Genetic correlation genetic architecture simulation
# # Load the simulation results
# load(file.path(result_dir, "popvar_gencor_space_simulation_results.RData"))
# # load(file.path(result_dir, "popvar_gencor_space_simulation_results_original.RData"))
# 
# 
# # Mutate the architecture space combinations
# sim_out1 <- popvar_corG_space_simulation_out %>% 
#   mutate(probcor = map(probcor, ~`names<-`(as.data.frame(.), c("dLinkage", "pLinkage")) %>% tail(., 1))) %>% # The tail is used to remove the probabilities of pleiotropy))
#   unnest(probcor) %>%
#   mutate(dLinkage = ifelse(pLinkage == 0, 0, dLinkage),
#          dLinkageFactor = ifelse(dLinkage == 0, "[0, 0]", paste0("(", round(dLinkage) - 5, ", ", dLinkage, "]"))) %>%
#   bind_cols(., as_data_frame(transpose(.$results))) %>%
#   select(-results)
# 
# ## Are there any missing combinations?
# sim_out1 %>%
#   # Remove pleiotrpic 
#   filter(pLinkage != 0) %>%
#   distinct(trait1_h2, trait2_h2, gencor, dLinkage, pLinkage, iter) %>%
#   mutate_all(as.factor) %>% 
#   mutate(obs = T) %>% 
#   complete(trait1_h2, trait2_h2, gencor, dLinkage, pLinkage, iter, fill = list(obs = F)) %>% 
#   filter(!obs)
# 
# ## Good!
# 
# # Tidy the results and extract the probability of linkage and the degree of linkage
# sim_results_tidy <- sim_out1 %>%
#   mutate_at(vars(trait1_h2:gencor, dLinkage, pLinkage), as.factor) %>%
#   mutate(dLinkageFactor = factor(dLinkageFactor, levels = c("[0, 0]", head(unique(dLinkageFactor), -1))))
# 
# 
# ## Extract each dataset
# correlation_tidy <- unnest(sim_results_tidy, other)
# predictions_tidy <- unnest(sim_results_tidy, summary)
# 
# 
# ## Extract the training population genetic correlation
# base_cor_summ <- correlation_tidy %>%
#   group_by(trait1_h2, trait2_h2, nQTL, tp_size, gencor, dLinkage, dLinkageFactor, pLinkage, variable) %>%
#   summarize_at(vars(value), funs(mean(., na.rm = T), sd(., na.rm = T), n())) %>%
#   ungroup()
# 
# base_cor_summ1 <- bind_rows(filter(base_cor_summ, pLinkage != 0),
#                             base_cor_summ %>% filter(pLinkage == 0) %>% select(-dLinkageFactor) %>% 
#                               left_join(., distinct(ungroup(base_cor_summ), gencor, dLinkageFactor)) )
# 
# 
# 
# 
# ## Plot
# g_base_cor <- base_cor_summ1 %>%
#   filter(variable == "tp_gencor") %>%
#   ggplot(aes(x = pLinkage, y = dLinkageFactor, fill = mean)) +
#   geom_tile() +
#   scale_fill_gradient2(low = "blue", high = "red") +
#   facet_wrap(~ gencor, nrow = 1) +
#   theme_acs()
# 
# # save
# ggsave(filename = "gencor_arch_space_base_corG.jpg", plot = g_base_cor, path = fig_dir,
#        height = 3, width = 6, dpi = 1000)
# 
# 
# 
# 
# ## Extract the prediction results
# # Summarize
# pred_results_summ <- predictions_tidy %>%
#   gather(variable, value, accuracy, bias) %>%
#   filter(!(variable == "bias" & abs(value) > 2)) %>%
#   group_by(trait1_h2, trait2_h2, nQTL, tp_size, gencor, dLinkage, dLinkageFactor, pLinkage, parameter, variable) %>%
#   summarize_at(vars(value), funs(mean(., na.rm = T), sd(., na.rm = T), n())) %>%
#   mutate(se = sd / sqrt(n), stat = se * qt(p = 1 - (alpha / 2), df = n - 1), lower = mean - stat, upper = mean + stat) %>%
#   ungroup()
# 
# 
# ## Plot results for genetic correlation
# g_pred_acc_corG <- pred_results_summ %>%
#   filter(parameter == "corG", variable == "accuracy") %>%
#   ggplot(aes(x = pLinkage, y = dLinkageFactor, fill = mean)) +
#   geom_tile() +
#   scale_fill_gradient(limits = c(0.40, 0.72), low = "white", high = "green", name = "Prediction\naccuracy") +
#   facet_grid(~ gencor) +
#   ylab("Maximum distance between QTL (cM)") +
#   xlab("Proportion of linked QTL pairs") +
#   theme_acs()
# 
# 
# # bias
# g_pred_bias_corG <- pred_results_summ1 %>%
#   filter(parameter == "corG", variable == "bias") %>%
#   ggplot(aes(x = pLinkage, y = dLinkageFactor, fill = mean)) +
#   geom_tile() +
#   scale_fill_gradient2(low = "red", high = "blue", name = "Bias") +
#   facet_grid(~ gencor) +
#   ylab("Maximum distance between QTL (cM)") +
#   xlab("Proportion of linked QTL pairs") +
#   theme_acs()
# 
# 
# # Combine
# g_pred_corG <- plot_grid(
#   g_pred_acc_corG + theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()), 
#   g_pred_bias_corG, 
#   ncol = 1, align = "hv")
# 
# ggsave(filename = "gencor_arch_space_pred.jpg", plot = g_pred_corG, path = fig_dir, width = 6, height = 5, dpi = 1000)
# 
# 
# 
# 
# # ## Fit a model for correlation
# # fit <- lm(value ~ gencor + dLinkage + pLinkage + gencor:pLinkage, data = correlation_tidy, subset = variable == "tp_gencor" & dLinkage != 0)
# # anova(fit)
# # plot(effects::allEffects(fit))
# # 
# # ## Notes
# # ## 1. It works
# # 
# predictions_tidy_tomodel <- predictions_tidy %>%
#   mutate(pLinkage = parse_number(pLinkage)) %>%
#   filter(parameter == "corG", dLinkage != 0)
# 
# 
# 
# 
# # For each gencor and dLinkage, fit a model regressing accuracy on pLinkage
# fit_list <- predictions_tidy_tomodel %>% 
#   group_by(gencor, dLinkage) %>% 
#   do(fit = lm(accuracy ~ pLinkage, data = .))
# 
# 
# 
# # Treat pLinkage as numeric
# fit <- lm(accuracy ~ gencor + dLinkage + pLinkage + gencor:pLinkage + dLinkage:pLinkage, data = predictions_tidy_tomodel,
#           subset = dLinkage != 50)
# anova(fit)
# plot(effects::allEffects(fit))
# 
# ## Note
# ## 1. Upward trend in prediction accuracy with decreasing pleiotropy, as expected.
# ## UPDATE: Negative trend in prediction accuracy with decreasing pleiotropy - not expected.
# 
# 
# predictions_tidy %>%
#   mutate(pLinkage = parse_number(pLinkage)) %>%
#   ggplot(aes(x = pLinkage, y = accuracy, color = gencor)) + 
#   geom_smooth(method = "lm") + 
#   facet_wrap(~ dLinkage, ncol = 5) +
#   theme_acs()
# 
# ggsave(filename = "gencor_plinkage_accuracy.jpg", plot = g_pLinkage_accuracy, path = fig_dir, width = 3, height = 3, dpi = 1000)
# 
# 
# # Plot the relationship between degree of linkage and accuracy
# # (assuming that 100% of the architecture is due to that degree of linkage)
# g_pred_linkage1 <- pred_results_summ %>% 
#   filter(pLinkage %in% c(0, 1), parameter == "corG", variable == "accuracy") %>%
#   ggplot(aes(x = dLinkageFactor, y = mean, color = gencor, group = gencor)) + 
#   geom_point() + 
#   geom_smooth(method = "lm", se = FALSE) + 
#   ylab("Accuracy") +
#   xlab("Interval between QTL for each trait (cM)") + 
#   scale_color_discrete(name = "Genetic\ncorrelation") +
#   theme_acs() +
#   theme(legend.position = c(0.80, 0.80))
# 
# ggsave(filename = "prediction_accuracy_space_linkage1.jpg", plot = g_pred_linkage1, path = fig_dir, width = 4.5, height = 4, dpi = 1000)
# 
# #
# 
# 
# 
# 
# 
