## PopVarValidation - analysis of predictions
##
## This script will look at the prediction output from PopVar for genetic correlations
## 
## Author: Jeff Neyhart
## Last modified: October 2, 2018
## 
## 


# Load the source script
repo_dir <- getwd()
source(file.path(repo_dir, "source.R"))

library(cowplot)

# Load the predictions
load(file.path(result_dir, "prediction_results.RData"))

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



### Plot distributions
### 

# Create trait combinations
trait_comb <- t(combn(x = traits, m = 2))


# Distribution of predicted corG
popvar_pred_toplot <- popvar_pred$realistic %>% 
  select(., parent1:trait, family_mean = pred_mu, variance = pred_varG, contains("cor")) %>% 
  gather(parameter, prediction, family_mean, variance, contains("cor"))  %>% 
  filter(!is.na(prediction))


## Pull out the family mean and variance predictions
popvar_pred_mu_varG <- popvar_pred_toplot %>% 
  filter(parameter %in% c("family_mean", "variance")) %>% 
  spread(parameter, prediction)


## Pull out the correlation data
popvar_pred_corG <- popvar_pred_toplot %>% 
  filter(!parameter %in% c("family_mean", "variance")) %>% 
  mutate(parameter = str_remove(parameter, "cor_")) %>% 
  select(parent1:family, trait1 = trait, trait2 = parameter, correlation = prediction)


## Add the mean and variance data back in
## Then calculate the covariance
popvar_pred_corG1 <- popvar_pred_corG %>% 
  left_join(., rename_at(popvar_pred_mu_varG, vars(family_mean, variance), ~paste0(., "1")), 
            by = c("parent1", "parent2", "family", "trait1" = "trait")) %>% 
  left_join(., rename_at(popvar_pred_mu_varG, vars(family_mean, variance), ~paste0(., "2")), 
            by = c("parent1", "parent2", "family", "trait2" = "trait")) %>%
  mutate(covariance = correlation * (sqrt(variance1) * sqrt(variance2))) %>%
  filter(trait1 %in% trait_comb[,1], trait2 %in% trait_comb[,2])



## Plot distributions of the correlations between traits
g_pred_cor_hist <- popvar_pred_corG1 %>%
  mutate(trait_pair = str_c(trait1, " / ", trait2),
         y = ifelse(is.na(family), NA, 10000)) %>%
  ggplot(aes(x = correlation)) +
  geom_histogram() +
  geom_point(aes(y = y, color = "Selected cross"), size = 2) +
  ylab("Count") +
  xlab(expression("Predicted"~r[A])) +
  scale_color_manual(values = neyhart_palette("umn1", 5)[3], name = NULL) +
  scale_y_continuous(breaks = pretty) +
  scale_x_continuous(breaks = pretty, limits = c(-1, 1)) +
  facet_grid(~ trait_pair) +
  theme_presentation2() +
  theme(legend.position = c(0.78, 0.87))

ggsave(filename = "pred_cor_hist.jpg", plot = g_pred_cor_hist, path = fig_dir, width = 10, height = 4, dpi = 1000)



## Plot an example for a presentation
g_pred_cor_hist_example <- popvar_pred_corG1 %>%
  filter(trait1 == "FHBSeverity", trait2 == "HeadingDate") %>%
  mutate_at(vars(trait1, trait2), funs(str_replace_all(., traits_replace))) %>%
  mutate(trait_pair = str_c(trait1, " / ", trait2),
         y = ifelse(is.na(family), NA, 10000)) %>%
  ggplot(aes(x = correlation)) +
  # geom_density() +
  geom_histogram() +
  facet_grid(~ trait_pair) +
  scale_y_continuous(breaks = pretty, name = "Number of potential crosses") +
  scale_x_continuous(breaks = pretty, limits = c(-1, 1), name = expression(Predicted~italic(r[G]))) + 
  theme_presentation2(base_size = 18)

# Save
ggsave(filename = "pred_cor_hist_example.jpg", plot = g_pred_cor_hist_example, path = fig_dir,
       height = 6, width = 6, dpi = 1000)

## Blank plots for highest and lowest correlation
popvar_pred_corG1_example <- popvar_pred_corG1 %>%
  filter(trait1 == "FHBSeverity", trait2 == "HeadingDate") %>%
  filter(correlation == max(correlation) | correlation == min(correlation)) %>%
  mutate_at(vars(trait1, trait2), funs(str_replace_all(., traits_replace)))

# For each population, simulate a bi-variate distribution of breeding values
popvar_pred_corG1_example_sim <- popvar_pred_corG1_example %>%
  group_by(parent1, parent2) %>%
  do(plot = {
    cross <- .
    
    mu <- c(cross$family_mean1, cross$family_mean2)
    sigma <- rbind(c(cross$variance1, cross$covariance), c(cross$covariance, cross$variance2))
    bv <- mvtnorm::rmvnorm(n = 150, mean = mu, sigma = sigma)
    
    ## Convert to df
    as_data_frame(bv) %>% 
      ggplot(aes(x = V1, y = V2)) +
      geom_point() +
      scale_x_continuous(name = cross$trait1, labels = NULL) +
      scale_y_continuous(name = cross$trait2, labels = NULL) +
      theme_classic(base_size = 16)
    
  })

## Save
g_example_sim <- plot_grid(plotlist = popvar_pred_corG1_example_sim$plot, nrow = 1)
ggsave(filename = "pred_cor_example.jpg", plot = g_example_sim, path = fig_dir, width = 8, height = 3, dpi = 1000)




## Plot trait1 mean versus trait2 mean versus correlation
g_pred_cor_mean <- popvar_pred_corG1 %>%
  filter(trait1 %in% trait_comb[,1], trait2 %in% trait_comb[,2]) %>%
  mutate(trait_pair = str_c(trait1, "_", trait2)) %>% 
  mutate_at(vars(contains("trait")), funs(str_replace_all(., traits_replace))) %>%
  split(.$trait_pair) %>%
  map(function(df) {
    df %>%
      # sample_n(10000) %>%
      ggplot(aes(x = family_mean1, y = family_mean2, color = correlation)) + 
      # geom_point(size = 0.5) +
      geom_point(size = 1) +
      scale_color_gradient2(name = expression("Predicted"~r[G]), limits = c(-1, 1)) +
      ylab(bquote(.(unique(df$trait2))~predicted~mu)) +
      xlab(bquote(.(unique(df$trait1))~predicted~mu)) +
      # theme_acs()
      theme_presentation2() + 
      theme(legend.position = "top", legend.direction = "horizontal", legend.justification = "right",
            legend.key.width = unit(1.5, "lines"))
      
  })


# Cowplot
g_pred_cor1 <- plot_grid(plotlist = map(g_pred_cor_mean, ~. + theme(legend.position = "none")), nrow = 1)
# g_pred_cor2 <- plot_grid(g_pred_cor1, get_legend(g_pred_cor_mean[[1]]), nrow = 1, rel_widths = c(1,0.15))
g_pred_cor2 <- plot_grid( get_legend(g_pred_cor_mean[[1]]), g_pred_cor1, ncol = 1, rel_heights = c(0.15,1))


# ggsave(filename = "realistic_prediction_mean_gencor.jpg", plot = g_pred_cor2, path = fig_dir, width = 8, height = 2.5, dpi = 1000)
ggsave(filename = "realistic_prediction_mean_gencor_presentation.jpg", plot = g_pred_cor2, path = fig_dir, width = 12, height = 4.5, dpi = 1000)




## Plot trait1 variance versus trait2 variance
g_pred_cov_var <- popvar_pred_corG1 %>%
  filter(trait1 %in% trait_comb[,1], trait2 %in% trait_comb[,2]) %>%
  mutate(trait_pair = str_c(trait1, "_", trait2)) %>%
  split(.$trait_pair) %>%
  map(function(df) {
    df %>%
      # sample_n(10000) %>%
      ggplot(aes(x = variance1, y = variance2, color = covariance)) +
      # ggplot(aes(x = variance1, y = variance2, color = correlation)) + 
      geom_point(size = 0.5) +
      scale_color_gradient2(name = expression("Predicted"~Cov[G]), limits = c(-1.5, 1.5)) +
      ylab(bquote(.(unique(df$trait2))~predicted~sigma[G]^2)) +
      xlab(bquote(.(unique(df$trait1))~predicted~sigma[G]^2)) +
      theme_acs()
  })

# Cowplot
g_pred_cov1 <- plot_grid(plotlist = map(g_pred_cov_var, ~. + theme(legend.position = "none")), nrow = 1)
g_pred_cov2 <- plot_grid(g_pred_cov1, get_legend(g_pred_cov_var[[1]]), nrow = 1, rel_widths = c(1,0.15))

ggsave(filename = "realistic_prediction_var_gencov.jpg", plot = g_pred_cov2, path = fig_dir, width = 8, height = 2.5, dpi = 1000)



## Look at the relationship between covariance and correlation.
## How much does variance or covariance explain correlation?
models <- popvar_pred_corG1 %>%
  group_by(trait1, trait2) %>%
  do({data_frame(
    fit1 = list(lm(correlation ~ variance1 + variance2 + covariance, data = .)),
    fit2 = list(lm(correlation ~ variance1 + variance2, data = .)),
    fit3 = list(lm(correlation ~ covariance, data = .))) })




## Plot FHB and heading date with associated predictions of variance
## First plot trait1 mean versus trait2 mean versus correlation
g_pred_cor_mean_example <- popvar_pred_corG1 %>%
  filter(trait1 == "FHBSeverity", trait2 == "HeadingDate") %>%
  mutate(trait_pair = str_c(trait1, "_", trait2)) %>% 
  mutate_at(vars(contains("trait")), funs(str_replace_all(., traits_replace))) %>%
  # sample_n(10000) %>%
  ggplot(aes(x = family_mean1, y = family_mean2, color = correlation)) + 
  # geom_point(size = 0.5) +
  geom_point(size = 1) +
  scale_color_gradient2(name = expression(hat(r)[G]), limits = c(-1, 1)) +
  ylab(expression(Heading~Date~hat(mu))) +
  xlab(expression(FHB~Severity~hat(mu))) +
  scale_x_continuous(breaks = pretty, position = "top") +
  scale_y_continuous(breaks = pretty) +
  # theme_acs()
  theme_presentation2() + 
  theme(legend.position = "left", legend.direction = "vertical", legend.justification = "center",
        legend.key.width = unit(1.5, "lines"))

## Plot mean versus variance for heading date
g_pred_var_mean_HD <- popvar_pred_mu_varG %>%
  filter(trait == "HeadingDate") %>%
  # sample_n(10000) %>%
  ggplot(aes(x = family_mean, y = variance)) +
  geom_point(size = 1) + 
  ylab(expression("Heading Date"~hat(sigma)[G]^2)) +
  scale_y_continuous(position = "right") +
  coord_flip() +
  # scale_y_reverse() +
  theme_presentation2() +
  theme(axis.title.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(),
        axis.line.y = element_blank(), panel.border = element_blank())

## Plot mean versus variance for FHB severity
g_pred_var_mean_FHB <- popvar_pred_mu_varG %>%
  filter(trait == "FHBSeverity") %>%
  # sample_n(10000) %>%
  ggplot(aes(x = family_mean, y = variance)) +
  geom_point(size = 1) + 
  ylab(expression("FHB Severity "~hat(sigma)[G]^2)) +
  scale_y_reverse() +
  theme_presentation2() +
  theme(axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank(),
        axis.line.x = element_blank(), panel.border = element_blank())
  

# Try patchwork
library(patchwork)

g_combine <- ((g_pred_cor_mean_example + g_pred_var_mean_HD) + g_pred_var_mean_FHB) + 
  plot_layout(nrow = 2, widths = c(1, 0.60), heights = c(1, 0.6))

# Save
ggsave(filename = "cor_mean_var_example_predictions.jpg", plot = g_combine, path = fig_dir, width = 8, height = 6, dpi = 1000)



# Mean and range of predictions
popvar_pred_summ <- popvar_pred_corG_toplot %>% 
  map(~group_by(., trait1, trait2) %>% 
        summarize_at(vars(prediction), funs(min, max, mean)))

# $`realistic`
# trait1      trait2         min   max   mean
# 1 FHBSeverity HeadingDate -0.938 0.521 -0.544
# 2 FHBSeverity PlantHeight -0.859 0.632 -0.258
# 3 HeadingDate FHBSeverity -0.938 0.521 -0.544
# 4 HeadingDate PlantHeight -0.730 0.865  0.237
# 5 PlantHeight FHBSeverity -0.859 0.632 -0.258
# 6 PlantHeight HeadingDate -0.730 0.865  0.237
# 
# $relevant
# trait1      trait2         min   max   mean
# 1 FHBSeverity HeadingDate -0.937 0.416 -0.569
# 2 FHBSeverity PlantHeight -0.865 0.654 -0.289
# 3 HeadingDate FHBSeverity -0.937 0.416 -0.569
# 4 HeadingDate PlantHeight -0.717 0.862  0.241
# 5 PlantHeight FHBSeverity -0.865 0.654 -0.289
# 6 PlantHeight HeadingDate -0.717 0.862  0.241


### Generate the same summary for the selected crosses
popvar_pred_summ_selected <- popvar_pred_corG_toplot %>% 
  map(~filter(., !is.na(family)) %>%
        group_by(., trait1, trait2) %>% 
        summarize_at(vars(prediction), funs(min, max, mean)))

# $`realistic`
# trait1      trait2         min    max   mean
# 1 FHBSeverity HeadingDate -0.793 -0.235 -0.563
# 2 FHBSeverity PlantHeight -0.710  0.227 -0.289
# 3 HeadingDate FHBSeverity -0.793 -0.235 -0.563
# 4 HeadingDate PlantHeight -0.287  0.696  0.216
# 5 PlantHeight FHBSeverity -0.710  0.227 -0.289
# 6 PlantHeight HeadingDate -0.287  0.696  0.216
# 
# $relevant
# trait1      trait2         min      max   mean
# 1 FHBSeverity HeadingDate -0.793 -0.420   -0.611
# 2 FHBSeverity PlantHeight -0.661  0.00646 -0.311
# 3 HeadingDate FHBSeverity -0.793 -0.420   -0.611
# 4 HeadingDate PlantHeight -0.241  0.664    0.244
# 5 PlantHeight FHBSeverity -0.661  0.00646 -0.311
# 6 PlantHeight HeadingDate -0.241  0.664    0.244



## Calculate the superior progeny mean and correlated superior progeny mean
popvar_pred_corG1_musp <- popvar_pred_corG1 %>% 
  mutate(pred_musp1 = family_mean1 - (k_sp * sqrt(variance1)), pred_musp2 = family_mean2 - (k_sp * sqrt(variance2)), 
         pred_musp1C = family_mean1 - (k_sp * correlation * sqrt(variance1)), pred_musp2C = family_mean2 - (k_sp * correlation * sqrt(variance2))) %>%
  filter(trait1 %in% trait_comb[,1], trait2 %in% trait_comb[,2]) %>%
  mutate(trait_pair = str_c(trait1, "_", trait2))



## Fit models
fit_muspC <- popvar_pred_corG1_musp %>% 
  split(.$trait_pair) %>%
  map(~lm(pred_musp2C ~ correlation, data = .)) %>%
  map(summary)





