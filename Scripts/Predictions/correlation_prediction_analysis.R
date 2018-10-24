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
traits <- sort(unique(popvar_pred_corG1$trait1))
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
  geom_point(aes(y = y), color = umn_palette(2)[3], size = 0.5) +
  ylab("Count") +
  xlab("Predicted genetic correlation") +
  xlim(c(-1, 1)) +
  facet_grid(~ trait_pair) +
  theme_acs()

ggsave(filename = "pred_cor_hist.jpg", plot = g_pred_cor_hist, path = fig_dir, width = 5, height = 3, dpi = 1000)


## Center and scale the means of each trait.
## Then create an index.
## Plot the index versus the predicted correlation.








## Plot trait1 mean versus trait2 mean versus correlation
g_pred_cor_mean <- popvar_pred_corG1 %>%
  filter(trait1 %in% trait_comb[,1], trait2 %in% trait_comb[,2]) %>%
  mutate(trait_pair = str_c(trait1, "_", trait2)) %>%
  split(.$trait_pair) %>%
  map(function(df) {
    df %>%
      # sample_n(10000) %>%
      ggplot(aes(x = family_mean1, y = family_mean2, color = correlation)) + 
      geom_point(size = 1) +
      scale_color_gradient2(name = "Predicted\ncorrelation") +
      ylab(paste(unique(df$trait2), "predicted family mean")) +
      xlab(paste(unique(df$trait1), "predicted family mean")) +
      theme_acs()
  })

# Cowplot
g_pred_cor1 <- plot_grid(plotlist = g_pred_cor_mean, nrow = 1)

ggsave(filename = "realistic_prediction_mean_gencor.jpg", plot = g_pred_cor1, path = fig_dir, width = 12, height = 4, dpi = 1000)


# ## Look at the relationship between covariance and correlation.
# ## How much does variance or covariance explain correlation?
# models <- popvar_pred_corG1 %>% 
#   group_by(trait1, trait2) %>% 
#   do({data_frame(
#     fit1 = list(lm(correlation ~ variance1 + variance2 + covariance, data = .)),
#     fit2 = list(lm(correlation ~ variance1 + variance2, data = .)),
#     fit3 = list(lm(correlation ~ covariance, data = .))) })


## Plot trait1 variance versus trait2 variance
g_pred_cov_var <- popvar_pred_corG1 %>%
  filter(trait1 %in% trait_comb[,1], trait2 %in% trait_comb[,2]) %>%
  mutate(trait_pair = str_c(trait1, "_", trait2)) %>%
  split(.$trait_pair) %>%
  map(function(df) {
    df %>%
      # sample_n(50000) %>%
      ggplot(aes(x = variance1, y = variance2, color = covariance)) +
      # ggplot(aes(x = variance1, y = variance2, color = correlation)) + 
      geom_point(size = 1) +
      scale_color_gradient2(name = "Predicted\ncovariance") +
      ylab(paste(unique(df$trait2), "predicted genetic variance")) +
      xlab(paste(unique(df$trait1), "predicted genetic variance")) +
      theme_acs()
  })

# Cowplot
g_pred_cov1 <- plot_grid(plotlist = g_pred_cov_var, nrow = 1)

ggsave(filename = "realistic_prediction_var_gencov.jpg", plot = g_pred_cov1, path = fig_dir, width = 12, height = 4, dpi = 1000)





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





