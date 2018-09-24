## PopVarValidation - analysis of predictions
##
## This script will look at the prediction output from PopVar for genetic correlations
## 
## Author: Jeff Neyhart
## Last modified: September 19, 2018
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

# Distribution of predicted corG
popvar_pred_toplot <- popvar_pred %>% 
  map(~select(., parent1:trait, family_mean = pred_mu,  contains("cor")) %>%
        gather(parameter, prediction, family_mean, contains("cor")) %>%
        filter(!is.na(prediction)))

## Separate the correlation data
popvar_pred_corG <- popvar_pred_toplot %>%
  map(~filter(., parameter != "family_mean") %>% 
        rename(trait1 = trait, trait2 = parameter) %>% 
        mutate(trait2 = str_replace(trait2, "cor_", "")) )

# Prep the family mean predictions for merging
popvar_pred_family_mean <- popvar_pred_toplot %>%
  map(~filter(., parameter == "family_mean") %>% spread(parameter, prediction))

# Merge
popvar_pred_corG_toplot <- list(popvar_pred_corG, popvar_pred_family_mean) %>%
  pmap(~left_join(x = .x, y = rename(.y, trait1_mean = family_mean), by = c("parent1", "parent2", "family", "trait1" = "trait")) %>%
         left_join(x = ., y = rename(.y, trait2_mean = family_mean), by = c("parent1", "parent2", "family", "trait2" = "trait")))

traits <- sort(unique(popvar_pred_corG_toplot$realistic$trait1))
trait_comb <- t(combn(x = traits, m = 2))

## Plot trait1 mean versus trait2 mean versus correlation
g_pred_cor <- popvar_pred_corG_toplot$realistic %>%
  filter(trait1 %in% trait_comb[1,], trait2 %in% trait_comb[2,]) %>%
  mutate(trait_pair = str_c(trait1, "_", trait2)) %>%
  split(.$trait_pair) %>%
  map(function(df) {
    df %>%
      filter(trait1 == trait1[1], trait2 == trait2[1]) %>%
      # sample_n(10000) %>%
      ggplot(aes(x = trait1_mean, y = trait2_mean, color = prediction)) + 
      geom_point(size = 1) +
      scale_color_gradient2(name = "Predicted\nCorrelation") +
      ylab(unique(df$trait2)) +
      xlab(unique(df$trait1)) +
      theme_acs()
  })

# Cowplot
g_pred_cor1 <- plot_grid(plotlist = g_pred_cor, nrow = 1)

ggsave(filename = "realistic_prediction_mean_gencor.jpg", plot = g_pred_cor1, path = fig_dir, width = 12, height = 4, dpi = 1000)



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
