## GenCorPrediction - analysis of predictions
##
## This script will look at the prediction output from PopVar for genetic correlations
## 
## Author: Jeff Neyhart
## Last modified: 22 July 2019
## 
## 


# Load the source script
repo_dir <- getwd()
source(file.path(repo_dir, "source.R"))

library(cowplot)
library(modelr)

# Load the predictions
load(file.path(result_dir, "prediction_results.RData"))

# Load genotypic and phenotypic data
load(file.path(geno_dir, "S2_genos_mat.RData")) # This is available from the GitHub repository
load(file.path(pheno_dir, "PVV_BLUE.RData")) # This is available from the Triticeae Toolbox


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



### Plot distributions
### 

# Create trait combinations
trait_comb <- t(combn(x = traits, m = 2))


# Distribution of predicted corG
popvar_pred_toplot <- popvar_pred$realistic %>% 
  select(., parent1:trait, family_mean = pred_mu, variance = pred_varG, musp = musp_low, contains("cor"), contains("muspC")) %>% 
  gather(parameter, prediction, family_mean, variance, contains("cor"), contains("muspC"))  %>% 
  filter(!is.na(prediction))





## Pull out the family mean and variance predictions
popvar_pred_mu_varG <- popvar_pred_toplot %>% 
  filter(parameter %in% c("family_mean", "variance")) %>% 
  spread(parameter, prediction)


## Pull out the correlation data
popvar_pred_corG <- popvar_pred_toplot %>% 
  filter(str_detect(parameter, "cor_")) %>% 
  mutate(parameter = str_remove(parameter, "cor_")) %>% 
  select(parent1:family, trait1 = trait, trait2 = parameter, correlation = prediction)

popvar_pred_muspC <- popvar_pred_toplot %>% 
  filter(str_detect(parameter, "muspC_")) %>% 
  mutate(parameter = str_remove(parameter, "muspC_")) %>% 
  select(parent1:family, trait1 = trait, trait2 = parameter, musp, muspC = prediction)


## Mean and range of predicted corG
popvar_pred_corG %>% group_by(trait1, trait2) %>% summarize_at(vars(correlation), funs(mean, min, max))

# trait1      trait2        mean    min   max
# 1 FHBSeverity HeadingDate -0.544 -0.938 0.521
# 2 FHBSeverity PlantHeight -0.258 -0.859 0.632
# 3 HeadingDate FHBSeverity -0.544 -0.938 0.521
# 4 HeadingDate PlantHeight  0.237 -0.730 0.865
# 5 PlantHeight FHBSeverity -0.258 -0.859 0.632
# 6 PlantHeight HeadingDate  0.237 -0.730 0.865


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
  ylab("Number of crosses") +
  xlab(expression("Predicted"~italic(r[G]))) +
  scale_color_manual(values = neyhart_palette("umn1", 5)[3], name = NULL) +
  scale_y_continuous(breaks = pretty) +
  scale_x_continuous(breaks = pretty, limits = c(-1, 1)) +
  facet_grid(~ trait_pair) +
  theme_genetics() +
  theme(legend.position = c(0.90, 0.90), strip.text = element_text(size = 8))

ggsave(filename = "pred_cor_hist.jpg", plot = g_pred_cor_hist, path = fig_dir, width = 10, height = 8, unit = "cm", dpi = 1000)




# Heat colors from wes anderson package
heat_colors <- wesanderson::wes_palette("Zissou1")

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
      geom_point(size = 0.5) +
      # geom_point(size = 1) +
      scale_color_gradient2(name = expression("Predicted"~italic(hat(r)[G])), limits = c(-1, 1), 
                            low = heat_colors[5], high = heat_colors[1]) +
      ylab(bquote(.(unique(df$trait2))~predicted~mu)) +
      xlab(bquote(.(unique(df$trait1))~predicted~mu)) +
      # theme_acs()
      theme_genetics() + 
      theme(legend.position = "bottom", legend.key.width = unit(1.5, "lines"))
      
  })


# Cowplot
g_pred_cor1 <- plot_grid(plotlist = map(g_pred_cor_mean, ~. + theme(legend.position = "none")), nrow = 1)
# g_pred_cor2 <- plot_grid(g_pred_cor1, get_legend(g_pred_cor_mean[[1]]), nrow = 1, rel_widths = c(1,0.15))
g_pred_cor2 <- plot_grid( g_pred_cor1, get_legend(g_pred_cor_mean[[1]]), ncol = 1, rel_heights = c(1,0.15))


ggsave(filename = "realistic_prediction_mean_gencor_paper.jpg", plot = g_pred_cor2, path = fig_dir, width = 20, height = 8, units = "cm", dpi = 1000)



