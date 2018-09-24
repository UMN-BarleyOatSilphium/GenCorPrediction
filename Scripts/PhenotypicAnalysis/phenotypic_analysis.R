## PopVarValidation
## Phenotypic analysis
##
## This script will run statistical analyses of the phenotypic data from the PopVarVal project. This
## will include:
## 1. Phenotypic analysis of training population data
## 2. Phenotypic analysis of validation family data
## 
## Author: Jeff Neyhart
## Last modified: August 19, 2018
## 

# Load the source script
repo_dir <- getwd()
source(file.path(repo_dir, "source.R"))

# Load the pbr library
library(pbr)
library(broom)
library(ggridges)
library(cowplot)


## Load the S2 BLUEs
load(file.path(gdrive_dir, "BarleyLab/Breeding/PhenotypicData/Final/MasterPhenotypes/S2_tidy_BLUE.RData"))




### Phenotypic analysis of training population data

# First gather data that would have been used to make the predictions
tp_prediction_tomodel <- s2_tidy_BLUE %>% 
  filter(trait %in% traits, line_name %in% tp, year %in% 2014:2015,
         location %in% c("STP", "CRM")) %>%
  arrange(year, location, trial, trait, line_name)

# Now gather data on the training population that would be relevant to the predictions, regardless of year
tp_relevant_tomodel <- s2_tidy_BLUE %>%
  filter(trait %in% traits, line_name %in% tp, year %in% 2014:2017, 
         location %in% c("STP", "CRM", "FND", "BCW"))


## Run models
tp_analysis <- tp_prediction_tomodel %>%
  group_by(trait) %>%
  do({
    df <- .
    print(unique(df$trait))
    summarize_pheno(data = df)
  })




## Look at variance components and heritability
(g_tp_prediction_h2 <- tp_analysis %>% 
  mutate(h2 = map_dbl(h2, "heritability")) %>% 
  qplot(x = trait, y = h2, geom = "col", fill = "blue", data = .) +
  geom_text(aes(label = str_c("Envs: ", n_e)), vjust = 2) +
  geom_text(aes(label = str_c("h2: ", round(h2, 3))), vjust = 4) + 
  scale_fill_discrete(guide = FALSE))

ggsave(filename = "tp_prediction_h2.jpg", plot = g_tp_prediction_h2, path = fig_dir, height = 5, width = 5, dpi = 1000)
  

tp_var_prop <- tp_analysis %>% 
  mutate(varcomp = map(h2, "var_comp")) %>%
  unnest(varcomp) %>% 
  group_by(trait) %>% 
  mutate(var_prop = variance / sum(variance)) 

# trait       source                variance  var_prop
# 1 FHBSeverity line_name:environment 49.2     0.825    
# 2 FHBSeverity line_name             10.4     0.175    
# 3 FHBSeverity Residual               0.00456 0.0000764
# 4 HeadingDate line_name:environment  0.725   0.0665   
# 5 HeadingDate line_name              9.20    0.844    
# 6 HeadingDate Residual               0.972   0.0892   
# 7 PlantHeight line_name:environment 19.5     0.642    
# 8 PlantHeight line_name             10.9     0.358    
# 9 PlantHeight Residual               0.00245 0.0000808

(g_tp_prediction_varprop <- tp_var_prop %>% 
  ggplot(aes(x = trait, y = var_prop, fill = source)) + 
  geom_col(position = "dodge"))

ggsave(filename = "tp_prediction_varprop.jpg", plot = g_tp_prediction_varprop, path = fig_dir, height = 5, width = 5, dpi = 1000)



tp_analysis %>%
  unnest(sig_test) %>%
  mutate(annotate = case_when(p_value <= 0.01 ~ "***", p_value <= 0.05 ~ "**", p_value <= 0.1 ~ "*", TRUE ~ ""))

# trait         n_e term     df statistic   p_value
# 1 FHBSeverity     4 g         1     26.5  2.59e-  7 ***     
# 2 FHBSeverity     4 ge        1     75.8  3.13e- 18 ***     
# 3 HeadingDate     4 g         1    721.   8.07e-159 ***     
# 4 HeadingDate     4 ge        1      5.87 1.54e-  2 **      
# 5 PlantHeight     2 g         1     24.9  6.00e-  7 ***     
# 6 PlantHeight     2 ge        1      4.67 3.07e-  2 ** 

## G and GxE are significant for all traits, though more so for FHB severity.

# Unnest the blues
tp_prediction_BLUE <- tp_analysis %>%
  unnest(BLUE) %>%
  ungroup() %>%
  select(-n_e)



# ## Do the same thing for the relevant tp data
# 
# ## Run models
# tp_analysis <- tp_relevant_tomodel %>%
#   group_by(trait) %>%
#   do({
#     df <- .
#     print(unique(df$trait))
#     summarize_pheno(data = df)
#   })
# 
# 
# ## Look at variance components and heritability
# tp_analysis %>% select(trait, h2, n_e) %>% mutate(h2 = map_dbl(h2, "heritability"))
# 
# # trait          h2   n_e
# # 1 FHBSeverity 0.446     4
# # 2 HeadingDate 0.972    10
# # 3 PlantHeight 0.746     9
# 
# (g_tp_relevant_h2 <- tp_analysis %>% 
#     mutate(h2 = map_dbl(h2, "heritability")) %>% 
#     qplot(x = trait, y = h2, geom = "col", fill = "blue", data = .) +
#     geom_text(aes(label = str_c("Envs: ", n_e)), vjust = 2) +
#     geom_text(aes(label = str_c("h2: ", round(h2, 3))), vjust = 4) + 
#     scale_fill_discrete(guide = FALSE))
# 
# ggsave(filename = "tp_relevant_h2.jpg", plot = g_tp_relevant_h2, path = fig_dir, height = 5, width = 5, dpi = 1000)
# 
# tp_var_prop <- tp_analysis %>% 
#   mutate(varcomp = map(h2, "var_comp")) %>%
#   unnest(varcomp) %>% 
#   group_by(trait) %>% 
#   mutate(var_prop = variance / sum(variance)) 
# 
# # trait         n_e source                   variance    var_prop
# # 1 FHBSeverity     4 line_name:environment 49.2        0.825      
# # 2 FHBSeverity     4 line_name             10.4        0.175      
# # 3 FHBSeverity     4 Residual               0.00456    0.0000764  
# # 4 HeadingDate    10 line_name:environment  2.98       0.222      
# # 5 HeadingDate    10 line_name             10.4        0.778      
# # 6 HeadingDate    10 Residual               0.00000856 0.000000638
# # 7 PlantHeight     9 line_name:environment 22.2        0.752      
# # 8 PlantHeight     9 line_name              7.33       0.248      
# # 9 PlantHeight     9 Residual               0.000507   0.0000172
# 
# 
# g_tp_relevant_varprop <- tp_var_prop %>% 
#   ggplot(aes(x = trait, y = var_prop, fill = source)) + 
#   geom_col(position = "dodge")
# 
# ggsave(filename = "tp_relevant_varprop.jpg", plot = g_tp_relevant_varprop, path = fig_dir, height = 5, width = 5, dpi = 1000)
# 
# 
# 
# tp_analysis %>%
#   unnest(sig_test)  %>%
#   mutate(annotate = case_when(p_value <= 0.01 ~ "***", p_value <= 0.05 ~ "**", p_value <= 0.1 ~ "*", TRUE ~ ""))
# 
# ## G and GxE are significant for all traits, though more so for FHB severity.
# 
# # trait         n_e term     df statistic   p_value annotate
# # 1 FHBSeverity     4 g         1      26.5 2.59e-  7 ***     
# # 2 FHBSeverity     4 ge        1      75.8 3.13e- 18 ***     
# # 3 HeadingDate    10 g         1    2066.  0.        ***     
# # 4 HeadingDate    10 ge        1     553.  3.27e-122 ***     
# # 5 PlantHeight     9 g         1     215.  1.30e- 48 ***     
# # 6 PlantHeight     9 ge        1     538.  5.42e-119 *** 
# 
# 
# 
# # Unnest the blues
# tp_relevant_BLUE <- tp_analysis %>%
#   unnest(BLUE) %>%
#   ungroup() %>%
#   select(-n_e)




### Explore FHB severity

## First fit GxYxL models
tp_FHB_tomodel <- tp_prediction_tomodel %>%
  filter(trait == "FHBSeverity")
  
control <- lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.nRE = "ignore")
wts <- tp_FHB_tomodel$std.error^2

# Vector of random effects
randos <- c("(1|line_name)", "(1|location)", "(1|year)", "(1|line_name:location)", "(1|line_name:year)", "(1|line_name:location:year)")

# Full model
full <- as.formula(paste("value ~ ", paste(randos, collapse = "+")))
# Drop1
dropped <- set_names(randos, randos) %>%
  map(~setdiff(randos, .)) %>%
  map(~paste(., collapse = "+")) %>%
  map(~as.formula(paste("value ~ ", paste(., collapse = "+")), env = .GlobalEnv))

forms <- c(full, dropped)

## Fit the models
fits <- modelr::fit_with(data = tp_FHB_tomodel, .f = lmer, .formulas = forms, control = control, weights = wts)
 
## Compare these reduced models to the full model and perform a LRT
lrt_out <- map(fits[-1], ~lr_test(model1 = fits[[1]], model2 = .)) %>% 
  list(., names(.)) %>% 
  pmap_df(~mutate(.x, term = .y) %>% select(term, names(.))) %>%
  mutate(p_value = formatC(p_value, digits = 2, format = "g")) %>%
  as_data_frame()

# term                        full_model    df    statistic p_value
# 1 (1|line_name)               model1         1   9.13       0.0025 
# 2 (1|location)                model1         1  15.8        7.1e-05
# 3 (1|year)                    model1         1   2.94       0.087  
# 4 (1|line_name:location)      model1         1   0.844      0.36   
# 5 (1|line_name:year)          model1         1  -0.00000833 "  1"  
# 6 (1|line_name:location:year) model1         1  59.2        1.4e-14


## Genotype, location, and GxYxL are all significant.
## The significance of location and the non-significance of year suggest that something
## about the location is driving changes, and this very well may be management.




#### Separate analysis for CRM and STP for FHB
## Run models
tp_analysis_FHB <- tp_prediction_tomodel %>%
  filter(trait %in% c("FHBSeverity", "HeadingDate")) %>%
  group_by(location, trait) %>%
  do({
    df <- .
    print(unique(df$trait))
    summarize_pheno(data = df)
  })


## Look at variance components and heritability
(g_tp_prediction_h2 <- tp_analysis_FHB %>% 
    filter(trait == "FHBSeverity") %>%
    mutate(h2 = map_dbl(h2, "heritability")) %>% 
    qplot(x = location, y = h2, geom = "col", fill = "blue", data = .) +
    geom_text(aes(label = str_c("Envs: ", n_e)), vjust = 2) +
    geom_text(aes(label = str_c("h2: ", round(h2, 3))), vjust = 4) + 
    scale_fill_discrete(guide = FALSE))

ggsave(filename = "tp_prediction_h2_FHB.jpg", plot = g_tp_prediction_h2, path = fig_dir, height = 5, width = 5, dpi = 1000)


## Calculate the correlation between FHB and HD
tp_FHB_HD <- tp_analysis_FHB %>% 
  unnest(BLUE) %>% 
  select(location, line_name, trait, value) %>% 
  spread(trait, value) %>% 
  group_by(location)



tp_FHB_HD %>% 
  summarize(cor = cor(FHBSeverity, HeadingDate))

# location    cor
# 1 CRM      -0.600
# 2 STP      -0.233

qplot(data = tp_FHB_HD, x = FHBSeverity, y = HeadingDate, facets = ~location) + theme_acs()



tp_var_prop <- tp_analysis_FHB %>% 
  mutate(varcomp = map(h2, "var_comp")) %>%
  unnest(varcomp) %>% 
  group_by(location) %>% 
  mutate(var_prop = variance / sum(variance)) 

# location   n_e source                variance var_prop
# 1 CRM          2 line_name:environment 45.4     0.687   
# 2 CRM          2 line_name             20.6     0.312   
# 3 CRM          2 Residual               0.0232  0.000352
# 4 STP          2 line_name:environment 47.3     0.890   
# 5 STP          2 line_name              5.87    0.110   
# 6 STP          2 Residual               0.00634 0.000119

(g_tp_prediction_varprop <- tp_var_prop %>% 
    ggplot(aes(x = location, y = var_prop, fill = source)) + 
    geom_col(position = "dodge"))

ggsave(filename = "tp_prediction_varprop_FHB.jpg", plot = g_tp_prediction_varprop, path = fig_dir, height = 5, width = 5, dpi = 1000)



tp_analysis_FHB %>%
  unnest(sig_test) %>%
  mutate(annotate = case_when(p_value <= 0.01 ~ "***", p_value <= 0.05 ~ "**", p_value <= 0.1 ~ "*", TRUE ~ ""))

# location   n_e term     df statistic      p_value annotate
# 1 CRM          2 g         1     17.6  0.0000280    ***     
# 2 CRM          2 ge        1     32.4  0.0000000124 ***     
# 3 STP          2 g         1      2.03 0.154        ""      
# 4 STP          2 ge        1     20.8  0.00000512   *** 

## G and GxE are significant for all traits, though more so for FHB severity.

# Unnest the blues
tp_prediction_BLUE_FHB <- tp_analysis_FHB %>%
  unnest(BLUE) %>%
  ungroup() %>%
  select(-n_e)

## Heritability
tp_analysis_FHB %>% 
  group_by(location) %>% 
  summarize(h2 = map_dbl(h2, "heritability"))

# location    h2
# 1 CRM      0.464
# 2 STP      0.189



## TP mean and range
tp_prediction_BLUE %>% 
  group_by(trait) %>% 
  summarize_at(vars(value), funs(min, max, mean))

# trait         min   max  mean
# 1 FHBSeverity  5.08  38.6  16.1
# 2 HeadingDate 43.9   56.5  50.2
# 3 PlantHeight 60.1   87.4  73.9

tp_prediction_BLUE_FHB %>% 
  group_by(location) %>% 
  summarize_at(vars(value), funs(min, max, mean))

# trait         min   max  mean
# 1 FHBSeverity  5.08  38.6  16.1
# 2 HeadingDate 43.9   56.5  50.2
# 3 PlantHeight 60.1   87.4  73.9




### Phenotypic analysis of validation family data
# Gather data

vp_family_tomodel <- s2_tidy_BLUE %>%
  filter(trait %in% traits, line_name %in% c(pot_pars, exper), str_detect(trial, "PVV"))


## Run models
vp_analysis <- vp_family_tomodel %>%
  group_by(trait) %>%
  do({
    df <- .
    print(unique(df$trait))
    summarize_pheno(data = df, blue.model = "sommer")
  })


## Look at variance components and heritability
vp_analysis %>% 
  select(trait, h2, n_e) %>% 
  mutate(h2 = map_dbl(h2, "heritability"))

# trait             h2   n_e
# 1 FHBSeverity 0.114     4
# 2 HeadingDate 0.784     4
# 3 PlantHeight 0.740     4

# Heritability is low for FHB, but adequate for HD and PH.

(g_vp_h2 <- vp_analysis %>% 
    mutate(h2 = map_dbl(h2, "heritability")) %>% 
    qplot(x = trait, y = h2, geom = "col", fill = "blue", data = .) +
    geom_text(aes(label = str_c("Envs: ", n_e)), vjust = 2) +
    geom_text(aes(label = str_c("h2: ", round(h2, 3))), vjust = 4) + 
    scale_fill_discrete(guide = FALSE))

ggsave(filename = "vp_h2.jpg", plot = g_vp_h2, path = fig_dir, height = 5, width = 5, dpi = 1000)


vp_var_prop <- vp_analysis %>% 
  mutate(varcomp = map(h2, "var_comp")) %>%
  unnest(varcomp) %>% 
  group_by(trait) %>% 
  mutate(var_prop = variance / sum(variance)) 

# trait         n_e source                 variance  var_prop
# 1 FHBSeverity     4 line_name:environment 137.      0.959    
# 2 FHBSeverity     4 line_name               5.83    0.0408   
# 3 FHBSeverity     4 Residual                0.00805 0.0000562
# 4 HeadingDate     4 line_name:environment   4.43    0.309    
# 5 HeadingDate     4 line_name               8.38    0.584    
# 6 HeadingDate     4 Residual                1.53    0.107    
# 7 PlantHeight     4 line_name:environment  19.5     0.465    
# 8 PlantHeight     4 line_name              22.5     0.535    
# 9 PlantHeight     4 Residual                0.00377 0.0000896


g_vp_varprop <- vp_var_prop %>% 
  ggplot(aes(x = trait, y = var_prop, fill = source)) + 
  geom_col(position = "dodge")

ggsave(filename = "vp_varprop.jpg", plot = g_vp_varprop, path = fig_dir, height = 5, width = 5, dpi = 1000)


vp_analysis %>%
  unnest(sig_test)  %>%
  mutate(annotate = case_when(p_value <= 0.01 ~ "***", p_value <= 0.05 ~ "**", p_value <= 0.1 ~ "*", TRUE ~ ""))

## G and GxE are significant for all traits

# trait         n_e term     df statistic   p_value annotate
# 1 FHBSeverity     4 g         1      10.4 1.28e-  3 ***
# 2 FHBSeverity     4 ge        1    1775.  0.        ***
# 3 HeadingDate     4 g         1    3499.  0.        ***
# 4 HeadingDate     4 ge        1      41.9 9.38e- 11 ***
# 5 PlantHeight     4 g         1    1698.  0.        ***
# 6 PlantHeight     4 ge        1     558.  2.40e-123 ***
# 


# Unnest the blues
vp_BLUE <- vp_analysis %>%
  unnest(BLUE) %>%
  ungroup() %>%
  select(-n_e)

## For FHB severity, try estimating variance components separately by location
vp_analysis_FHB <- vp_family_tomodel %>%
  filter(trait == "FHBSeverity") %>%
  group_by(location) %>%
  do(summarize_pheno(data = ., blue.model = "sommer"))

# Look at heritability
vp_analysis_FHB %>% 
  mutate(h2 = map_dbl(h2, "heritability")) %>%
  select(location, n_e, h2)

# location   n_e     h2
# 1 CRM          2 0.267 
# 2 STP          2 0.0417


# Is GE still significant when analyzing locations independently?
vp_analysis_FHB %>% unnest(sig_test)

# location   n_e term     df statistic   p_value
# 1 CRM          2 g         1    33.2   8.46e-  9
# 2 CRM          2 ge        1   332.    4.26e- 74
# 3 STP          2 g         1     0.755 3.85e-  1
# 4 STP          2 ge        1   850.    7.16e-187

## Yes


# Extract the BLUEs
vp_FHB_BLUE <- vp_analysis_FHB %>%
  unnest(BLUE) %>% 
  select(-n_e)


## Save the BLUEs
save("tp_prediction_BLUE", "tp_relevant_BLUE", "tp_prediction_BLUE_FHB", "vp_BLUE", "vp_FHB_BLUE", 
     file = file.path(data_dir, "PVV_BLUE.RData"))

# Load
load(file.path(data_dir, "PVV_BLUE.RData"))










## Other plots


# Plot the distributions of each trait and for each family
vp_family_BLUE <- vp_BLUE %>%
  filter(line_name %in% exper) %>%
  mutate(family = str_extract(line_name, "4[0-9]{3}")) %>%
  group_by(family, trait) %>%
  mutate(family_mean_est = mean(value)) %>%
  ungroup()


## Calculate the mean, min, and max for all traits for the VP
vp_family_BLUE %>% 
  group_by(trait) %>% 
  summarize(mean = mean(value), min = min(value), max = max(value))

# trait        mean   min   max
# 1 FHBSeverity  21.1  4.19  56.3
# 2 HeadingDate  53.5 44.2   66.6
# 3 PlantHeight 101.  82.2  117.




# Combine
vp_family_BLUE_toplot <- bind_rows(
  mutate(vp_family_BLUE, type = "By family"), 
  mutate(vp_family_BLUE, family = "All", type = "All"))

# Plot per trait
g_vp_family_density <- vp_family_BLUE_toplot %>%
  split(.$trait) %>%
  map(~{
    temp <- .
    # Order the family based on the family mean
    family_order <- distinct(temp, family, family_mean_est) %>% 
      filter(family != "All") %>% 
      arrange(family_mean_est) %>% 
      pull(family)
    
    # Convert family to factor
    temp$family <- factor(temp$family, levels = c(family_order, "All"))
    
    # Create a color scheme
    family_color <- setNames(c(all_colors(n = nlevels(temp$family) - 1), "grey"), levels(temp$family))
    
    # Plot
    temp %>%
      ggplot(aes(x = value, y = type, fill = family)) +
      geom_density_ridges(alpha = 0.2) +
      facet_wrap(~ trait, ncol = 1, scale = "free") +
      scale_fill_manual(values = family_color, name = "Family", guide = FALSE) +
      theme_acs() +
      theme(axis.title = element_blank())
    
  })

# Cowplot
g_density_plot <- plot_grid(plotlist = g_vp_family_density, ncol = 1, align = "hv")
ggsave(filename = "vp_pheno_mean_density.jpg", plot = g_density_plot, path = fig_dir,
       height = 7, width = 4, dpi = 1000)









