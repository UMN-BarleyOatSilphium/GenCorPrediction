## GenCorPrediction
## 
## Compare predictions using PopVar with predictions based on the deterministic formula
## 
## Author: Jeff Neyhart
## Last modified: March 27, 2018
## 
##

# Run on a local machine
repo_dir <- getwd()
source(file.path(repo_dir, "source.R"))
# Additional libraries
library(pbsim)
library(pbsimData)
library(cowplot)



## Fixed parameters
sim_pop_size <- 150
n_iter <- 50
n_env <- 3
n_rep <- 1
n_crosses <- 100
k_sp <- 1.76
i_sp <- 0.2


map_sim <- s2_snp_info %>%
  split(.$chrom) %>%
  map(~setNames(.$cM_pos, .$rs)) %>%
  map(~structure(., class = "A"))  %>%
  structure(., class = "map") %>%
  # Jitter
  qtl::jittermap(.) %>%
  `names<-`(., seq_along(.))

# Create the base genome and the pedigree - these are fixed for all simulations
genome <- sim_genome(map = map_sim)
# Pedigree to accompany the crosses
ped <- sim_pedigree(n.ind = sim_pop_size, n.selfgen = Inf)

probcor_list <- data_frame(arch = c("pleio", "close_link", "loose_link"),
                           input = list(cbind(0, 1), cbind(5, 1), rbind(c(25, 0), c(35, 1)) ))


trait1_h2 <- 0.6
trait2_h2 <- 0.3
L <- 100
maxL <- 100
tp_size <- 600
gencor <- 0.5
probcor <- probcor_list[1,]$input[[1]]

maxL <- ifelse(probcor[1] == 0, maxL + (L / 2), maxL)


# Simulate QTL
qtl_model <- replicate(n = 2, matrix(NA, ncol = 4, nrow = L), simplify = FALSE)
genome1 <- sim_multi_gen_model(genome = genome, qtl.model = qtl_model, add.dist = "geometric", max.qtl = maxL,
                               corr = gencor, prob.corr = probcor)

## Adjust the genetic architecture - only if pleiotropy is not present
if (probcor[1] != 0) {
  genome1 <- adj_multi_gen_model(genome = genome1, geno = s2_cap_genos, gencor = gencor)
}




# Create the TP by random selection
tp1 <- create_pop(genome = genome1, geno = s2_cap_genos[sort(sample(nrow(s2_cap_genos), size = tp_size)), ]) %>% 
  # Phenotype the base population
  sim_phenoval(pop = ., h2 = c(trait1_h2, trait2_h2), n.env = n_env, n.rep = n_rep) %>%
  # Generate marker effects
  pred_mar_eff(genome = genome1, training.pop = ., method = "RRBLUP")




# Randomly create crosses from the cycle1 individuals
crossing_block <- sim_crossing_block(parents = indnames(tp1), n.crosses = n_crosses)




## Predict genetic variance and correlation
system.time({pred_out <- pred_genvar(genome = genome1, pedigree = ped, training.pop = tp1, founder.pop = tp1,
                                     crossing.block = crossing_block, method = "RRBLUP") %>%
  mutate(pred_musp = pred_mu + (k_sp * sqrt(pred_varG))) })

## Predicting 100 crosses took 1.5 seconds




## Check the predictions of correlation versus PopVar
# Convert items for PopVar
# Pheno
pheno_use <- tp1$pheno_val$pheno_mean
geno_use <- genotype(genome1, tp1)
# map_use <- genome1$gen_model$trait1 %>%
#   select(qtl_name, chr, pos)
map_use <- map_to_popvar(genome1) %>% 
  filter(marker %in% colnames(geno_use))

mar_eff_mat <- tp1$mar_eff %>% column_to_rownames("marker") %>% as.matrix()
mar_eff_mat1 <- mixed.solve(y = tp1$pheno_val$pheno_mean$trait1, Z = geno_use)$u
mar_eff_mat2 <- mixed.solve(y = tp1$pheno_val$pheno_mean$trait2, Z = geno_use)$u

## Predict genotypic values
pgv <- as.data.frame(geno_use %*% cbind(trait1 = mar_eff_mat1, trait2 = mar_eff_mat2)) %>%
  rownames_to_column("line_name") %>%
  gather(trait, pgv, -line_name)

## Predict cross means
cross_means <- left_join(crossing_block, pgv, by = c("parent1" = "line_name")) %>% 
  left_join(., pgv, by = c("parent2" = "line_name", "trait")) %>%
  group_by(parent1, parent2, trait) %>%
  mutate(cross_mean = (pgv.x + pgv.y) / 2)

## Correlate cross means with predictions from the deterministic formula
full_join(pred_out, cross_means) %>% group_by(trait) %>% summarize(cor = cor(pred_mu, cross_mean))
 
# trait    cor
# 1 trait1     1
# 2 trait2     1



# # Pass to PopVar
# Convert genotypes into something useable for PopVar
geno_use1 <- as.data.frame(cbind( c("", row.names(geno_use)), rbind(colnames(geno_use), geno_use)) )

system.time({pred_out_pv <- PopVar::pop.predict(G.in = geno_use1, y.in = pheno_use, map.in = map_use,
                                                crossing.table = crossing_block, tail.p = i_sp, nInd = 1000,
                                                min.maf = 0, mkr.cutoff = 1, entry.cutoff = 1, remove.dups = FALSE,
                                                nSim = 1, nCV.iter = 2, models = "rrBLUP", impute = "pass")})


## Predicting 100 crosses with 25 iterations of 150 lines each took 200 seconds
## Predicting 100 crosses with 1 iteration of 1000 lines each took ~35 seconds
## 


pred_out_tidy <- tidy.popvar(pred_out_pv) %>%
  select(trait, parent1 = Par1, parent2 = Par2, pred.mu, pred.varG, mu.sp_high, contains("cor")) %>%
  left_join(., filter(select(., parent1, parent2, pred.corG = `cor_w/_trait2`), !is.na(pred.corG))) %>%
  select(-contains("cor_"))
  


## Correlate predictions of genetic variance and correlation
pred_combined <- pred_out %>%
  left_join(., pred_out_tidy)

pred_combined_summ <- pred_combined %>% 
  group_by(trait) %>% 
  summarize(corMean = cor(pred_mu, pred.mu),
            corVarG = cor(pred_varG, pred.varG),
            corCorG = cor(pred_corG, pred.corG))

## They correlate perfectly, but there is bias
pred_combined %>% group_by(trait) %>% summarize(bias = (mean(pred.mu) - mean(pred_mu)) / mean(pred_mu)) 



## What is the predicted mean in the first two families of the crossing block?
family1 <- sim_family_cb(genome = genome1, pedigree = sim_pedigree(n.ind = 1000), founder.pop = tp1, crossing.block = crossing_block[1:2,]) %>%
  pred_geno_val(genome = genome1, training.pop = tp1, candidate.pop = .)


family1$pred_val %>% 
  mutate(family = str_extract(ind, "[0-9]{4}")) %>% 
  group_by(family) %>% 
  summarize_at(vars(trait1, trait2), mean)
head(pred_out)


g_pred_mu <- pred_combined %>% 
  ggplot(aes(x = pred_mu, y = pred.mu)) +
  geom_abline(slope = 1, intercept = 0) +
  geom_point() +
  facet_grid(~ trait) +
  xlab("Predicted cross mean\n(deterministic)") +
  ylab("Predicted cross mean\n(PopVar)") +
  theme_presentation2()

g_pred_varG <- pred_combined %>% 
  ggplot(aes(x = pred_varG, y = pred.varG)) +
  geom_abline(slope = 1, intercept = 0) +
  geom_point() +
  facet_grid(~ trait) +
  xlab("Predicted cross\ngenetic variance\n(deterministic)") +
  ylab("Predicted cross\ngenetic variance\n(PopVar)") +
  theme_presentation2()

g_pred_corG <- pred_combined %>% 
  filter(trait == "trait1") %>%
  ggplot(aes(x = pred_corG, y = pred.corG)) +
  geom_abline(slope = 1, intercept = 0) +
  geom_point() +
  xlab("Predicted cross\ngenetic correlation\n(deterministic)") +
  ylab("Predicted cross\ngenetic correlation\n(PopVar)") +
  theme_presentation2()


## Combine plots
g_combine <- plot_grid(g_pred_mu, g_pred_varG, g_pred_corG, ncol = 1, labels = LETTERS[1:3], rel_heights = c(0.9, 0.9, 1), align = "hv", axis = "tblr")
ggsave(filename = "popvar_equation_compare.jpg", plot = g_combine, path = fig_dir, height = 10, width = 5, dpi = 1000)

