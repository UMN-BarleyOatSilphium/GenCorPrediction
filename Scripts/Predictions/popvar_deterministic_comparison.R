## GenCorPrediction
## 
## Compare predictions using PopVar with predictions based on the deterministic formula
## 
##

# Run on a local machine
repo_dir <- getwd()
source(file.path(repo_dir, "source.R"))
# Additional libraries
library(pbsim)
library(cowplot)

load(file.path(gdrive_dir, "BarleyLab/Projects/SideProjects/Resources/s2_cap_simulation_data.RData"))


# Remove monomorphic SNPs
# Filter for breeding programs relevant to my data
s2_cap_genos <- s2_cap_genos[str_detect(string = row.names(s2_cap_genos), pattern = "AB|BA|WA|N2|MT"),]
s2_cap_genos <- s2_cap_genos[,!colMeans(s2_cap_genos) %in% c(0, 2)]
s2_snp_info <- subset(s2_snp_info, rs %in% colnames(s2_cap_genos))



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





# # Pass to PopVar
# Convert genotypes into something useable for PopVar
geno_use1 <- as.data.frame(cbind( c("", row.names(geno_use)), rbind(colnames(geno_use), geno_use)) )

system.time({pred_out_pv <- PopVar::pop.predict(G.in = geno_use1, y.in = pheno_use, map.in = map_use,
                                                crossing.table = crossing_block, tail.p = i_sp, nInd = sim_pop_size,
                                                min.maf = 0, mkr.cutoff = 1, entry.cutoff = 1, remove.dups = FALSE,
                                                nSim = 25, nCV.iter = 1, models = "rrBLUP", impute = "pass")})


## Predicting 100 crosses with 25 iterations of 150 lines each took 200 seconds

                          
pred_out_tidy <- tidy.popvar(pred_out_pv) %>%
  select(trait, parent1 = Par1, parent2 = Par2, pred.mu, pred.varG, mu.sp_high, contains("cor")) %>%
  left_join(., filter(select(pred_out_tidy, parent1, parent2, pred.corG = `cor_w/_trait2`), !is.na(pred.corG))) %>%
  select(-contains("cor_"))
  


## Correlate predictions of genetic variance and correlation
pred_combined <- pred_out %>%
  left_join(., pred_out_tidy)


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
  xlab("Predicted cross variance\n(deterministic)") +
  ylab("Predicted cross variance\n(PopVar)") +
  theme_presentation2()

g_pred_corG <- pred_combined %>% 
  filter(trait == "trait1") %>%
  ggplot(aes(x = pred_corG, y = pred.corG)) +
  geom_abline(slope = 1, intercept = 0) +
  geom_point() +
  xlab("Predicted cross genetic\n(deterministic)") +
  ylab("Predicted cross genetic\n(PopVar)") +
  theme_presentation2()


## Combine plots
g_combine <- plot_grid(g_pred_mu, g_pred_varG, g_pred_corG, ncol = 1, labels = LETTERS[1:3])
ggsave(filename = "popvar_equation_compare.jpg", plot = g_combine, path = fig_dir, height = 10, width = 5, dpi = 1000)

