## PopVarVal simulations of prediction accuracy for the genetic
## correlation between two traits
## 
## # Demonstration simulations
## 
## 
## Author: Jeff Neyhart
## Last modified: October 12, 2018
## 

# Run on a local machine
repo_dir <- getwd()
source(file.path(repo_dir, "source.R"))
# Additional libraries
library(pbsim)
library(PopVar)

# Simulate a genome with 5 QTL, where 4 will be considered markers
gen_map <- list(`1` = c("M1" = 0, "M2" = 5, "Q" = 7.5, "M3" = 10, "M4" = 15))
genome <- sim_genome(map = gen_map)

genos <- replicate(length(names(gen_map[[1]])), c(2, 0)) %>%
  `dimnames<-`(list(c("founder1", "founder2"), names(gen_map[[1]])))


ped <- sim_pedigree(n.ind = 100)
cb <- sim_crossing_block(indnames(parents), , 1)

sim_out_rr <- replicate(n = 500, expr = {

  qtl.model <- replicate(n = 2, cbind(1, gen_map[[1]][-3], sample(rep(c(0.25, -0.25), each = 2)), 0), simplify = F)

  # Create a genome with the markers
  genome_markers <- sim_multi_gen_model(genome, qtl.model = qtl.model, prob.corr = cbind(0, 1))
  genome_markers$gen_model <- genome_markers$gen_model %>% map(~mutate(., qtl_name = names(gen_map[[1]])[-3]))
  
  # Simulate founders
  parents <- create_pop(genome = genome_markers, geno = genos)
  
  # Simulate a family and calculate the correlation
  cor(sim_family(genome = genome_markers, pedigree = ped, founder.pop = parents)$geno_val[,-1])[1,2]
  
})


sim_out_vc <- replicate(n = 500, expr = {
  
  qtl.model <- replicate(n = 2, cbind(1, gen_map[[1]][-3], c(0, sample(c(0.5, -0.5)), 0), 0), simplify = F)
  
  # Create a genome with the markers
  genome_markers <- sim_multi_gen_model(genome, qtl.model = qtl.model, prob.corr = cbind(0, 1))
  genome_markers$gen_model <- genome_markers$gen_model %>% map(~mutate(., qtl_name = names(gen_map[[1]])[-3]))
  
  # Simulate founders
  parents <- create_pop(genome = genome_markers, geno = genos)

  # Simulate a family and calculate the correlation
  cor(sim_family(genome = genome_markers, pedigree = ped, founder.pop = parents)$geno_val[,-1])[1,2]
  
})


## Plot
sim_out_df <- data.frame(AllMarkers = sim_out_rr, VariableSelection = sim_out_vc) %>%
  gather(method, correlation)

g_example_cor <- sim_out_df %>% 
  ggplot(aes(x = correlation, fill = method)) + 
  geom_histogram(position = "dodge", bins = 10) +
  ylab("Count") +
  xlab("Predicted correlation") +
  scale_fill_discrete(name = NULL) + 
  theme_poster() +
  theme(legend.position = c(0.5, 0.5))

ggsave(filename = "example_marker_selection_correlation.jpg", plot = g_example_cor, path = fig_dir, width = 4, height = 5, dpi = 1000)





