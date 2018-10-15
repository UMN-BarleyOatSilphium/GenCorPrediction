## GenCorPrediction 
## 
## Simulation checking
## 
## Run some checks on the genetic architecture of simulated genomes
## 
## 
## Author: Jeff Neyhart
## Last modified: October 15, 2018
## 





# Run on a local machine
repo_dir <- getwd()
source(file.path(repo_dir, "source.R"))
# Additional libraries
library(pbsim)

# Load the two-row simulation genotypes
load(file.path(gdrive_dir, "BarleyLab/Projects/SideProjects/Resources/s2_cap_simulation_data.RData"))


# Filter for breeding programs relevant to my data
s2_cap_genos <- s2_cap_genos[str_detect(string = row.names(s2_cap_genos), pattern = "AB|BA|WA|N2|MT"),]
s2_cap_genos <- s2_cap_genos[,!colMeans(s2_cap_genos) %in% c(0, 2)]
s2_snp_info <- subset(s2_snp_info, rs %in% colnames(s2_cap_genos))




##### Recurrent selection experiment

## Fixed parameters
tp_size <- 600

L <- 100
n_iter <- 100


## Outline the parameters to perturb
probcor_list <- list(cbind(5, 1), rbind(c(25, 0), c(35, 1)) )


map_sim <- s2_snp_info %>%
  split(.$chrom) %>%
  map(~setNames(.$cM_pos, .$rs)) %>%
  map(~structure(., class = "A"))  %>%
  structure(., class = "map") %>%
  # Jitter
  qtl::jittermap(.) %>%
  `names<-`(., seq_along(.))

# Create the base genome - this is fixed for all simulations
genome <- sim_genome(map = map_sim)


i <- 2
probcor <- probcor_list[[i]]
# probcor <- probcor[2,,drop = F]

# Simulate QTL
qtl_model <- replicate(n = 2, matrix(NA, ncol = 4, nrow = L), simplify = FALSE)



dist_mean <- replicate(100, {
  genome1 <- sim_multi_gen_model(genome = genome, qtl.model = qtl_model, add.dist = "geometric", max.qtl = L,
                                 corr = 0.5, prob.corr = probcor)
  
  
  
  ## What is the range of genetic distances between trait1 QTL and trait2 QTL
  qtl <- genome1$gen_model[[2]] %>% 
    split(., .$chr) %>%
    map(~list(.$qtl_name, .$qtl1_pair))
  
  qtl_dist <- map2(.x = map_sim, .y = qtl, ~{
    qtl1 <- .y[[1]]
    qtl2 <- .y[[2]]
    
    distance <- matrix(0, nrow = length(qtl1), ncol = length(qtl2), dimnames = list(qtl1, qtl2))
    for (i in seq(nrow(distance))) {
      for (j in seq(nrow(distance))) {
        distance[i,j] <- abs(diff(.x[c(qtl1[i], qtl2[j])]))
      }
    }
    
    distance
    
  })
    
  qtl_dist %>% 
    map(diag) %>% 
    unlist() %>% 
    mean()
  
})


  




























