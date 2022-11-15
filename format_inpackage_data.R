## This file is used to generate the data contains in the `DCATS` package

library(Seurat)
library(tidyverse)
library(DCATS)

## the simulation.Rdata
cluster_num = 4  # numbers of clusters
concentration = 100 # indicate how simulated proportion far away from true, the larger the closer
setresolu = 0.6
rep1 = 2
rep2 = 1
simulation_times = 30
sample_size1 = 1000
sample_size2 = 1000
cell_pool = "splatter"  # cell pool used, can be selected from ("SPAR", "splatter", "realWorld")
more_negative = ""
probC1 = c(0.3, 0.3, 0.2, 0.2)  # True proportions
probC2 = c(0.2, 0.3, 0.2, 0.3)
truthRes = c("N", "P", "P", "N")
source("functionsV2.r")
source("glm.R")     # for scDC
options(future.globals.maxSize = 15000 * 1024^2)

cluster_num = 8
if (cell_pool == "splatter"){
  #load(str_c("D:/#Data/DCATS/cells_pool", as.character(cluster_num), "_splatter.RData"))
  load(str_c("/storage/holab/linxy/DCATS/cells_pool", as.character(cluster_num), "_splatter.RData"))
} else if (cell_pool == "SPAR") {
  load(str_c("D:/Data/DCATS/cells_pool", as.character(cluster_num), "_SPARSim.RData"))
  
} else if (cell_pool == "realWorld") {
  load(str_c("D:/Data/DCATS/cells_pool", as.character(cluster_num), "_realWorld.RData"))
}

set.seed(123)
simulation = simulator_noInt(totals1 = runif(rep1, sample_size1, sample_size2), totals2 = runif(rep2, sample_size1, sample_size2), probC1, probC2, setresolu, sim_mat)
## sim_mat will be loaded from cell_pool

true_count = rbind(cond1 = simulation$numb_cond1, cond2 = simulation$numb_cond2)
## data convert
numb_cond1 = simulation$dfRes %>% 
  filter(condition == "Cond1") %>%
  group_by(clusterRes, batch) %>% 
  summarise(n = n()) %>% 
  pivot_wider(names_from = "clusterRes", values_from = "n") %>% 
  column_to_rownames(var = "batch")
numb_cond2 = simulation$dfRes %>% 
  filter(condition == "Cond2") %>%
  group_by(clusterRes, batch) %>% 
  summarise(n = n()) %>% 
  pivot_wider(names_from = "clusterRes", values_from = "n") %>% 
  column_to_rownames(var = "batch")

graphs = simulation$integratedSamples@graphs$RNA_snn
labels = Idents(simulation$integratedSamples)
knn_mat = knn_simMat(graphs, labels)
DimPlot(simulation$integratedSamples)

simulation = list(numb_cond1 = numb_cond1, numb_cond2 = numb_cond2, knn_mat = knn_mat, knnGraphs = graphs, labels = labels)
save(simulation, file = "/storage/holab/linxy/DCATS/package_data/current_version/simulation.RData")