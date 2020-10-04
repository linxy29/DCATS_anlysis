library(splatter)
library(Seurat)
library(speckle)
library(DCATS)
library(ggplot2)
library(tidyverse)

source("functions.r")

set.seed(123)
probNor = c(0.1, 0.2, 0.1, 0.2, 0.2, 0.2)
probMut = c(0.1, 0.2, 0.05, 0.1, 0.25, 0.3)
setresolu = 0.5
batch_size = 1000
de_prob = c(0.5,0.5,0.1, 0.1, 0.05, 0.05)

sim_list = simualtion(probNor, probMut, de_prob, batch_size)
integratedSamples = runSeurat(sim_list, batch_size, setresolu)

time = rep(NA,3)
dfRes = data.frame(clusterRes = integratedSamples@active.ident, batch = integratedSamples$batch, condition = integratedSamples$condition) %>% 
  tibble::rownames_to_column("cellID")

res_scDC <- scDC_noClustering(cellTypes = dfRes$clusterRes, dfRes$batch, calCI = TRUE, 
                              calCI_method = c("percentile", "BCa", "multinom"),
                              nboot = 50)

source("glm.R")
res_GLM <- fitGLM_fixed(res_scDC, c(rep("cond1",18),rep("cond2",18)), pairwise = FALSE)
print(res_GLM)