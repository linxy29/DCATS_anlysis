### This is the file used for simulation of extreme cases
library(Seurat)
library(SeuratWrappers)
library(speckle)
library(DCATS)
library(ggplot2)
library(tidyverse)
library(diffcyt)
library(MCMCpack)
#library(scdney)
library(scDC)
library(tidymodels)  ## for DCATS
library(miloR)
library(SingleCellExperiment)

## simulation setting
cluster_num = 8  # numbers of clusters
setresolu = 0.6
rep1 = 3
rep2 = 3
simulation_times = 30
concentration = 30
sample_size1 = 3000
sample_size2 = 3000
cell_pool = "splatter"  # cell pool used, can be selected from ("SPAR", "splatter", "realWorld")
more_negative = ""
countC1 = rep(8, 8)
countC2 = c(12, 12, rep(8,6))
probC1 = countC1/sum(countC1)
probC2 = countC2/sum(countC2)
truthRes = c("T", "T", "F", "F", "F", "F", "F", "F")

## load data and function
load("/storage/holab/linxy/DCATS/cells_pool8_splatter.RData")
source("functionsV2.r")
source("glm.R")     # for scDC
options(future.globals.maxSize = 15000 * 1024^2)

## only simulate count
resDF = data.frame()
for (i in 1:50){
  simulation = simulator_base(rep(3000,3), rep(3000,3), probC1*70, probC2*70, diag(nrow = 8, ncol = 8))
  sim_design = data.frame(condition = c(rep("g1", rep1), rep("g2", rep2)))
  print(simulation$numb_cond1/rowSums(simulation$numb_cond1))
  print(simulation$numb_cond2/rowSums(simulation$numb_cond2))
  ref_order = detect_reference(rbind(simulation$numb_cond1, simulation$numb_cond2), sim_design)
  ref_res1 = dcats_new(rbind(simulation$numb_cond1, simulation$numb_cond2), sim_design, reference = ref_order[1:2])
  ref_res2 = dcats_new(rbind(simulation$numb_cond1, simulation$numb_cond2), sim_design, reference = ref_order[1])
  org_res = dcats_GLM(rbind(simulation$numb_cond1, simulation$numb_cond2), sim_design)
  sub_res = data.frame(truth = truthRes, ref_pval1 = as.vector(ref_res1$LRT_pvals), ref_pval2 = as.vector(ref_res2$LRT_pvals), org_pval = as.vector(org_res$LRT_pvals))
  resDF = rbind(resDF, sub_res)
}




