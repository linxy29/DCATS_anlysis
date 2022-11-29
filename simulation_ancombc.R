## This file use ancombc in the simulation data

library(tidyverse)
library(ANCOMBC)
library(TreeSummarizedExperiment)

## load data
load("/Users/linxy29/Documents/Data/DCATS/simulation/current/replicates2&2_K8_con100_splatter3000&3000para.RData")

## ancombc test
count_info = rbind(oper[[1]]$seurat_count$cond1, oper[[1]]$seurat_count$cond2)
assay_data= count_info %>% 
col_date()

filenames = list.files("/Users/linxy29/Documents/Data/DCATS/simulation/current", full.names = TRUE)