library(splatter)
library(Seurat)
library(speckle)
library(DCATS)
library(ggplot2)
library(tidyverse)
library(diffcyt)

## HCC data
HCC_seurat <- readRDS("/storage/holab/hcc_data/RDS_files/SCT/HCC.final.SCT_pc10_res0.1.rds")
count_HCC = HCC_seurat@assays$RNA@counts
params <- splatEstimate(as.matrix(count_HCC))
params
saveRDS(params, "./data/HCCparams.rds")



## mnnCT7
mnnCT7_seurat <- readRDS("/storage/holab/Nelson/mnnCT7.rds")
count_mnnCT7 = mnnCT7_seurat@assays$RNA@counts
params <- splatEstimate(as.matrix(count_mnnCT7))
params
saveRDS(params, "./data/mnnCT7params.rds")

ovarian_raw = read.csv("./data/real_word/ovarian_raw.csv.gz", )
params <- splatEstimate(as.matrix(ovarian_raw[,-1]))
params
saveRDS(params, "./data/ovarian_params.rds")

