## This file use ancombc in the simulation data

library(tidyverse)
library(ANCOMBC)
library(TreeSummarizedExperiment)

## load data
source("functionsV2.r")

## simulation 2
filenames = list.files("/Users/linxy29/Documents/Data/DCATS/simulation/current", full.names = TRUE)
for (file in filenames[1]) {
  load(file)
  for (i in 1:length(oper)) {
    if (is.na(oper[[i]])[1]) next
    print(i)
    count_info = rbind(oper[[i]]$seurat_count$cond1, oper[[i]]$seurat_count$cond2)
    #count_info = rbind(oper[[i]]$seurat_count$cond1, oper[[i]]$seurat_count$cond2) %>%  mutate(condition = str_replace(condition, "Condition ", "cond"))
    assay_data = count_info %>% 
      #select(-condition) %>% 
      as.matrix() %>% t()
    col_data = data.frame(condition = rownames(count_info) %>% str_remove("s[1-9]"))
    t9start = Sys.time()
    ancombcRes = getANCOMBC(assay_data, col_data)
    oper[[i]]$time = c(oper[[i]]$time, Sys.time() - t9start)
    ## ancombc test with bias correction
    dcats_bc_knn = dcats_bc(assay_data %>% t(), oper[[i]]$knn_matrix)
    dcats_bc_true = dcats_bc(assay_data %>% t(), oper[[i]]$true_matrix)
    dcats_bc_svm = dcats_bc(assay_data %>% t(), oper[[i]]$svm_matrix)
    ancombc_knn = getANCOMBC(dcats_bc_knn %>% t(), col_data)
    ancombc_true = getANCOMBC(dcats_bc_true %>% t(), col_data)
    ancombc_svm = getANCOMBC(dcats_bc_svm %>% t(), col_data)
    ## merge results
    ancombcRes = ancombcRes %>% 
      dplyr::rename(cluster = taxon) %>% 
      dplyr::rename(ancombc_pvals = q_conditioncond2) %>% 
      mutate(ancombc_knn_pvals = ancombc_knn$q_conditioncond2,
             ancombc_true_pvals = ancombc_true$q_conditioncond2,
             ancombc_svm_pvals = ancombc_svm$q_conditioncond2)
    oper[[i]]$clusterDF = oper[[i]]$clusterDF %>% 
      merge(ancombcRes)
  }
  save(oper, file = file)
}

## simulation 3
file = filenames[2]
load(file)
for (i in 1:length(oper)) {
  if (is.na(oper[[i]])[1]) next
  print(i)
  count_info = oper[[i]]$seurat_count
  #count_info = rbind(oper[[i]]$seurat_count$cond1, oper[[i]]$seurat_count$cond2) %>%  mutate(condition = str_replace(condition, "Condition ", "cond"))
  assay_data = count_info %>% 
    #select(-condition) %>% 
    as.matrix() %>% t()
  col_data = design_mat
  t9start = Sys.time()
  ancombc_nbc = getANCOMBC(assay_data, col_data)
  oper[[i]]$time = c(oper[[i]]$time, Sys.time() - t9start)
  ## ancombc test with bias correction
  dcats_bc_knn = dcats_bc(assay_data %>% t(), oper[[i]]$knn_matrix)
  dcats_bc_true = dcats_bc(assay_data %>% t(), oper[[i]]$true_matrix)
  dcats_bc_svm = dcats_bc(assay_data %>% t(), oper[[i]]$svm_matrix)
  ancombc_knn = getANCOMBC(dcats_bc_knn %>% t(), col_data)
  ancombc_true = getANCOMBC(dcats_bc_true %>% t(), col_data)
  ancombc_svm = getANCOMBC(dcats_bc_svm %>% t(), col_data)
  ## merge results
  ancombcRes = ancombc_nbc %>% 
    dplyr::rename(cluster = taxon) %>% 
    dplyr::rename(ancombc_pvals = q_conditioncond2) %>% 
    dplyr::select(cluster, ancombc_pvals) %>% 
    mutate(ancombc_knn_pvals = ancombc_knn$q_conditioncond2,
           ancombc_true_pvals = ancombc_true$q_conditioncond2,
           ancombc_svm_pvals = ancombc_svm$q_conditioncond2)
  oper[[i]]$conditionDF = oper[[i]]$conditionDF %>% 
    merge(ancombcRes)
  ancombcRes = ancombc_nbc %>% 
    dplyr::rename(cluster = taxon) %>% 
    dplyr::rename(ancombc_pvals = q_gendermale) %>% 
    dplyr::select(cluster, ancombc_pvals) %>% 
    mutate(ancombc_knn_pvals = ancombc_knn$q_gendermale,
           ancombc_true_pvals = ancombc_true$q_gendermale,
           ancombc_svm_pvals = ancombc_svm$q_gendermale)
  oper[[i]]$genderDF = oper[[i]]$genderDF %>% 
    merge(ancombcRes)
}
save(oper, file = file)



