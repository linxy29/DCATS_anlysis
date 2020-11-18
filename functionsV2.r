# functions used start from 'simulation_plus5'(included)

if(FALSE){
  library(splatter)
  library(Seurat)
  library(speckle)
  library(DCATS)
  library(ggplot2)
  library(tidyverse)
  library(diffcyt)
}

## function1: cell selection
cellSelect = function(sim_mat, origLabels, ppC1, ppC2, sample_size){
  set.seed(123)
  cell_numC1 = (sample_size*ppC1) %>% ceiling()
  cell_numC2 = (sample_size*ppC2) %>% ceiling()
  ## condition 1
  group1_idx = c(1:length(origLabels))[origLabels == "Group1"] %>% sample(cell_numC1[1])
  group2_idx = c(1:length(origLabels))[origLabels == "Group2"] %>% sample(cell_numC1[2])
  group3_idx = c(1:length(origLabels))[origLabels == "Group3"] %>% sample(cell_numC1[3])
  indexC1 = c(group1_idx, group2_idx, group3_idx) %>% sort()
  summary(indexC1)
  simC1_mat = sim_mat[,indexC1]
  origLabelsC1 = origLabels[indexC1]
  ## condition 2
  group1_idx = c(1:length(origLabels))[origLabels == "Group1"] %>% sample(cell_numC2[1])
  group2_idx = c(1:length(origLabels))[origLabels == "Group2"] %>% sample(cell_numC2[2])
  group3_idx = c(1:length(origLabels))[origLabels == "Group3"] %>% sample(cell_numC2[3])
  indexC2 = c(group1_idx, group2_idx, group3_idx) %>% sort()
  summary(indexC2)
  simC2_mat = sim_mat[,indexC2]
  origLabelsC2 = origLabels[indexC2]
  return(cell_sltL = list(simC1_mat = simC1_mat, simC2_mat = simC2_mat, origLabelsC1 = origLabelsC1, origLabelsC2 = origLabelsC2))
}

## function2: Seurat process(might need few minutes)
runSeurat = function(cell_sltL, setresolu){
  # set input
  simC1_mat = cell_sltL$simC1_mat
  simC2_mat = cell_sltL$simC2_mat
  
  # pre-process
  seuratC1 <- CreateSeuratObject(counts = simC1_mat, project="Splatter")
  seuratC1 <- AddMetaData(object = seuratC1, metadata = rep("Cond1", dim(simC1_mat)[2]), col.name = 'condition')
  seuratC2 <- CreateSeuratObject(counts = simC2_mat, project="Splatter")
  seuratC2 <- AddMetaData(object = seuratC2, metadata = rep("Cond2", dim(simC2_mat)[2]), col.name = 'condition')
  
  listSamples = list(cond2 = seuratC1, cond2 = seuratC2)
  # log-normalization and identify variable features
  for (i in 1:length(listSamples)) {
    listSamples[[i]] <- NormalizeData(listSamples[[i]], verbose = FALSE)
    listSamples[[i]] <- FindVariableFeatures(listSamples[[i]], selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  }
  
  # integrate all batches
  anchors <- FindIntegrationAnchors(object.list = listSamples, dims = 1:30, verbose = FALSE)
  integratedSamples <- IntegrateData(anchorset = anchors, dims = 1:30, verbose = FALSE)
  DefaultAssay(integratedSamples) <- "integrated"
  
  # clustering
  integratedSamples <- ScaleData(integratedSamples, verbose = FALSE)
  integratedSamples <- RunPCA(integratedSamples, npcs = 30, verbose = FALSE)
  integratedSamples <- FindNeighbors(integratedSamples, dims = 1:30, verbose = FALSE)
  integratedSamples <- FindClusters(integratedSamples, resolution = setresolu, algorithm=2, verbose = FALSE)
  integratedSamples <- RunUMAP(integratedSamples, reduction = "pca", dims = 1:30, verbose = FALSE)
  
  # change labels to A, B, C
  integratedSamples@active.ident = integratedSamples@active.ident %>% 
    #plyr::mapvalues(from = c("0", "1", "2", "3", "4", "5"), to = c("A", "B", "C", "D", "E", "F"))
    plyr::mapvalues(from = c(0:11), to = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L"))
  
  return(integratedSamples)
}

## function3: rewrite betabinLRT to add bias correction
betabinLRT_rw <- function(counts1, counts2, pseudo_count=NULL, binom_only=FALSE, bias_corr=TRUE, similarity_mat=NULL) {
  ## Check counts1 and counts2 shape
  if (length(counts1) == 1 || is.null(dim(counts1)) ||
      length(dim(counts1)) < 2) {
    counts1 = matrix(counts1, nrow=1)
  }
  if (length(counts2) == 1 || is.null(dim(counts2)) ||
      length(dim(counts2)) < 2) {
    counts2 = matrix(counts2, nrow=1)
  }
  
  ## add base pseudo count if zero cells for all replicate
  if (is.null(pseudo_count)) {
    if (any(colMeans(counts1) == 0) || any(colMeans(counts2) == 0) ) {
      print(paste("Empty cell type exists in at least one conidtion;",
                  "adding replicate & condition specific pseudo count:"))
      print(0.01 * rowMeans(counts1))
      print(0.01 * rowMeans(counts2))
      
      # counts1 = counts1 + 0.01 * rowMeans(counts1)
      # counts2 = counts2 + 0.01 * rowMeans(counts2)
      counts1 = counts1 + 1
      counts2 = counts2 + 1
    }
  } else {
    counts1 = counts1 + pseudo_count
    counts2 = counts2 + pseudo_count
  }
  prop1 <- counts1 / rowSums(counts1)
  prop2 <- counts2 / rowSums(counts2)
  
  ## using estimated the latent cell counts
  counts1_latent = counts1
  counts2_latent = counts2
  if (bias_corr) {
    for (i in seq_len(nrow(counts1))) {
      counts1_latent[i, ] <- sum(counts1[i, ]) *
        multinom_EM(counts1[i, ], similarity_mat, verbose = FALSE)$mu
      counts1_latent[i, ] = ceiling(counts1_latent[i, ])
    }
    for (i in seq_len(nrow(counts2))) {
      counts2_latent[i, ] <- sum(counts2[i, ]) *
        multinom_EM(counts2[i, ], similarity_mat, verbose = FALSE)$mu
      counts2_latent[i, ] = ceiling(counts2_latent[i, ])
    }
  }
  
  ## number of cell types
  K <- ncol(counts1_latent)
  
  ## Binomial regression for each sampling
  intercept_val <- rep(NA, K)
  intercept_err <- rep(NA, K)
  
  coeffs_val <- rep(NA, K)
  coeffs_err <- rep(NA, K)
  LR_vals <- rep(NA, K)
  LRT_pvals <- rep(NA, K)
  
  totals <- c(rowSums(counts1_latent), rowSums(counts2_latent))
  labels <- c(rep(1, nrow(counts1_latent)), rep(0, nrow(counts2_latent)))
  
  if (binom_only) {
    for (i in seq_len(K)) {
      n1 <- c(counts1_latent[, i], counts2_latent[, i])
      df_tmp <- data.frame(n1 = n1, n2 = totals - n1, label=labels)
      model0 <- glm(cbind(n1, n2) ~ 1, family = binomial(),
                    data = df_tmp)
      model1 <- glm(cbind(n1, n2) ~ label + 1, family = binomial(),
                    data = df_tmp)
      
      intercept_val[i] <- summary(model1)$coefficients[1, 1]
      intercept_err[i] <- summary(model1)$coefficients[1, 2]
      
      coeffs_val[i] <- summary(model1)$coefficients[2, 1]
      coeffs_err[i] <- summary(model1)$coefficients[2, 2]
      
      LR_vals[i] <- model0$deviance - model1$deviance
      LRT_pvals[i] <-  pchisq(LR_vals[i], df=1, lower.tail = FALSE,
                              log.p = FALSE)
    }
  } else{
    for (i in seq_len(K)) {
      n1 <- c(counts1_latent[, i], counts2_latent[, i])
      df_tmp <- data.frame(n1 = n1, n2 = totals - n1, label=labels)
      
      fm0 <- aod::betabin(cbind(n1, n2) ~ 1, ~ 1, data = df_tmp)
      fm1 <- aod::betabin(cbind(n1, n2) ~ label, ~ 1, data = df_tmp)
      
      intercept_val[i] <- fm1@param[1]      # summary(fm1)@Coef[1, 1]
      intercept_err[i] <- fm1@varparam[1, 1]
      
      coeffs_val[i] <- fm1@param[2]
      coeffs_err[i] <-fm1@varparam[2, 2]
      
      LR_vals[i] <- fm0@dev - fm1@dev
      LRT_pvals[i] <-  pchisq(LR_vals[i], df=1, lower.tail = FALSE,
                              log.p = FALSE)
    }
  }
  
  data.frame(
    "prop1_mean" = colMeans(prop1),
    "prop1_std"  = matrixStats::colSds(prop1),
    "prop2_mean" = colMeans(prop2),
    "prop2_std"  = matrixStats::colSds(prop2),
    "coeff_mean" = coeffs_val,
    "coeff_std"  = sqrt(coeffs_err),
    "intecept_mean" = intercept_val,
    "intecept_std"  = sqrt(intercept_err),
    "pvals" = LRT_pvals,
    "LR" = LR_vals,
    row.names = colnames(counts1_latent)
  )
}