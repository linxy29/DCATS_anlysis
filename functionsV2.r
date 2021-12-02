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

#print("You are really a smart girl!!!!")

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

## function2: Seurat process(might need few minutes) for one duplicate
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
    plyr::mapvalues(from = c(0:15), to = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P"))
  
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
  
  
  for (i in seq_len(K)) {
    n1 <- c(counts1_latent[, i], counts2_latent[, i])
    df_tmp <- data.frame(n1 = n1, n2 = totals - n1, label=labels)
    if (binom_only) {
      model0 <- glm(cbind(n1, n2) ~ 1, family = binomial(), data = df_tmp)
      model1 <- glm(cbind(n1, n2) ~ label + 1, family = binomial(), data = df_tmp)
      
      intercept_val[i] <- summary(model1)$coefficients[1, 1]
      intercept_err[i] <- summary(model1)$coefficients[1, 2]
      
      coeffs_val[i] <- summary(model1)$coefficients[2, 1]
      coeffs_err[i] <- summary(model1)$coefficients[2, 2]
      
      LR_vals[i] <- model0$deviance - model1$deviance
      LRT_pvals[i] <-  pchisq(LR_vals[i], df=1, lower.tail = FALSE,
                              log.p = FALSE)
    } else{
      fm0 <- aod::betabin(cbind(n1, n2) ~ 1, ~ 1, data = df_tmp)
      fm1 <- aod::betabin(cbind(n1, n2) ~ label, ~ 1, data = df_tmp)
      
      if (any(is.na(fm1@varparam))) {
        print(fm1)
        print(df)
      } else {
        coeffs_err[i] <-fm1@varparam[2, 2]
      }
      
      intercept_val[i] <- fm1@param[1]      # summary(fm1)@Coef[1, 1]
      intercept_err[i] <- fm1@varparam[1, 1]
      
      coeffs_val[i] <- fm1@param[2]
      
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

## function4: cell selection for each duplicate
cell_slt_dup = function(cell_num, sim_mat, origLabels){
  K = length(cell_num)
  cellname = origLabels %>% as.factor() %>% levels()
  index = numeric()
  for (j in 1:K){
    group_idx = c(1:length(origLabels))[origLabels == cellname[j]] %>% 
      sample(cell_num[j])
    index = c(index, group_idx)
  }
  sub_sim_mat = sim_mat[,index]
  sub_origLabels = origLabels[index]
  return(list(sub_sim_mat = sub_sim_mat, sub_origLabels = sub_origLabels))
}

# function 5: Seurat process(might need few minutes) 
runSeurat_mul = function(slt_sim_matC1, slt_sim_matC2, slt_batchesC1, slt_batchesC2, setresolu){
  # set input
  simNor_mat = slt_sim_matC1
  simMut_mat = slt_sim_matC2
  batchNor = slt_batchesC1
  batchMut = slt_batchesC2
  
  # pre-process
  seuratNor <- CreateSeuratObject(counts = simNor_mat, project="Splatter")
  seuratNor <- AddMetaData(object = seuratNor, metadata = batchNor, col.name = 'batch')
  seuratNor <- AddMetaData(object = seuratNor, metadata = rep("Cond1", length(batchNor)), col.name = 'condition')
  seuratMut <- CreateSeuratObject(counts = simMut_mat, project="Splatter")
  seuratMut <- AddMetaData(object = seuratMut, metadata = batchMut, col.name = 'batch')
  seuratMut <- AddMetaData(object = seuratMut, metadata = rep("Cond2", length(batchMut)), col.name = 'condition')
  
  listNor = SplitObject(seuratNor, split.by = "batch")
  listMut = SplitObject(seuratMut, split.by = "batch")
  listSamples = c(listNor, listMut)
  
  # log-normalization and identify variable features
  for (i in 1:length(listSamples)) {
    listSamples[[i]] <- NormalizeData(listSamples[[i]], verbose = FALSE)
    listSamples[[i]] <- FindVariableFeatures(listSamples[[i]], selection.method = "vst", 
                                             nfeatures = 2000, verbose = FALSE)
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
    plyr::mapvalues(from = c(0:15), to = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P"))
  
  return(integratedSamples)
}

# function 6: fastMNN(useless)
fastMNN = function(slt_sim_matC1, slt_sim_matC2, slt_batchesC1, slt_batchesC2, setresolu){
  # set input
  simNor_mat = slt_sim_matC1
  simMut_mat = slt_sim_matC2
  batchNor = slt_batchesC1
  batchMut = slt_batchesC2
  
  # pre-process
  seuratObj <- CreateSeuratObject(counts = cbind(simNor_mat, simMut_mat), project="Splatter")
  seuratObj <- AddMetaData(object = seuratObj, metadata = c(batchNor, batchMut), col.name = 'batch')
  seuratObj <- AddMetaData(object = seuratObj, metadata = c(rep("Cond1", length(batchNor)), rep("Cond2", length(batchMut))), col.name = 'condition')
  seuratObj <- NormalizeData(seuratObj, verbose = FALSE)
  #seuratObj <- FindVariableFeatures(seuratObj, verbose = FALSE)
  integratedSamples <- RunFastMNN(object.list = SplitObject(seuratObj, split.by = "batch"), verbose = FALSE)
  integratedSamples <- FindVariableFeatures(integratedSamples, verbose = FALSE)
  integratedSamples <- ScaleData(integratedSamples, verbose = FALSE)
  integratedSamples <- RunPCA(integratedSamples, npcs = 30, verbose = FALSE)
  integratedSamples <- RunUMAP(integratedSamples, reduction = "mnn", dims = 1:30, verbose = FALSE)
  integratedSamples <- FindNeighbors(integratedSamples, reduction = "mnn", dims = 1:30, verbose = FALSE)
  integratedSamples <- FindClusters(integratedSamples, resolution = setresolu, verbose = FALSE)
  
  # change labels to A, B, C
  integratedSamples@active.ident = integratedSamples@active.ident %>% 
    #plyr::mapvalues(from = c("0", "1", "2", "3", "4", "5"), to = c("A", "B", "C", "D", "E", "F"))
    plyr::mapvalues(from = c(0:15), to = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P"))
  
  return(integratedSamples)
}

# function 7: simulator with fastMNN
simulator_fastMNN = function(totals1, totals2, probC1, probC2, setresolu, sim_mat){
  K = length(probC1)
  n_rep1 = length(totals1)
  n_rep2 = length(totals2)
  
  prop_cond1 = matrix(0, n_rep1, K)
  prop_cond2 = matrix(0, n_rep2, K)
  numb_cond1 <- matrix(0, n_rep1, K)
  numb_cond2 <- matrix(0, n_rep2, K)
  
  # cell selection
  slt_sim_matC1 = matrix(NA, ncol = 0, nrow = dim(sim_mat)[1])
  slt_origLabelsC1 = vector()
  slt_batchesC1 = character()
  for (rep in seq_len(n_rep1)) {
    prop_cond1[rep, ] <- MCMCpack::rdirichlet(1, probC1 * concentration)
    #prop_cond1[rep, ] = probC1
    numb_cond1[rep, ] <- rmultinom(1, totals1[rep], prop_cond1[rep, ])
    #numb_cond1[rep, ] = (totals1[rep] * prop_cond1[rep, ]) %>% ceiling()
    cell_slt = cell_slt_dup(numb_cond1[rep,], sim_mat, origLabels)
    slt_sim_matC1 = cbind(slt_sim_matC1, cell_slt$sub_sim_mat)
    slt_origLabelsC1 = c(slt_origLabelsC1, cell_slt$sub_origLabels)
    slt_batchesC1 = c(slt_batchesC1, rep(str_c("cond1s", as.character(rep)), length(cell_slt$sub_origLabels)))
  }
  
  slt_sim_matC2 = matrix(NA, ncol = 0, nrow = dim(sim_mat)[1])
  slt_origLabelsC2 = factor()
  slt_batchesC2 = character()
  for (rep in seq_len(n_rep2)) {
    prop_cond2[rep, ] <- MCMCpack::rdirichlet(1, probC2 * concentration)
    #prop_cond2[rep, ] = probC2
    numb_cond2[rep, ] <- rmultinom(1, totals2[rep], prop_cond2[rep, ])
    cell_slt = cell_slt_dup(numb_cond2[rep,], sim_mat, origLabels)
    slt_sim_matC2 = cbind(slt_sim_matC2, cell_slt$sub_sim_mat)
    slt_origLabelsC2 = c(slt_origLabelsC2, cell_slt$sub_origLabels)
    slt_batchesC2 = c(slt_batchesC2, rep(str_c("cond2s", as.character(rep)), length(cell_slt$sub_origLabels)))
  }
  
  print(numb_cond1)
  print(numb_cond2)
  
  ## seurat process
  ### pre-process
  seuratObj <- CreateSeuratObject(counts = cbind(slt_sim_matC1, slt_sim_matC2), project="Splatter")
  seuratObj <- AddMetaData(object = seuratObj, metadata = c(slt_batchesC1, slt_batchesC2), col.name = 'batch')
  seuratObj <- AddMetaData(object = seuratObj, metadata = c(rep("Cond1", length(slt_batchesC1)), rep("Cond2", length(slt_batchesC2))), col.name = 'condition')
  seuratObj <- NormalizeData(seuratObj, verbose = FALSE)
  integratedSamples <- RunFastMNN(object.list = SplitObject(seuratObj, split.by = "batch"), verbose = FALSE)
  integratedSamples <- FindVariableFeatures(integratedSamples, verbose = FALSE)
  integratedSamples <- ScaleData(integratedSamples, verbose = FALSE)
  integratedSamples <- RunPCA(integratedSamples, npcs = 30, verbose = FALSE)
  #integratedSamples <- RunUMAP(integratedSamples, reduction = "mnn", dims = 1:30, verbose = FALSE)
  #integratedSamples <- FindNeighbors(integratedSamples, reduction = "mnn", dims = 1:30, verbose = FALSE)
  integratedSamples <- RunUMAP(integratedSamples, dims = 1:30, verbose = FALSE)
  integratedSamples <- FindNeighbors(integratedSamples, dims = 1:30, verbose = FALSE)
  integratedSamples <- FindClusters(integratedSamples, resolution = setresolu, verbose = FALSE)
  
  ## decide resolution
  Kprep = integratedSamples@active.ident %>% as.factor() %>% summary() %>% length()
  str_c("Kprep: ", as.character(Kprep)) %>% print()
  str_c("setresolu: ", as.character(setresolu)) %>% print()
  while (Kprep != K & setresolu > 0.03) {
    if (Kprep > K){
      setresolu = setresolu - 0.03
      integratedSamples <- FindClusters(integratedSamples, resolution = setresolu, verbose = FALSE)
      Kprep = integratedSamples@active.ident %>% as.factor() %>% summary() %>% length()
      str_c("Kprep: ", as.character(Kprep)) %>% print()
      str_c("setresolu: ", as.character(setresolu)) %>% print()
    } else {
      setresolu = setresolu + 0.01
      integratedSamples <- FindClusters(integratedSamples, resolution = setresolu, verbose = FALSE)
      Kprep = integratedSamples@active.ident %>% as.factor() %>% summary() %>% length()
      str_c("Kprep: ", as.character(Kprep)) %>% print()
      str_c("setresolu: ", as.character(setresolu)) %>% print()
    }
  }
  ### change labels to A, B, C
  integratedSamples@active.ident = integratedSamples@active.ident %>% 
    plyr::mapvalues(from = c(0:15), to = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P"))
  
  # get count data
  dfRes = data.frame(clusterRes = integratedSamples@active.ident, batch = integratedSamples$batch, condition = integratedSamples$condition) %>% 
    tibble::rownames_to_column("cellID")
  if(Kprep == K){
    Res = list(integratedSamples = integratedSamples, dfRes = dfRes, trueLabels = c(slt_origLabelsC1, slt_origLabelsC2), numb_cond1 = numb_cond1, numb_cond2 = numb_cond2, prop_cond1 = prop_cond1, prop_cond2 = prop_cond2)
    return(Res)
  } else {
    return(NA)
  }

}

# function 8: calculate p-values of fisher's exact test
getFisher = function(counts1, counts2){
  ### Check counts1 and counts2 shape
  if (length(counts1) == 1 || is.null(dim(counts1)) || length(dim(counts1)) < 2) {
    counts1 = matrix(counts1, nrow=1)
  }
  if (length(counts2) == 1 || is.null(dim(counts2)) ||
      length(dim(counts2)) < 2) {
    counts2 = matrix(counts2, nrow=1)
  }
  ### set up
  totalCounts1 = colSums(counts1)
  totalCounts2 = colSums(counts2)
  count_mtrx = rbind(totalCounts1, totalCounts2)
  ### test 
  fisher_pvals = rep(NA,ncol(count_mtrx))
  for (i in 1:ncol(count_mtrx)) {
    tested_table = cbind(count_mtrx[, i], rowSums(count_mtrx)-count_mtrx[, i])
    fisher_pvals[i] = fisher.test(tested_table)$p.value
  }
  return(fisher_pvals)
}

# function 9: get p_values from diffcyt
getDiffcyt = function(numb_cond1, numb_cond2, dfRes){
  d_random <- function(batch_size, mean = 0, sd = 1, ncol = 20, cofactor = 5) {
    d <- sinh(matrix(rnorm(20*batch_size, mean, sd), ncol = ncol)) * cofactor
    colnames(d) <- paste0("marker", sprintf("%02d", 1:ncol))
    d
  }
  
  # Create random data (without differential signal)
  total_cond1 = apply(numb_cond1, 1, sum)
  total_cond2 = apply(numb_cond2, 1, sum)
  
  sample_numC1 = length(total_cond1)
  sample_numC2 = length(total_cond2)
  
  d_input <- vector(mode = "list", length = sample_numC1 + sample_numC2)
  for (i in 1:sample_numC1){
    d_input[[i]] = d_random(total_cond1[i])
  }
  for (i in 1:sample_numC2){
    d_input[[i+sample_numC1]] = d_random(total_cond2[i])
  }
  
  experiment_info <- data.frame(
    sample_id = dfRes$batch %>% unique() %>% as.factor(), 
    group_id = factor(c(rep("Cond1", sample_numC1), rep("Cond2", sample_numC2))), 
    stringsAsFactors = FALSE)
  marker_info <- data.frame(
    channel_name = paste0("channel", sprintf("%03d", 1:20)), 
    marker_name = paste0("marker", sprintf("%02d", 1:20)), 
    marker_class = factor(c(rep("type", 10), rep("state", 10)), 
                          levels = c("type", "state", "none")), 
    stringsAsFactors = FALSE)
  
  # Prepare data
  d_se <- prepareData(d_input, experiment_info, marker_info)
  # Transform data
  d_se <- transformData(d_se)
  # Generate clusters
  d_se <- generateClusters(d_se)
  # Replace with our simulated data info
  d_se@elementMetadata$sample_id = dfRes$batch
  d_se@elementMetadata$group_id = dfRes$condition
  d_se@elementMetadata$cluster_id = dfRes$clusterRes
  # Calculate counts
  d_counts <- calcCounts(d_se)
  # Create design matrix
  design <- createDesignMatrix(experiment_info, cols_design = "group_id")
  # Create contrast matrix
  contrast <- createContrast(c(0, 1))
  # Test for differential abundance (DA) of clusters
  res_DA <- testDA_edgeR(d_counts, design, contrast)
  #diffcytP = topTable(res_DA, format_vals = TRUE) %>% 
  diffcytP = res_DA@elementMetadata %>% 
    as.data.frame() %>% 
    dplyr::rename(cluster = cluster_id) %>% 
    dplyr::rename(diffcyt_pvals = p_adj) %>% 
    dplyr::select(cluster, diffcyt_pvals)
  
  return(diffcytP)
}

# function 10: for sending email
sendEmail = function(subject = "The job is done"){
  library(emayili)
  email <- envelope() %>%
    from("4sendEmail29@gmail.com") %>%
    to("xl2836@outlook.com") %>%
    subject(subject) %>%
    text("Yeah!!!!!")
  smtp <- server(host = "smtp.gmail.com", port = 465, username = "4sendEmail29@gmail.com", password = "dummy_acc29")
  smtp(email, verbose = FALSE)
}

# function 12: get similarity matrix from rf and svm
get_similarity_matRW = function(K, confuse_rate, method = "uniform", df){
  library(tidymodels)
  if (method == "uniform"){
    simil_mat = diag(K) * (1 - confuse_rate) + confuse_rate * (1 - diag(K))/(K - 1)
  } else if (method == "rf"| method == "svm"){
    cv <- vfold_cv(df, v = 5)
    recipe <- recipe(clusterRes ~ ., data = df)
    if (method == "rf"){
      model <- 
        rand_forest() %>%
        set_engine("ranger", importance = "impurity") %>%
        set_mode("classification") 
    } else {
      model <-
        svm_rbf() %>% 
        set_mode("classification") %>%
        set_engine("kernlab")
    }
    workflow <- workflow() %>%
      add_recipe(recipe) %>%
      add_model(model)
    predDF = data.frame()
    for (i in 1:5) {
      onefold_split = cv[[1]][[i]]
      fit <- workflow %>%
        last_fit(onefold_split)
      pred = fit %>% 
        collect_predictions() %>% 
        mutate(pred = .pred_class) %>% 
        dplyr::select(pred, clusterRes)
      predDF = rbind(predDF, pred)
    }
    conf.mat <- table(predDF$clusterRes, predDF$pred)
    simil_mat <- t(t(conf.mat)/apply(conf.mat,2,sum))
  }
  return(simil_mat)
}

# function 13: plot ROC plot
getROC = function(truth, pred){
  ## sensitivity and specificity
  library(ggplot2)
  thresholds = unique(c(0, pred))
  sensitivity = rep(NA, length(thresholds))
  specificity = rep(NA, length(thresholds))
  TP = rep(NA, length(thresholds))
  TN = rep(NA, length(thresholds))
  FP = rep(NA, length(thresholds))
  FN = rep(NA, length(thresholds))
  for (j in 1:length(thresholds)) {
    pred_res = ifelse(pred <= thresholds[j], "P", "N")
    TP[j] <- sum(pred_res=="P"&truth=="P")
    TN[j] <- sum(pred_res=="N"&truth=="N")
    FP[j] <- sum(pred_res=="P"&truth=="N")
    FN[j] <- sum(pred_res=="N"&truth=="P")
    sensitivity[j] = TP[j]/(TP[j]+FN[j])
    specificity[j] = TN[j]/(TN[j]+FP[j])
  }
  eval = data.frame(thresholds = thresholds, TP = TP, TN = TN, FP = FP, FN = FN)
  df = data.frame(sensitivity = sensitivity, specificity = specificity, x = 1-specificity) %>% 
    arrange(x)
  ## plot
  plot = df %>% 
    ggplot(aes(x = x, y = sensitivity)) +
    geom_line(orientation = "y") +
    geom_abline(slope=1, color="gray", linetype="dashed")
  ## auc
  auc = 0
  for (n in 2:length(thresholds)){
    subArea = (df$sensitivity[n-1]+df$sensitivity[n])*(df$x[n]-df$x[n-1])/2
    auc = auc + subArea
  }
  return(list(plot = plot, auc = auc, df = df, eval = eval))
}

# function 13: plot PRC plot
getPRC = function(truth, pred){
  ## precision and recall
  library(ggplot2)
  thresholds = unique(c(0, round(pred, 3)))
  precision = rep(NA, length(thresholds))
  recall = rep(NA, length(thresholds))
  TP = rep(NA, length(thresholds))
  TN = rep(NA, length(thresholds))
  FP = rep(NA, length(thresholds))
  FN = rep(NA, length(thresholds))
  for (j in 1:length(thresholds)) {
    pred_res = ifelse(pred <= thresholds[j], "P", "N")
    TP[j] <- sum(pred_res=="P"&truth=="P")
    TN[j] <- sum(pred_res=="N"&truth=="N")
    FP[j] <- sum(pred_res=="P"&truth=="N")
    FN[j] <- sum(pred_res=="N"&truth=="P")
    precision[j] = TP[j]/(TP[j]+FP[j])
    recall[j] = TP[j]/(TP[j]+FN[j])
  }
  eval = data.frame(thresholds = thresholds, TP = TP, TN = TN, FP = FP, FN = FN)
  df = data.frame(precision = precision, recall = recall) %>% 
    arrange(desc(precision))
  df[is.na(df)] <- 0
  df = df[rowSums(df)!=0,]
  ## plot
  plot = df %>% 
    ggplot(aes(x = recall, y = precision)) +
    geom_line(orientation = "x")
  ## prauc
  prauc = 0
  for (n in 2:nrow(df)){
    subArea = (df$precision[n-1]-df$precision[n])*(df$recall[n-1]+df$recall[n])/2
    prauc = prauc + subArea
  }
  prauc = prauc + df$precision[n] * 1 # add precidion till 0
  return(list(plot = plot, prauc = prauc, df = df, eval = eval))
}

# function 14: get Phi for each cluster
eachPhi <- function(count_mat, design_mat, similarity_mat=NULL, n_samples=NULL, pseudo_count=NULL,  base_model='NULL', fix_phi=NULL) {
  # Output matrices
  coeffs     <- matrix(NA, ncol(count_mat), ncol(design_mat))
  coeffs_err <- matrix(NA, ncol(count_mat), ncol(design_mat))
  LR_vals    <- matrix(NA, ncol(count_mat), ncol(design_mat))
  LRT_pvals  <- matrix(NA, ncol(count_mat), ncol(design_mat))
  pvals  <- matrix(NA, ncol(count_mat), ncol(design_mat))
  LRT_fdr    <- matrix(NA, ncol(count_mat), ncol(design_mat))
  fm0_phi <- matrix(NA, ncol(count_mat), ncol(design_mat))
  fm1_phi <- matrix(NA, ncol(count_mat), ncol(design_mat))
  
  # Check colnames
  if (is.null(colnames(count_mat)))
    colnames(count_mat) <- paste0('cell_type_', seq(ncol(count_mat)))
  if (is.null(colnames(design_mat)))
    colnames(design_mat) <- paste0('factor_', seq(ncol(design_mat)))
  
  # Add rownames and colnames
  rownames(fm0_phi) <- rownames(fm1_phi) <- rownames(LR_vals) <- rownames(LRT_pvals) <- rownames(pvals) <- rownames(LRT_fdr) <- rownames(coeffs) <- rownames(coeffs_err) <- colnames(count_mat)
  colnames(fm0_phi) <- colnames(fm1_phi) <- colnames(LR_vals) <- colnames(LRT_pvals) <- colnames(pvals) <- colnames(LRT_fdr) <-
    colnames(coeffs) <- colnames(coeffs_err) <- colnames(design_mat)
  
  
  ## using estimated the latent cell counts
  count_latent = count_mat
  if(!is.null(similarity_mat)) {
    for (i in seq_len(nrow(count_mat))) {
      count_latent[i, ] <- sum(count_mat[i, ]) *
        multinom_EM(count_mat[i, ], similarity_mat, verbose = FALSE)$mu
    }
  }
  
  K <- ncol(count_mat) ## number of cell types
  if (is.null(similarity_mat)) {
    n_samples <- 1
  }
  
  if (!is.null(n_samples) && !is.null(similarity_mat)) {
    count_use <- matrix(0, nrow(count_mat) * n_samples, K)
    for (i in seq_len(nrow(count_mat))) {
      idx <- seq((i - 1) * n_samples + 1, i * n_samples)
      for (j in seq_len(K)) {
        count_use[idx, ] <- (
          count_use[idx, ] + t(rmultinom(n_samples, count_latent[i, j], similarity_mat[j, ])))
      }
    }
  } else{
    count_use <- count_mat
    n_samples <- 1
  }
  
  # adding pseudo counts
  if (is.null(pseudo_count)) {
    if (any(colMeans(count_mat) == 0)) {
      print(paste("Empty cell type exists in at least one conidtion;",
                  "adding replicate & condition specific pseudo count:"))
      count_use <- count_use + 1
    }
  } else {
    count_use = count_use + pseudo_count
  }
  
  count_use = round(count_use)
  
  # Test each factor
  for (k in seq_len(ncol(design_mat))) {    ## for each factor
    sub_LR_val <- matrix(NA, n_samples, K)
    sub_coeffs_val <- matrix(NA, n_samples, K)
    sub_coeffs_err <- matrix(NA, n_samples, K)
    sub_fm0phi <- matrix(NA, n_samples, K)
    sub_fm1phi <- matrix(NA, n_samples, K)
    for (ir in seq_len(n_samples)) {          ## for each sampling
      for (m in seq_len(ncol(count_use))) {       ## for each cluster
        idx <- seq(1, nrow(count_use), n_samples) + ir - 1
        
        df_use <- data.frame(n1 = count_use[, m], total=rowSums(count_use))[idx,]
        df_use <- cbind(df_use, design_mat)
        df_tmp <- df_use[!is.na(design_mat[, k]), ]
        
        ## model fitting using betabin
        if (base_model=='NULL' | ncol(design_mat) == 1) {
          formula_fm0 <- as.formula('cbind(n1, total-n1) ~ 1')
          formula_fm1 <- as.formula(paste0('cbind(n1, total-n1)', '~ 1+', colnames(design_mat)[k], sep=''))
        } else if (base_model=='FULL') {
          fm0_right <- paste(colnames(design_mat)[-k], collapse = " + ")
          fm1_right <- paste(colnames(design_mat), collapse = " + ")
          formula_fm0 <- as.formula(paste0('cbind(n1, total-n1)', ' ~ 1 + ', fm0_right, sep=''))
          formula_fm1 <- as.formula(paste0('cbind(n1, total-n1)', ' ~ 1 + ', fm1_right, sep=''))
        }
        fm0 <- aod::betabin(formula_fm0, ~ 1, data = df_tmp, warnings = FALSE)
        fm1 <- aod::betabin(formula_fm1, ~ 1, data = df_tmp, warnings = FALSE)
        if (!is.null(fix_phi)){
          fm0 <- aod::betabin(formula_fm0, ~ 1, data = df_tmp, warnings = FALSE, fixpar = list(fm0@nbpar, fix_phi))
          fm1 <- aod::betabin(formula_fm1, ~ 1, data = df_tmp, warnings = FALSE, fixpar = list(fm1@nbpar, fix_phi))
        }
        
        ## ignore the fitting if the hessian matrix is singular
        if (length(fm1@varparam) < 4 || is.na(fm1@varparam[2, 2])) {next}
        
        
        sub_LR_val[ir, m] <- fm0@dev - fm1@dev
        sub_fm0phi[ir, m] <- tail(fm0@param, n=1)
        sub_fm1phi[ir, m] <- tail(fm1@param, n=1)
        parID <- grep(colnames(design_mat)[k], names(fm1@param))
        if(length(parID) > 1) stop("Please check the design matrix, make sure all factors are continous or categorical with only two levels.")
        sub_coeffs_val[ir, m] <- fm1@param[parID]
        if(is.null(fix_phi))
          sub_coeffs_err[ir, m] <-fm1@varparam[parID, parID]
      }
    }
    
    coeff_val_mean <- colMeans(sub_coeffs_val, na.rm = TRUE)
    if (is.null(fix_phi)){
      ## averaging the estimation to get the final result
      if (is.null(n_samples) || is.null(similarity_mat) || n_samples == 1){
        sub_coeff_err_pool <- colMeans(sub_coeffs_err**2, na.rm = TRUE)
      } else {
        sub_coeff_err_pool <- colMeans(sub_coeffs_err**2, na.rm = TRUE) +
          matrixStats::colSds(sub_coeffs_val) +
          matrixStats::colSds(sub_coeffs_val) / n_samples
      }
      # p values with Ward test: https://en.wikipedia.org/wiki/Wald_test
      pvals[,k] <- pnorm(-abs(coeff_val_mean) / sqrt(sub_coeff_err_pool))  * 2
      coeffs_err[, k] <- sqrt(sub_coeff_err_pool)
      fm0_phi[,k] <- colMeans(sub_fm0phi, na.rm = TRUE)
      fm1_phi[,k] <- colMeans(sub_fm1phi, na.rm = TRUE)
    }
    
    
    LR_median = robustbase::colMedians(sub_LR_val, na.rm = TRUE)
    LR_vals[, k] <- LR_median
    LRT_pvals[, k] <- pchisq(LR_median, df=1, lower.tail = FALSE, log.p = FALSE)
    coeffs[, k] <- coeff_val_mean
  }
  
  # Return list
  LRT_fdr[,] <- p.adjust(LRT_pvals, method = 'fdr')
  res <- list('ceoffs'=coeffs, 'coeffs_err'=coeffs_err, 'pvals' = pvals, 'LR_vals'=LR_vals, 'LRT_pvals'=LRT_pvals, 'fdr'=LRT_fdr, 'fm0_phi' = fm0_phi, 'fm1_phi' = fm1_phi)
  res
}

# function 15: get adjust Phi and average Phi
getPhi2 = function(count_mat, design_mat, similarity_mat=NULL, n_samples=NULL, pseudo_count=NULL,  base_model='NULL', fix_phi=NULL){
  res = eachPhi(count_mat, design_mat, similarity_mat, n_samples, pseudo_count,  base_model, fix_phi)
  avrgPhiV = res$fm0_phi %>% rowMeans()
  avrgPhi = mean(avrgPhiV)
  data = data.frame(phi = avrgPhiV, count = colSums(count_mat))
  fm = lm(phi~count, data)
  adjPhiV = (avrgPhiV + avrgPhi)/2
  return(list(adjPhiV = adjPhiV, avrgPhi = avrgPhi))
}

# function 16: new version of dcats that allow different input phi for different clusters
dcats_GLM_mltPhi <- function(count_mat, design_mat, similarity_mat=NULL, n_samples=50, pseudo_count=NULL,  base_model='NULL', phi=NULL) {
  if(!is.null(phi) & length(phi)!=ncol(count_mat)) stop("Please check the phi vector, make sure you provide over dispersion for all cell types.")
  
  # Output matrices
  coeffs     <- matrix(NA, ncol(count_mat), ncol(design_mat))
  coeffs_err <- matrix(NA, ncol(count_mat), ncol(design_mat))
  LR_vals    <- matrix(NA, ncol(count_mat), ncol(design_mat))
  LRT_pvals  <- matrix(NA, ncol(count_mat), ncol(design_mat))
  pvals  <- matrix(NA, ncol(count_mat), ncol(design_mat))
  LRT_fdr    <- matrix(NA, ncol(count_mat), ncol(design_mat))
  
  # Check colnames
  if (is.null(colnames(count_mat)))
    colnames(count_mat) <- paste0('cell_type_', seq(ncol(count_mat)))
  if (is.null(colnames(design_mat)))
    colnames(design_mat) <- paste0('factor_', seq(ncol(design_mat)))
  
  # Add rownames and colnames
  rownames(LR_vals) <- rownames(LRT_pvals) <- rownames(pvals) <- rownames(LRT_fdr) <-
    rownames(coeffs) <- rownames(coeffs_err) <- colnames(count_mat)
  colnames(LR_vals) <- colnames(LRT_pvals) <- colnames(pvals) <- colnames(LRT_fdr) <-
    colnames(coeffs) <- colnames(coeffs_err) <- colnames(design_mat)
  
  
  ## using estimated the latent cell counts
  count_latent = count_mat
  if(!is.null(similarity_mat)) {
    for (i in seq_len(nrow(count_mat))) {
      count_latent[i, ] <- sum(count_mat[i, ]) *
        multinom_EM(count_mat[i, ], similarity_mat, verbose = FALSE)$mu
    }
  }
  
  K <- ncol(count_mat) ## number of cell types
  if (is.null(similarity_mat)) {
    n_samples <- 1
  }
  
  if (!is.null(n_samples) && !is.null(similarity_mat)) {
    count_use <- matrix(0, nrow(count_mat) * n_samples, K)
    for (i in seq_len(nrow(count_mat))) {
      idx <- seq((i - 1) * n_samples + 1, i * n_samples)
      for (j in seq_len(K)) {
        count_use[idx, ] <- (
          count_use[idx, ] + t(rmultinom(n_samples, count_latent[i, j], similarity_mat[j, ])))
      }
    }
  } else{
    count_use <- count_latent
    n_samples <- 1
  }
  
  # adding pseudo counts
  if (is.null(pseudo_count)) {
    if (any(colMeans(count_mat) == 0)) {
      print(paste("Empty cell type exists in at least one conidtion;",
                  "adding replicate & condition specific pseudo count:"))
      count_use <- count_use + 1
    }
  } else {
    count_use = count_use + pseudo_count
  }
  
  count_use = round(count_use)
  
  #print(count_use)
  # Test each factor
  for (k in seq_len(ncol(design_mat))) {    ## for each factor
    sub_LR_val <- matrix(NA, n_samples, K)
    sub_coeffs_val <- matrix(NA, n_samples, K)
    sub_coeffs_err <- matrix(NA, n_samples, K)
    for (ir in seq_len(n_samples)) {          ## for each sampling
      for (m in seq_len(ncol(count_use))) {       ## for each cluster
        idx <- seq(1, nrow(count_use), n_samples) + ir - 1
        
        df_use <- data.frame(n1 = count_use[, m], total=rowSums(count_use))[idx,]
        df_use <- cbind(df_use, design_mat)
        df_tmp <- df_use[!is.na(design_mat[, k]), ]
        
        ## model fitting using betabin
        if (base_model=='NULL' | ncol(design_mat) == 1) {
          formula_fm0 <- as.formula('cbind(n1, total-n1) ~ 1')
          
          formula_fm1 <- as.formula(paste0('cbind(n1, total-n1)', '~ 1+', colnames(design_mat)[k], sep=''))
        } else if (base_model=='FULL') {
          fm0_right <- paste(colnames(design_mat)[-k], collapse = " + ")
          fm1_right <- paste(colnames(design_mat), collapse = " + ")
          formula_fm0 <- as.formula(paste0('cbind(n1, total-n1)', ' ~ 1 + ', fm0_right, sep=''))
          formula_fm1 <- as.formula(paste0('cbind(n1, total-n1)', ' ~ 1 + ', fm1_right, sep=''))
        }
        fm0 <- aod::betabin(formula_fm0, ~ 1, data = df_tmp, warnings = FALSE)
        fm1 <- aod::betabin(formula_fm1, ~ 1, data = df_tmp, warnings = FALSE)
        if (!is.null(phi)){
          fm0 <- aod::betabin(formula_fm0, ~ 1, data = df_tmp, warnings = FALSE, fixpar = list(fm0@nbpar, phi[m]))
          fm1 <- aod::betabin(formula_fm1, ~ 1, data = df_tmp, warnings = FALSE, fixpar = list(fm1@nbpar, phi[m]))
        }
        
        ## ignore the fitting if the hessian matrix is singular
        if (length(fm1@varparam) < 4 || is.na(fm1@varparam[2, 2])) {next}
        
        sub_LR_val[ir, m] <- fm0@dev - fm1@dev
        parID <- grep(colnames(design_mat)[k], names(fm1@param))
        if(length(parID) > 1) stop("Please check the design matrix, make sure all factors are continous or categorical with only two levels.")
        sub_coeffs_val[ir, m] <- fm1@param[parID]
        if(is.null(phi))
          sub_coeffs_err[ir, m] <-fm1@varparam[parID, parID]
      }
    }
    
    coeff_val_mean <- colMeans(sub_coeffs_val, na.rm = TRUE)
    if (is.null(phi)){
      ## averaging the estimation to get the final result
      if (is.null(n_samples) || is.null(similarity_mat) || n_samples == 1){
        sub_coeff_err_pool <- colMeans(sub_coeffs_err**2, na.rm = TRUE)
      } else {
        sub_coeff_err_pool <- colMeans(sub_coeffs_err**2, na.rm = TRUE) +
          matrixStats::colSds(sub_coeffs_val) +
          matrixStats::colSds(sub_coeffs_val) / n_samples
      }
      # p values with Ward test: https://en.wikipedia.org/wiki/Wald_test
      pvals[,k] <- pnorm(-abs(coeff_val_mean) / sqrt(sub_coeff_err_pool))  * 2
      coeffs_err[, k] <- sqrt(sub_coeff_err_pool)
    }
    
    
    LR_median = robustbase::colMedians(sub_LR_val, na.rm = TRUE)
    LR_vals[, k] <- LR_median
    LRT_pvals[, k] <- pchisq(LR_median, df=1, lower.tail = FALSE, log.p = FALSE)
    coeffs[, k] <- coeff_val_mean
  }
  
  # Return list
  LRT_fdr[,] <- p.adjust(LRT_pvals, method = 'fdr')
  res <- list('ceoffs'=coeffs, 'coeffs_err'=coeffs_err, 'pvals' = pvals, 'LR_vals'=LR_vals, 'LRT_pvals'=LRT_pvals, 'fdr'=LRT_fdr)
  res
}

# function 17: inverse of %in%
'%!in%' <- function(x,y) !('%in%'(x,y))

# function 18: simulator for milo(without removing batch effect)
simulator_noInt = function(totals1, totals2, probC1, probC2, setresolu, sim_mat){
  K = length(probC1)
  n_rep1 = length(totals1)
  n_rep2 = length(totals2)
  
  prop_cond1 = matrix(0, n_rep1, K)
  prop_cond2 = matrix(0, n_rep2, K)
  numb_cond1 <- matrix(0, n_rep1, K)
  numb_cond2 <- matrix(0, n_rep2, K)
  
  # cell selection
  slt_sim_matC1 = matrix(NA, ncol = 0, nrow = dim(sim_mat)[1])
  slt_origLabelsC1 = vector()
  slt_batchesC1 = character()
  for (rep in seq_len(n_rep1)) {
    prop_cond1[rep, ] <- MCMCpack::rdirichlet(1, probC1 * concentration)
    #prop_cond1[rep, ] = probC1
    numb_cond1[rep, ] <- rmultinom(1, totals1[rep], prop_cond1[rep, ])
    #numb_cond1[rep, ] = (totals1[rep] * prop_cond1[rep, ]) %>% ceiling()
    cell_slt = cell_slt_dup(numb_cond1[rep,], sim_mat, origLabels)
    slt_sim_matC1 = cbind(slt_sim_matC1, cell_slt$sub_sim_mat)
    slt_origLabelsC1 = c(slt_origLabelsC1, cell_slt$sub_origLabels)
    slt_batchesC1 = c(slt_batchesC1, rep(str_c("cond1s", as.character(rep)), length(cell_slt$sub_origLabels)))
  }
  
  slt_sim_matC2 = matrix(NA, ncol = 0, nrow = dim(sim_mat)[1])
  slt_origLabelsC2 = factor()
  slt_batchesC2 = character()
  for (rep in seq_len(n_rep2)) {
    prop_cond2[rep, ] <- MCMCpack::rdirichlet(1, probC2 * concentration)
    #prop_cond2[rep, ] = probC2
    numb_cond2[rep, ] <- rmultinom(1, totals2[rep], prop_cond2[rep, ])
    cell_slt = cell_slt_dup(numb_cond2[rep,], sim_mat, origLabels)
    slt_sim_matC2 = cbind(slt_sim_matC2, cell_slt$sub_sim_mat)
    slt_origLabelsC2 = c(slt_origLabelsC2, cell_slt$sub_origLabels)
    slt_batchesC2 = c(slt_batchesC2, rep(str_c("cond2s", as.character(rep)), length(cell_slt$sub_origLabels)))
  }
  
  print(numb_cond1)
  print(numb_cond2)
  
  ## seurat process
  ### pre-process
  seuratObj <- CreateSeuratObject(counts = cbind(slt_sim_matC1, slt_sim_matC2), project="Splatter")
  seuratObj <- AddMetaData(object = seuratObj, metadata = c(slt_batchesC1, slt_batchesC2), col.name = 'batch')
  seuratObj <- AddMetaData(object = seuratObj, metadata = c(rep("Cond1", length(slt_batchesC1)), rep("Cond2", length(slt_batchesC2))), col.name = 'condition')
  seuratObj <- NormalizeData(seuratObj, verbose = FALSE)
  #integratedSamples <- RunFastMNN(object.list = SplitObject(seuratObj, split.by = "batch"), verbose = FALSE)
  integratedSamples <- FindVariableFeatures(seuratObj, verbose = FALSE)
  integratedSamples <- ScaleData(integratedSamples, verbose = FALSE)
  integratedSamples <- RunPCA(integratedSamples, npcs = 30, verbose = FALSE)
  #integratedSamples <- RunUMAP(integratedSamples, reduction = "mnn", dims = 1:30, verbose = FALSE)
  #integratedSamples <- FindNeighbors(integratedSamples, reduction = "mnn", dims = 1:30, verbose = FALSE)
  integratedSamples <- RunUMAP(integratedSamples, dims = 1:30, verbose = FALSE)
  integratedSamples <- FindNeighbors(integratedSamples, dims = 1:30, verbose = FALSE)
  integratedSamples <- FindClusters(integratedSamples, resolution = setresolu, verbose = FALSE)
  
  ## decide resolution
  Kprep = integratedSamples@active.ident %>% as.factor() %>% summary() %>% length()
  str_c("Kprep: ", as.character(Kprep)) %>% print()
  str_c("setresolu: ", as.character(setresolu)) %>% print()
  while (Kprep != K & setresolu > 0.03) {
    if (Kprep > K){
      setresolu = setresolu - 0.03
      integratedSamples <- FindClusters(integratedSamples, resolution = setresolu, verbose = FALSE)
      Kprep = integratedSamples@active.ident %>% as.factor() %>% summary() %>% length()
      str_c("Kprep: ", as.character(Kprep)) %>% print()
      str_c("setresolu: ", as.character(setresolu)) %>% print()
    } else {
      setresolu = setresolu + 0.01
      integratedSamples <- FindClusters(integratedSamples, resolution = setresolu, verbose = FALSE)
      Kprep = integratedSamples@active.ident %>% as.factor() %>% summary() %>% length()
      str_c("Kprep: ", as.character(Kprep)) %>% print()
      str_c("setresolu: ", as.character(setresolu)) %>% print()
    }
  }
  ### change labels to A, B, C
  integratedSamples@active.ident = integratedSamples@active.ident %>% 
    plyr::mapvalues(from = c(0:15), to = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P"))
  
  # get count data
  dfRes = data.frame(clusterRes = integratedSamples@active.ident, batch = integratedSamples$batch, condition = integratedSamples$condition) %>% 
    tibble::rownames_to_column("cellID")
  if(Kprep == K){
    Res = list(integratedSamples = integratedSamples, dfRes = dfRes, trueLabels = c(slt_origLabelsC1, slt_origLabelsC2), numb_cond1 = numb_cond1, numb_cond2 = numb_cond2, prop_cond1 = prop_cond1, prop_cond2 = prop_cond2)
    return(Res)
  } else {
    return(NA)
  }
  
}
## function to get milo results
getMilo = function(integratedSamples){
  # create object
  sceObj <- as.SingleCellExperiment(integratedSamples)
  milo <- Milo(sceObj)
  milo <- buildGraph(milo, k = 10, d = 30)
  milo <- makeNhoods(milo, prop = 0.1, k = 10, d=30, refined = TRUE)
  # get milo result
  milo <- countCells(milo, meta.data = data.frame(colData(milo)), sample="batch")
  if ("gender" %in% colnames(colData(milo))) {
    designDF <- data.frame(colData(milo))[,c("batch", "condition", "gender", "age")]
    designDF <- distinct(designDF)
    milo <- calcNhoodDistance(milo, d=30)
    rownames(designDF) <- designDF$batch
    da_results <- testNhoods(milo, design = ~ age + gender + condition, design.df = designDF)
  } else {
    designDF <- data.frame(colData(milo))[,c("batch", "condition")]
    designDF <- distinct(designDF)
    milo <- calcNhoodDistance(milo, d=30)
    rownames(designDF) <- designDF$batch
    da_results <- testNhoods(milo, design = ~ condition, design.df = designDF)
  }
  da_results <- annotateNhoods(milo, da_results, coldata_col = "ident")
  counts <- nhoodCounts(milo)
  # get nhoods
  anno_vec <- colData(milo) %>% rownames()
  if (!is.factor(anno_vec)) {
    anno_vec <- factor(anno_vec, levels=unique(anno_vec))
  }
  
  ## Count occurrence of labels in each nhood
  n.levels <- length(levels(anno_vec))
  nhood_counts <- vapply(seq_len(ncol(nhoods(milo))), FUN=function(n)
    table(anno_vec[which(nhoods(milo)[,n]==1)]),FUN.VALUE=numeric(n.levels))
  colnames(nhood_counts) = str_c("hoods", 1:ncol(nhood_counts))
  select_nhood = nhood_counts[rowSums(nhood_counts) > 0,]
  nhoodsID = sapply(seq_len(nrow(select_nhood)), FUN = function(n)
    sample(colnames(select_nhood)[select_nhood[n,] != 0], size = 1))
  names(nhoodsID) = rownames(select_nhood)
  #nhoodsInfo = data.frame(cellID = rownames(select_nhood), nhoodsID = hoodsID)
  return(list(da_results = da_results, counts = counts, nhoodsID = nhoodsID))
}

## function to select cells for simulation with multi-covariates
simulator_propM = function(totalV, prop_mat, setresolu, sim_mat, design_mat){
  K = ncol(prop_mat)
  prop_sim = matrix(NA, nrow = nrow(prop_mat), ncol = ncol(prop_mat))
  numb_sim = matrix(NA, nrow = nrow(prop_mat), ncol = ncol(prop_mat))
  
  # cell selection
  slt_sim_mat = matrix(NA, ncol = 0, nrow = dim(sim_mat)[1])
  slt_origLabels = vector()
  slt_batches = character()
  slt_conditions = character()
  slt_genders = character()
  slt_ages = vector()
  for (idx in 1:nrow(prop_mat)) {
    prop_sim[idx, ] <- MCMCpack::rdirichlet(1, prop_mat[idx,] * concentration)
    numb_sim[idx, ] <- rmultinom(1, totalV[idx], prop_sim[idx, ])
    cell_slt = cell_slt_dup(numb_sim[idx, ], sim_mat, origLabels)
    slt_sim_mat = cbind(slt_sim_mat, cell_slt$sub_sim_mat)
    slt_origLabels = c(slt_origLabels, cell_slt$sub_origLabels)
    slt_batches = c(slt_batches, rep(str_c("rep", as.character(idx)), length(cell_slt$sub_origLabels)))
    slt_conditions = c(slt_conditions, rep(design_mat[idx,1], length(cell_slt$sub_origLabels)))
    slt_genders = c(slt_genders, rep(design_mat[idx,2], length(cell_slt$sub_origLabels)))
    slt_ages = c(slt_ages, rep(design_mat[idx,3], length(cell_slt$sub_origLabels)))
  }
  
  print(numb_sim)
  
  ## seurat process
  ### pre-process
  seuratObj <- CreateSeuratObject(counts = slt_sim_mat, project="Splatter")
  seuratObj <- AddMetaData(object = seuratObj, metadata = slt_batches, col.name = 'batch')
  seuratObj <- AddMetaData(object = seuratObj, metadata = slt_conditions, col.name = 'condition')
  seuratObj <- AddMetaData(object = seuratObj, metadata = slt_genders, col.name = 'gender')
  seuratObj <- AddMetaData(object = seuratObj, metadata = slt_ages, col.name = 'age')
  seuratObj <- NormalizeData(seuratObj, verbose = FALSE)
  integratedSamples <- FindVariableFeatures(seuratObj, verbose = FALSE)
  integratedSamples <- ScaleData(integratedSamples, verbose = FALSE)
  integratedSamples <- RunPCA(integratedSamples, npcs = 30, verbose = FALSE, features = VariableFeatures(object = integratedSamples))
  integratedSamples <- RunUMAP(integratedSamples, dims = 1:30, verbose = FALSE)
  integratedSamples <- FindNeighbors(integratedSamples, dims = 1:30, verbose = FALSE)
  integratedSamples <- FindClusters(integratedSamples, resolution = setresolu, verbose = FALSE)
  
  ## decide resolution
  Kprep = integratedSamples@active.ident %>% as.factor() %>% summary() %>% length()
  str_c("Kprep: ", as.character(Kprep)) %>% print()
  str_c("setresolu: ", as.character(setresolu)) %>% print()
  while (Kprep != K & setresolu > 0.03) {
    if (Kprep > K){
      setresolu = setresolu - 0.03
      integratedSamples <- FindClusters(integratedSamples, resolution = setresolu, verbose = FALSE)
      Kprep = integratedSamples@active.ident %>% as.factor() %>% summary() %>% length()
      str_c("Kprep: ", as.character(Kprep)) %>% print()
      str_c("setresolu: ", as.character(setresolu)) %>% print()
    } else {
      setresolu = setresolu + 0.01
      integratedSamples <- FindClusters(integratedSamples, resolution = setresolu, verbose = FALSE)
      Kprep = integratedSamples@active.ident %>% as.factor() %>% summary() %>% length()
      str_c("Kprep: ", as.character(Kprep)) %>% print()
      str_c("setresolu: ", as.character(setresolu)) %>% print()
    }
  }
  ### change labels to A, B, C
  integratedSamples@active.ident = integratedSamples@active.ident %>% 
    plyr::mapvalues(from = c(0:15), to = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P"))
  
  # get count data
  dfRes = data.frame(clusterRes = integratedSamples@active.ident, batch = integratedSamples$batch, condition = integratedSamples$condition) %>% 
    tibble::rownames_to_column("cellID")
  if(Kprep == K){
    Res = list(integratedSamples = integratedSamples, dfRes = dfRes, trueLabels = slt_origLabels, numb_sim = numb_sim, prop_sim = prop_sim)
    return(Res)
  } else {
    return(NA)
  }
}