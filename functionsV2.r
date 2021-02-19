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
  cell_numO = sort(cell_num, decreasing = TRUE)
  cellname = as.factor(origLabels) %>% 
    fct_infreq() %>% 
    levels()
  K = length(cell_numO)
  index = numeric()
  for (j in 1:K){
    group_idx = c(1:length(origLabels))[origLabels == cellname[j]] %>% 
      sample(cell_numO[j])
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
    plyr::mapvalues(from = c(0:11), to = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L"))
  
  return(integratedSamples)
}

# function 6: fastMNN
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
    plyr::mapvalues(from = c(0:11), to = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L"))
  
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
  for (i in seq_len(n_rep1)) {
    prop_cond1[i, ] <- MCMCpack::rdirichlet(1, probC1 * concentration)
    numb_cond1[i, ] <- rmultinom(1, totals1[i], prop_cond1[i, ])
    #numb_cond1[i, ] = (totals1[i] * prop_cond1[i, ]) %>% ceiling()
    cell_slt = cell_slt_dup(numb_cond1[i,], sim_mat, origLabels)
    slt_sim_matC1 = cbind(slt_sim_matC1, cell_slt$sub_sim_mat)
    slt_origLabelsC1 = c(slt_origLabelsC1, cell_slt$sub_origLabels)
    slt_batchesC1 = c(slt_batchesC1, rep(str_c("cond1s", as.character(i)), length(cell_slt$sub_origLabels)))
  }
  
  slt_sim_matC2 = matrix(NA, ncol = 0, nrow = dim(sim_mat)[1])
  slt_origLabelsC2 = factor()
  slt_batchesC2 = character()
  for (i in seq_len(n_rep2)) {
    prop_cond2[i, ] <- MCMCpack::rdirichlet(1, probC2 * concentration)
    numb_cond2[i, ] <- rmultinom(1, totals2[i], prop_cond2[i, ])
    cell_slt = cell_slt_dup(numb_cond2[i,], sim_mat, origLabels)
    slt_sim_matC2 = cbind(slt_sim_matC2, cell_slt$sub_sim_mat)
    slt_origLabelsC2 = c(slt_origLabelsC2, cell_slt$sub_origLabels)
    slt_batchesC2 = c(slt_batchesC2, rep(str_c("cond2s", as.character(i)), length(cell_slt$sub_origLabels)))
  }
  
  print(numb_cond1)
  print(numb_cond2)
  integratedSamples = fastMNN(slt_sim_matC1, slt_sim_matC2, slt_batchesC1, slt_batchesC2, setresolu)
  Kprep = integratedSamples@active.ident %>% as.factor() %>% summary() %>% length()
  str_c("Kprep: ", as.character(Kprep)) %>% print()
  str_c("setresolu: ", as.character(setresolu)) %>% print()
  while (Kprep != K & setresolu > 0.03) {
  #while (Kprep != K) {
    if (Kprep > K){
      setresolu = setresolu - 0.03
      #if (setresolu <= 0) {setresolu = setresolu + 0.02}
      integratedSamples = fastMNN(slt_sim_matC1, slt_sim_matC2, slt_batchesC1, slt_batchesC2, setresolu)
      Kprep = integratedSamples@active.ident %>% as.factor() %>% summary() %>% length()
      str_c("Kprep: ", as.character(Kprep)) %>% print()
      str_c("setresolu: ", as.character(setresolu)) %>% print()
    } else {
      setresolu = setresolu + 0.01
      integratedSamples = fastMNN(slt_sim_matC1, slt_sim_matC2, slt_batchesC1, slt_batchesC2, setresolu)
      Kprep = integratedSamples@active.ident %>% as.factor() %>% summary() %>% length()
      str_c("Kprep: ", as.character(Kprep)) %>% print()
      str_c("setresolu: ", as.character(setresolu)) %>% print()
    }
  }
  
  # get count data
  dfRes = data.frame(clusterRes = integratedSamples@active.ident, batch = integratedSamples$batch, condition = integratedSamples$condition) %>% 
    tibble::rownames_to_column("cellID")
  if(Kprep == K){
    Res = list(integratedSamples = integratedSamples, dfRes = dfRes, trueLabels = c(slt_origLabelsC1, slt_origLabelsC2), numb_cond1 = numb_cond1, numb_cond2 = numb_cond2)
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
  diffcytP = topTable(res_DA, format_vals = TRUE) %>% 
    as.data.frame() %>% 
    rename(cluster = cluster_id) %>% 
    rename(diffcyt_pvals = p_adj) %>% 
    dplyr::select(cluster, diffcyt_pvals)
  
  return(diffcytP)
}

# function 10: for sending email
sendEmail = function(subject = "The job is done"){
  library(emayili)
  email <- envelope() %>%
    from("dummy_acc29@outlook.com") %>%
    to("xl2836@outlook.com") %>%
    subject(subject) %>%
    text("Yeah!!!!!")
  smtp <- server(host = "smtp.gmail.com", port = 465, username = "4sendEmail29@gmail.com", password = "dummy_acc29")
  smtp(email, verbose = FALSE)
}

# function 11: rewrite dcats_betabin to find out matrix causes the @fm1 error
dcats_betabinRW <- function(counts1, counts2, similarity_mat=NULL, n_samples=50,
                          pseudo_count=NULL, binom_only=TRUE) {
  ## Check counts1 and counts2 shape
  if (length(counts1) == 1 || is.null(dim(counts1)) ||
      length(dim(counts1)) < 2) {
    counts1 = matrix(counts1, nrow=1)
  }
  if (length(counts2) == 1 || is.null(dim(counts2)) ||
      length(dim(counts2)) < 2) {
    counts2 = matrix(counts2, nrow=1)
  }
  
  prop1 <- counts1 / rowSums(counts1)
  prop2 <- counts2 / rowSums(counts2)
  
  ## using estimated the latent cell counts
  counts1_latent = counts1
  counts2_latent = counts2
  for (i in seq_len(nrow(counts1))) {
    counts1_latent[i, ] <- sum(counts1[i, ]) *
      multinom_EM(counts1[i, ], similarity_mat, verbose = FALSE)$mu
  }
  for (i in seq_len(nrow(counts2))) {
    counts2_latent[i, ] <- sum(counts2[i, ]) *
      multinom_EM(counts2[i, ], similarity_mat, verbose = FALSE)$mu
  }
  
  ## number of cell types
  K <- ncol(counts1)
  if (is.null(similarity_mat)) {
    n_samples <- 1
  }
  
  ##------------------------ for debugging only ----------------
  #set.seed(123)
  ##------------------------------------------------------------
  
  if (!is.null(n_samples) && !is.null(similarity_mat)) {
    counts1_use <- matrix(0, nrow(counts1) * n_samples, K)
    counts2_use <- matrix(0, nrow(counts2) * n_samples, K)
    for (i in seq_len(nrow(counts1))) {
      idx <- seq((i - 1) * n_samples + 1, i * n_samples)
      for (j in seq_len(K)) {
        counts1_use[idx, ] <- (
          counts1_use[idx, ] + t(rmultinom(n_samples,
                                           counts1_latent[i, j],
                                           similarity_mat[j, ])))
      }
    }
    for (i in seq_len(nrow(counts2))) {
      idx <- seq((i - 1) * n_samples + 1, i * n_samples)
      for (j in seq_len(K)) {
        counts2_use[idx, ] <- (
          counts2_use[idx, ] + t(rmultinom(n_samples,
                                           counts2_latent[i, j],
                                           similarity_mat[j, ])))
      }
    }
  } else{
    counts1_use <- counts1
    counts2_use <- counts2
  }
  
  # adding pseudo counts
  if (is.null(pseudo_count)) {
    if (any(colMeans(counts1) == 0) || any(colMeans(counts2) == 0) ) {
      print(paste("Empty cell type exists in at least one conidtion;",
                  "adding replicate & condition specific pseudo count:"))
      counts1_use <- counts1_use + 1
      counts2_use <- counts2_use + 1
    }
  } else {
    counts2_use = counts2_use + pseudo_count
    counts2_use = counts2_use + pseudo_count
  }
  
  ## Binomial regression for each sampling
  LR_val <- matrix(NA, n_samples, K)
  coeffs_val <- matrix(NA, n_samples, K)
  coeffs_err <- matrix(NA, n_samples, K)
  intercept_val <- matrix(NA, n_samples, K)
  intercept_err <- matrix(NA, n_samples, K)
  total_all <- c(rowSums(counts1_use), rowSums(counts2_use))
  label_all <- c(rep(1, nrow(counts1_use)), rep(0, nrow(counts2_use)))
  for (ir in seq_len(n_samples)) {
    idx <- seq(1, length(total_all), n_samples) + ir - 1
    for (i in seq_len(K)) {
      n1 <- c(counts1_use[, i], counts2_use[, i])[idx]
      df <- data.frame(n1 = n1, n2 = total_all[idx] - n1,
                       label = label_all[idx])
      if (binom_only) {
        model1 <- glm(cbind(n1, n2) ~ label + 1,
                      family = binomial(), data = df)
        coeffs_val[ir, i] <- summary(model1)$coefficients[2, 1]
        coeffs_err[ir, i] <- summary(model1)$coefficients[2, 2]
        intercept_val[ir, i] <- summary(model1)$coefficients[1, 1]
        intercept_err[ir, i] <- summary(model1)$coefficients[1, 2]
        
        model0 <- glm(cbind(n1, n2) ~ 1,
                      family = binomial(), data = df)
        LR_val[ir, i] <- model0$deviance - model1$deviance
      } else{
        
        fm1 <- aod::betabin(cbind(n1, n2) ~ label, ~ 1, data = df)
        if (any(is.na(fm1@varparam))) {
          print(fm1)
          print(df)
        } else {
          coeffs_err[ir, i] <-fm1@varparam[2, 2]
        }
        
        coeffs_val[ir, i] <- fm1@param[2]
        intercept_val[ir, i] <- fm1@param[1] # summary(fm1)@Coef[1, 1]
        intercept_err[ir, i] <- fm1@varparam[1, 1]
        
        fm0 <- aod::betabin(cbind(n1, n2) ~ 1, ~ 1, data = df)
        LR_val[ir, i] <- fm0@dev - fm1@dev
      }
      
    }
  }
  
  ## Averaging the coeffcients errors
  if (is.null(n_samples) || is.null(similarity_mat) || n_samples == 1) {
    coeff_val_mean <- colMeans(coeffs_val)
    coeff_err_pool <- colMeans(coeffs_err**2)
    intercept_val_mean <- colMeans(intercept_val)
    intercept_err_pool <- colMeans(intercept_err**2)
  } else{
    coeff_val_mean <- colMeans(coeffs_val)
    coeff_err_pool <- colMeans(coeffs_err**2) +
      matrixStats::colSds(coeffs_val) +
      matrixStats::colSds(coeffs_val) / n_samples
    
    intercept_val_mean <- colMeans(intercept_val)
    intercept_err_pool <- colMeans(intercept_err**2) +
      matrixStats::colSds(intercept_val) +
      matrixStats::colSds(intercept_val) / n_samples
  }
  
  
  # p values with Ward test: https://en.wikipedia.org/wiki/Wald_test
  pvals <- pnorm(-abs(coeff_val_mean) / sqrt(coeff_err_pool))  * 2
  
  LR_median = robustbase::colMedians(LR_val)
  LRT_pvals <-  pchisq(LR_median, df=1, lower.tail = FALSE, log.p = FALSE)
  
  data.frame(
    "prop1_mean" = colMeans(prop1),
    "prop1_std"  = matrixStats::colSds(prop1),
    "prop2_mean" = colMeans(prop2),
    "prop2_std"  = matrixStats::colSds(prop2),
    "coeff_mean" = coeff_val_mean,
    "coeff_std"  = sqrt(coeff_err_pool),
    "intecept_mean" = intercept_val_mean,
    "intecept_std"  = sqrt(intercept_err_pool),
    "pvals" = pvals,
    "LRT_pvals" = LRT_pvals,
    row.names = colnames(counts1)
  )
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
