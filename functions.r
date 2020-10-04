if(FALSE){
  library(splatter)
  library(Seurat)
  library(speckle)
  library(DCATS)
  library(ggplot2)
  library(tidyverse)
  library(diffcyt)
}

# function 1: simulation
simualtion = function(probNor, probMut, de_prob, batch_size){
  set.seed(12345)
  
  # simulate normal
  param.groups <- newSplatParams(batchCells = c(batch_size, batch_size, batch_size), nGenes = 100)
  simNor <- splatSimulateGroups(param.groups, group.prob = probNor, de.prob = de_prob, verbose = FALSE)
  simNor@colData@rownames = str_replace(simNor@colData@rownames, "Cell", "NorCell")
  simNor_mat <- counts(simNor)
  
  # simulate mutate
  simMut <- splatSimulateGroups(param.groups, group.prob = probMut, de.prob = de_prob, verbose = FALSE)
  simMut@colData@rownames = str_replace(simMut@colData@rownames, "Cell", "MutCell")
  simMut_mat <- counts(simMut)
  
  # batch information
  batchNor = simNor@colData@listData$Batch %>% 
    str_replace("Batch", "Nor")
  batchMut = simMut@colData@listData$Batch %>% 
    str_replace("Batch", "Mut")
  
  origLabels = c(simNor@colData@listData$Group, simMut@colData@listData$Group)
  
  return(list(simNor_mat = simNor_mat, simMut_mat = simMut_mat, batchNor = batchNor, batchMut = batchMut, origLabels = origLabels))
}

# function 2: Seurat process(might need few minutes)
runSeurat = function(sim_list, batch_size, setresolu){
  # set input
  simNor_mat = sim_list$simNor_mat
  simMut_mat = sim_list$simMut_mat
  batchNor = sim_list$batchNor
  batchMut = sim_list$batchMut
  
  # pre-process
  seuratNor <- CreateSeuratObject(counts = simNor_mat, project="Splatter")
  seuratNor <- AddMetaData(object = seuratNor, metadata = batchNor, col.name = 'batch')
  seuratNor <- AddMetaData(object = seuratNor, metadata = rep("Normal", batch_size*3), col.name = 'condition')
  seuratMut <- CreateSeuratObject(counts = simMut_mat, project="Splatter")
  seuratMut <- AddMetaData(object = seuratMut, metadata = batchMut, col.name = 'batch')
  seuratMut <- AddMetaData(object = seuratMut, metadata = rep("Mutate", batch_size*3), col.name = 'condition')
  
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
  integratedSamples <- FindNeighbors(integratedSamples, dims = 1:5, verbose = FALSE)
  integratedSamples <- FindClusters(integratedSamples, resolution = setresolu, algorithm=2, verbose = FALSE)
  integratedSamples <- RunUMAP(integratedSamples, reduction = "pca", dims = 1:30, verbose = FALSE)
  
  # change labels to A, B, C
  integratedSamples@active.ident = integratedSamples@active.ident %>% 
    plyr::mapvalues(from = c("0", "1", "2", "3", "4", "5"), to = c("A", "B", "C", "D", "E", "F"))
  
  return(integratedSamples)
}

# function 3: function to get p-values and time(Fisher, speckle, dcats)
getPandTimeFSD = function(integratedSamples, sim_list){
  time = rep(NA,3)
  dfRes = data.frame(clusterRes = integratedSamples@active.ident, batch = integratedSamples$batch, condition = integratedSamples$condition) %>% 
    tibble::rownames_to_column("cellID")
  
  ## Fisher's exact test
  dfCount = dfRes %>% 
    group_by(condition, clusterRes) %>% 
    summarise(n = n()) %>% 
    pivot_wider(names_from = "clusterRes", values_from = "n") %>% 
    mutate(nonA = batch_size*6 - A,
           nonB = batch_size*6 - B,
           nonC = batch_size*6 - C,
           nonD = batch_size*6 - D,
           nonE = batch_size*6 - E,
           nonF = batch_size*6 - F)
  t1start = Sys.time()
  fisher_pvals = rep(NA,3)
  for (i in 1:6){
    fisher_pvals[i] = fisher.test(dfCount[,c(i+1,i+7)])$p.value
  }
  time[1] = Sys.time() - t1start
  fisherP = data.frame(cluster = c("A", "B", "C", "D", "E", "F"), fisher_pvals = fisher_pvals)
  
  ## speckle
  t2start = Sys.time()
  speckleRes = propeller(clusters = dfRes$clusterRes, sample = dfRes$batch, group = dfRes$condition)
  time[2] = Sys.time() - t2start
  #print(speckleRes)
  speckleP = data.frame(cluster = speckleRes$BaselineProp.clusters, speckle_pvals = speckleRes$FDR)
  
  ## DCATS
  celllabels_orig = sim_list$origLabels
  conf.mat<-table(Idents(integratedSamples), celllabels_orig)
  true.conf<-t(t(conf.mat)/apply(conf.mat,2,sum))
  library(clue)
  order = as.vector(solve_LSAP(true.conf, maximum = TRUE))
  #print(true.conf)
  condition = integratedSamples@meta.data$condition
  condNor<-Idents(integratedSamples)[condition == "Normal"]
  condMut<-Idents(integratedSamples)[condition == "Mutate"]
  countNor = table(sim_list$batchNor, condNor)
  countMut = table(sim_list$batchMut, condMut)
  t3start = Sys.time()
  dcatsResT = dcats_fit(countNor, countMut, true.conf[,order])  # true similarity matrix
  time[3] = Sys.time() - t3start
  dcatsResI = dcats_fit(countNor, countMut, diag(6))  # identity matrix
  dcatsResU = dcats_fit(countNor, countMut, get_similarity_mat(6, 0.05)) # uniform matrix
  graphs = integratedSamples@graphs$integrated_snn
  labels = Idents(integratedSamples)
  knn_mat = KNN_transition(graphs, labels)
  order = colnames(countNor)
  dcatsResK = dcats_fit(countNor, countMut, knn_mat[order, order])
  #print(dcatsRes)
  dcatsPT = data.frame(cluster = rownames(dcatsResT), dcats_pvalsT = dcatsResT$pvals)
  dcatsPI = data.frame(cluster = rownames(dcatsResI), dcats_pvalsI = dcatsResI$pvals)
  dcatsPU = data.frame(cluster = rownames(dcatsResU), dcats_pvalsU = dcatsResU$pvals)
  dcatsPK = data.frame(cluster = rownames(dcatsResK), dcats_pvalsK = dcatsResK$pvals)
  
  ## results
  Res_df = merge(dcatsPT, speckleP, by = "cluster") %>% 
    merge(dcatsPI, by = "cluster") %>% 
    merge(dcatsPU, by = "cluster") %>% 
    merge(dcatsPK, by = "cluster") %>%
    merge(fisherP, by = "cluster")
  time_df = data.frame(methods = c("fisher", "sepckle", "dcats"), time = time)
  return(list(Res_df = Res_df, time_df = time_df))
}

# function 4: function to get p-values and time(Fisher, speckle, dcats, diffcyt)
getPandTimeFSDD = function(integratedSamples, sim_list, batch_size = 1000){
  time = rep(NA,4)
  dfRes = data.frame(clusterRes = integratedSamples@active.ident, batch = integratedSamples$batch, condition = integratedSamples$condition) %>% 
    tibble::rownames_to_column("cellID")
  
  ## Fisher's exact test
  dfCount = dfRes %>% 
    group_by(condition, clusterRes) %>% 
    summarise(n = n()) %>% 
    pivot_wider(names_from = "clusterRes", values_from = "n") %>% 
    mutate(nonA = batch_size*6 - A,
           nonB = batch_size*6 - B,
           nonC = batch_size*6 - C,
           nonD = batch_size*6 - D,
           nonE = batch_size*6 - E,
           nonF = batch_size*6 - F)
  t1start = Sys.time()
  fisher_pvals = rep(NA,3)
  for (i in 1:6){
    fisher_pvals[i] = fisher.test(dfCount[,c(i+1,i+7)])$p.value
  }
  time[1] = Sys.time() - t1start
  fisherP = data.frame(cluster = c("A", "B", "C", "D", "E", "F"), fisher_pvals = fisher_pvals)
  
  ## speckle
  t2start = Sys.time()
  speckleRes = propeller(clusters = dfRes$clusterRes, sample = dfRes$batch, group = dfRes$condition)
  time[2] = Sys.time() - t2start
  #print(speckleRes)
  speckleP = data.frame(cluster = speckleRes$BaselineProp.clusters, speckle_pvals = speckleRes$FDR)
  
  ## DCATS
  celllabels_orig = sim_list$origLabels
  conf.mat<-table(Idents(integratedSamples), celllabels_orig)
  true.conf<-t(t(conf.mat)/apply(conf.mat,2,sum))
  library(clue)
  order = as.vector(solve_LSAP(true.conf, maximum = TRUE))
  #print(true.conf)
  condition = integratedSamples@meta.data$condition
  condNor<-Idents(integratedSamples)[condition == "Normal"]
  condMut<-Idents(integratedSamples)[condition == "Mutate"]
  countNor = table(sim_list$batchNor, condNor)
  countMut = table(sim_list$batchMut, condMut)
  t3start = Sys.time()
  dcatsResT = dcats_fit(countNor, countMut, true.conf[,order])  # true similarity matrix
  time[3] = Sys.time() - t3start
  dcatsResI = dcats_fit(countNor, countMut, diag(6))  # identity matrix
  dcatsResU = dcats_fit(countNor, countMut, get_similarity_mat(6, 0.05)) # uniform matrix
  graphs = integratedSamples@graphs$integrated_snn
  labels = Idents(integratedSamples)
  knn_mat = KNN_transition(graphs, labels)
  order = colnames(countNor)
  dcatsResK = dcats_fit(countNor, countMut, knn_mat[order, order])
  #print(dcatsRes)
  dcatsPT = data.frame(cluster = rownames(dcatsResT), dcats_pvalsT = dcatsResT$pvals)
  dcatsPI = data.frame(cluster = rownames(dcatsResI), dcats_pvalsI = dcatsResI$pvals)
  dcatsPU = data.frame(cluster = rownames(dcatsResU), dcats_pvalsU = dcatsResU$pvals)
  dcatsPK = data.frame(cluster = rownames(dcatsResK), dcats_pvalsK = dcatsResK$pvals)
  
  ## diffcyt
  # Function to create random data (one sample)
  d_random <- function(n = 20*batch_size, mean = 0, sd = 1, ncol = 20, cofactor = 5) {
    d <- sinh(matrix(rnorm(n, mean, sd), ncol = ncol)) * cofactor
    colnames(d) <- paste0("marker", sprintf("%02d", 1:ncol))
    d
  }
  # Create random data (without differential signal)
  set.seed(123)
  d_input <- list(
    Nor1 = d_random(), 
    Nor2 = d_random(), 
    Nor3 = d_random(), 
    Mut1 = d_random(),
    Mut2 = d_random(),
    Mut3 = d_random()
  )
  experiment_info <- data.frame(
    sample_id = factor(c("Nor1", "Nor2", "Nor3", "Mut1", "Mut2", "Mut3")), 
    group_id = factor(c(rep("Normal", 3), rep("Mutate", 3))), 
    stringsAsFactors = FALSE
  )
  marker_info <- data.frame(
    channel_name = paste0("channel", sprintf("%03d", 1:20)), 
    marker_name = paste0("marker", sprintf("%02d", 1:20)), 
    marker_class = factor(c(rep("type", 10), rep("state", 10)), 
                          levels = c("type", "state", "none")), 
    stringsAsFactors = FALSE
  )
  # Prepare data
  d_se <- prepareData(d_input, experiment_info, marker_info)
  # Transform data
  d_se <- transformData(d_se)
  # Generate clusters
  d_se <- generateClusters(d_se)
  # Replace with our simulated data info
  d_se@elementMetadata$sample_id = integratedSamples$batch
  d_se@elementMetadata$group_id = integratedSamples$condition
  d_se@elementMetadata$cluster_id = integratedSamples@active.ident
  # Calculate counts
  d_counts <- calcCounts(d_se)
  # Create design matrix
  design <- createDesignMatrix(experiment_info, cols_design = "group_id")
  # Create contrast matrix
  contrast <- createContrast(c(0, 1))
  # Test for differential abundance (DA) of clusters
  t4start = Sys.time()
  res_DA <- testDA_edgeR(d_counts, design, contrast)
  time[4] = Sys.time() - t4start
  diffcytP = topTable(res_DA, format_vals = TRUE) %>% 
    as.data.frame() %>% 
    rownames_to_column("cluster") %>% 
    rename(diffcyt_pvals = p_adj) %>% 
    select(cluster, diffcyt_pvals)
  
  ## results
  Res_df = merge(dcatsPT, speckleP, by = "cluster") %>% 
    merge(dcatsPI, by = "cluster") %>% 
    merge(dcatsPU, by = "cluster") %>% 
    merge(dcatsPK, by = "cluster") %>%
    merge(fisherP, by = "cluster") %>% 
    merge(diffcytP, by = "cluster")
  time_df = data.frame(methods = c("fisher", "sepckle", "dcats", "diffcyt"), time = time)
  return(list(Res_df = Res_df, time_df = time_df))
}

# function 5: function to get p-values and time(Fisher, speckle, dcats, diffcyt) for three clusters
getPandTimeFSDD3 = function(integratedSamples, sim_list, batch_size = 1000){
  time = rep(NA,4)
  dfRes = data.frame(clusterRes = integratedSamples@active.ident, batch = integratedSamples$batch, condition = integratedSamples$condition) %>% 
    tibble::rownames_to_column("cellID")
  
  ## Fisher's exact test
  dfCount = dfRes %>% 
    group_by(condition, clusterRes) %>% 
    summarise(n = n()) %>% 
    pivot_wider(names_from = "clusterRes", values_from = "n") %>% 
    mutate(nonA = batch_size*6 - A,
           nonB = batch_size*6 - B,
           nonC = batch_size*6 - C)
  t1start = Sys.time()
  fisher_pvals = rep(NA,3)
  for (i in 1:3){
    fisher_pvals[i] = fisher.test(dfCount[,c(i+1,i+4)])$p.value
  }
  time[1] = Sys.time() - t1start
  fisherP = data.frame(cluster = c("A", "B", "C"), fisher_pvals = fisher_pvals)
  
  ## speckle
  t2start = Sys.time()
  speckleRes = propeller(clusters = dfRes$clusterRes, sample = dfRes$batch, group = dfRes$condition)
  time[2] = Sys.time() - t2start
  #print(speckleRes)
  speckleP = data.frame(cluster = speckleRes$BaselineProp.clusters, speckle_pvals = speckleRes$FDR)
  
  ## DCATS
  celllabels_orig = sim_list$origLabels
  conf.mat<-table(Idents(integratedSamples), celllabels_orig)
  true.conf<-t(t(conf.mat)/apply(conf.mat,2,sum))
  library(clue)
  order = as.vector(solve_LSAP(true.conf, maximum = TRUE))
  #print(true.conf)
  condition = integratedSamples@meta.data$condition
  condNor<-Idents(integratedSamples)[condition == "Normal"]
  condMut<-Idents(integratedSamples)[condition == "Mutate"]
  countNor = table(sim_list$batchNor, condNor)
  countMut = table(sim_list$batchMut, condMut)
  t3start = Sys.time()
  dcatsResT = dcats_fit(countNor, countMut, true.conf[,order])  # true similarity matrix
  time[3] = Sys.time() - t3start
  dcatsResI = dcats_fit(countNor, countMut, diag(3))  # identity matrix
  dcatsResU = dcats_fit(countNor, countMut, get_similarity_mat(3, 0.05)) # uniform matrix
  graphs = integratedSamples@graphs$integrated_snn
  labels = Idents(integratedSamples)
  knn_mat = KNN_transition(graphs, labels)
  order = colnames(countNor)
  dcatsResK = dcats_fit(countNor, countMut, knn_mat[order, order])
  #print(dcatsRes)
  dcatsPT = data.frame(cluster = rownames(dcatsResT), dcats_pvalsT = dcatsResT$pvals)
  dcatsPI = data.frame(cluster = rownames(dcatsResI), dcats_pvalsI = dcatsResI$pvals)
  dcatsPU = data.frame(cluster = rownames(dcatsResU), dcats_pvalsU = dcatsResU$pvals)
  dcatsPK = data.frame(cluster = rownames(dcatsResK), dcats_pvalsK = dcatsResK$pvals)
  
  ## diffcyt
  # Function to create random data (one sample)
  d_random <- function(n = 20*batch_size, mean = 0, sd = 1, ncol = 20, cofactor = 5) {
    d <- sinh(matrix(rnorm(n, mean, sd), ncol = ncol)) * cofactor
    colnames(d) <- paste0("marker", sprintf("%02d", 1:ncol))
    d
  }
  # Create random data (without differential signal)
  set.seed(123)
  d_input <- list(
    Nor1 = d_random(), 
    Nor2 = d_random(), 
    Nor3 = d_random(), 
    Mut1 = d_random(),
    Mut2 = d_random(),
    Mut3 = d_random()
  )
  experiment_info <- data.frame(
    sample_id = factor(c("Nor1", "Nor2", "Nor3", "Mut1", "Mut2", "Mut3")), 
    group_id = factor(c(rep("Normal", 3), rep("Mutate", 3))), 
    stringsAsFactors = FALSE
  )
  marker_info <- data.frame(
    channel_name = paste0("channel", sprintf("%03d", 1:20)), 
    marker_name = paste0("marker", sprintf("%02d", 1:20)), 
    marker_class = factor(c(rep("type", 10), rep("state", 10)), 
                          levels = c("type", "state", "none")), 
    stringsAsFactors = FALSE
  )
  # Prepare data
  d_se <- prepareData(d_input, experiment_info, marker_info)
  # Transform data
  d_se <- transformData(d_se)
  # Generate clusters
  d_se <- generateClusters(d_se)
  # Replace with our simulated data info
  d_se@elementMetadata$sample_id = integratedSamples$batch
  d_se@elementMetadata$group_id = integratedSamples$condition
  d_se@elementMetadata$cluster_id = integratedSamples@active.ident
  # Calculate counts
  d_counts <- calcCounts(d_se)
  # Create design matrix
  design <- createDesignMatrix(experiment_info, cols_design = "group_id")
  # Create contrast matrix
  contrast <- createContrast(c(0, 1))
  # Test for differential abundance (DA) of clusters
  t4start = Sys.time()
  res_DA <- testDA_edgeR(d_counts, design, contrast)
  time[4] = Sys.time() - t4start
  diffcytP = topTable(res_DA, format_vals = TRUE) %>% 
    as.data.frame() %>% 
    rownames_to_column("cluster") %>% 
    rename(diffcyt_pvals = p_adj) %>% 
    select(cluster, diffcyt_pvals)
  
  ## results
  Res_df = merge(dcatsPT, speckleP, by = "cluster") %>% 
    merge(dcatsPI, by = "cluster") %>% 
    merge(dcatsPU, by = "cluster") %>%
    merge(dcatsPK, by = "cluster") %>%
    merge(fisherP, by = "cluster") %>% 
    merge(diffcytP, by = "cluster")
  time_df = data.frame(methods = c("fisher", "sepckle", "dcats", "diffcyt"), time = time)
  return(list(Res_df = Res_df, time_df = time_df))
}

# function 6: function to get p-values and time(test different uniform matrix)
getPandTimeTestU = function(integratedSamples, sim_list){
  dfRes = data.frame(clusterRes = integratedSamples@active.ident, batch = integratedSamples$batch, condition = integratedSamples$condition) %>% 
    tibble::rownames_to_column("cellID")
  condition = integratedSamples@meta.data$condition
  condNor<-Idents(integratedSamples)[condition == "Normal"]
  condMut<-Idents(integratedSamples)[condition == "Mutate"]
  countNor = table(sim_list$batchNor, condNor)
  countMut = table(sim_list$batchMut, condMut)
  dcatsResU1 = dcats_fit(countNor, countMut, get_similarity_mat(6, 0.05)) # uniform matrix
  dcatsResU2 = dcats_fit(countNor, countMut, get_similarity_mat(6, 0.1)) # uniform matrix
  dcatsResU3 = dcats_fit(countNor, countMut, get_similarity_mat(6, 0.15)) # uniform matrix
  #print(dcatsRes)
  dcatsPU1 = data.frame(cluster = rownames(dcatsResU1), dcats_pvalsU1 = dcatsResU1$pvals)
  dcatsPU2 = data.frame(cluster = rownames(dcatsResU2), dcats_pvalsU2 = dcatsResU2$pvals)
  dcatsPU3 = data.frame(cluster = rownames(dcatsResU3), dcats_pvalsU3 = dcatsResU3$pvals)
  
  ## results
  Res_df = merge(dcatsPU1, dcatsPU2, by = "cluster") %>% 
    merge(dcatsPU3, by = "cluster") 
  return(Res_df)
}

# 0.1, 0.2 keep the same, 0.1->0.05, 0.2->0.1, 0.2->0.25, 0.2->0.3

if(FALSE)
{
  ## test
  set.seed(123)
  probNor = c(0.1, 0.2, 0.1, 0.2, 0.2, 0.2)
  probMut = c(0.1, 0.2, 0.05, 0.1, 0.25, 0.3)
  setresolu = 0.5
  batch_size = 1000
  de_prob = c(0.5,0.5,0.1, 0.1, 0.05, 0.05)
  
  sim_list = simualtion(probNor, probMut, de_prob, batch_size)
  integratedSamples = runSeurat(sim_list, batch_size, setresolu)
  conf.mat<-table(Idents(integratedSamples), sim_list$origLabels)
  print(conf.mat)
  plot = DimPlot(integratedSamples, ncol = 3, reduction = "umap", split.by = "batch") + ggtitle("Cluster Results of Different Clusters in Different Samples")
  print(plot)
  Res = getPandTimeFSDD(integratedSamples, sim_list)
  print(Res)
}