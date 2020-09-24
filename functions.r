if(FALSE){
  library(splatter)
  library(Seurat)
  library(speckle)
  library(DCATS)
  library(ggplot2)
  library(tidyverse)
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
    plyr::mapvalues(from = c("1", "0", "2"), to = c("A", "B", "C"))
  
  return(integratedSamples)
}

getPandTime = function(integratedSamples, sim_list){
  time = rep(NA,3)
  dfRes = data.frame(clusterRes = integratedSamples@active.ident, batch = integratedSamples$batch, condition = integratedSamples$condition) %>% 
    tibble::rownames_to_column("cellID")
  ## Fisher's exact test
  dfCount = dfRes %>% 
    group_by(condition, clusterRes) %>% 
    summarise(n = n()) %>% 
    pivot_wider(names_from = "clusterRes", values_from = "n") %>% 
    mutate(nonA = B + C,
           nonB = A + C,
           nonC = A + B) %>% 
    select(A, B, C:nonC)
  t1start = Sys.time()
  fisher_pvals = rep(NA,3)
  for (i in 1:3){
    fisher_pvals[i] = fisher.test(dfCount[,c(i+1,i+4)])$p.value
  }
  time[1] = Sys.time() - t1start
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
  #print(true.conf)
  condition = integratedSamples@meta.data$condition
  condNor<-Idents(integratedSamples)[condition == "Normal"]
  condMut<-Idents(integratedSamples)[condition == "Mutate"]
  countNor = table(sim_list$batchNor, relevel(condNor, "A"))
  countMut = table(sim_list$batchMut, relevel(condMut, "A"))
  t3start = Sys.time()
  dcatsRes1 = dcats_fit(countNor, countMut, true.conf)  # true similarity matrix
  time[3] = Sys.time() - t3start
  dcatsRes2 = dcats_fit(countNor, countMut, diag(3))  # identity matrix
  #print(dcatsRes)
  dcatsP1 = data.frame(cluster = rownames(dcatsRes1), dcats_pvals1 = dcatsRes1$pvals)
  Res_df = merge(dcatsP1, speckleP, by = "cluster") %>% 
    mutate(dcats_pvals2 = dcatsRes2$pvals,
           fisher_pvals = fisher_pvals)
  time_df = data.frame(methods = c("fisher", "sepckle", "dcats"), time = time)
  return(list(Res_df = Res_df, time_df = time_df))
}


if(FALSE)
{
  ## test
  set.seed(123)
  probNor = c(1/3,1/3,1/3)
  probMut = c(1/2,1/6,1/3)
  setresolu = 0.5
  batch_size = 600
  de_prob = c(0.06,0.06,0.5)
  
  sim_list = simualtion(probNor, probMut, de_prob, batch_size)
  integratedSamples = runSeurat(sim_list, batch_size, setresolu)
  Res = getPandTime(integratedSamples, sim_list)
  print(Res)
}
