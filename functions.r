library(splatter)
library(Seurat)
library(speckle)
library(DCATS)
library(ggplot2)
library(tidyverse)

probNor = c(1/3,1/3,1/3)
probMut = c(1/2,1/6,1/3)
setresolu = 0.5
de_prob = c(0.1,0.1,0.5)
batch_size = 600

# fucntion 1: simulation
simualtion = function(probNor, probMut, resolution, de_prob, batch_size){
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

## test function
probNor = c(1/3,1/3,1/3)
probMut = c(1/2,1/6,1/3)
setresolu = 0.5
de_prob = c(0.1,0.1,0.5)
batch_size = 600

sim_list = simualtion(probNor, probMut, resolution, de_prob, batch_size)
head(sim_list$simNor_mat)
head(sim_list$simMut_mat)
print(sim_list$batchNor)
print(sim_list$batchMut)

# function 2: Seurat process(might need few minutes)
runSeurat = function(sim_list, batch_size){
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
  
  return(integratedSamples)
}

## test function
probNor = c(1/3,1/3,1/3)
probMut = c(1/2,1/6,1/3)
setresolu = 0.5
de_prob = c(0.1,0.1,0.5)
batch_size = 600

sim_list = simualtion(probNor, probMut, resolution, de_prob, batch_size)
integratedSamples = runSeurat(sim_list, batch_size)