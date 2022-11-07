#https://www.bioconductor.org/packages/release/bioc/html/splatter.html
library(splatter)
library(Seurat)
#https://github.com/huangyh09/DCATS/
library(DCATS)


param.groups <- newSplatParams(batchCells = 600, nGenes = 100)

sim1 <- splatSimulateGroups(param.groups, group.prob = c(1/3,1/3,1/3), de.prob = c(0.02,0.02,0.5), verbose = FALSE)
celllabels_orig <- sim1@colData@listData$Group
sim1count_mat <- counts(sim1)
seuratobject <- CreateSeuratObject(counts = sim1count_mat, project="Splatter")

seuratobject <- NormalizeData(seuratobject)
seuratobject <- FindVariableFeatures(seuratobject, selection.method = "vst", nfeatures = 50)
all.genes <- rownames(seuratobject)
seuratobject <- ScaleData(seuratobject,features = all.genes)
seuratobject <- RunPCA(seuratobject,features = VariableFeatures(object = seuratobject))

seuratobject<-FindNeighbors(seuratobject, k.params=3, dims = 1:5)
seuratobject<-FindClusters(seuratobject, resolution = 0.7, algorithm=2)



DimPlot(seuratobject, reduction = "pca", group.by='ident')

seuratobject <- RunUMAP(seuratobject, dims = 1:10)
DimPlot(seuratobject, reduction = "umap", group.by='ident')


Idents(seuratobject)


DimHeatmap(seuratobject, dims = 1:3, cells = 300, balanced = TRUE)


# Confusion matrix
conf.mat<-table(Idents(seuratobject), celllabels_orig)
conf.mat;

true.conf<-t(t(conf.mat)/apply(conf.mat,2,sum))

ind.cond1<-sample(1:600,300)
ind.cond2<-setdiff(1:600,ind.cond1);

cond1<-Idents(seuratobject)[ind.cond1];
cond2<-Idents(seuratobject)[ind.cond2];

cell.count.mat<-cbind(table(cond1),table(cond2))


# Fisher's exact test
fisher.test(cell.count.mat)$p.value

# DCATS
dcats_fit(table(cond1), table(cond2), diag(3), n_samples = 1)
dcats_fit(table(cond1), table(cond2), true.conf, n_samples = 1)


# other methods for comparisons:
# https://github.com/Oshlack/speckle
# https://www.bioconductor.org/packages/release/bioc/html/diffcyt.html
# Good to cehck: https://f1000research.com/articles/6-748/v4
# https://github.com/KlugerLab/DAseq
# Good to check out: https://academic-oup-com.eproxy.lib.hku.hk/gigascience/article/8/9/giz107/5572529
# https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-019-3211-9









