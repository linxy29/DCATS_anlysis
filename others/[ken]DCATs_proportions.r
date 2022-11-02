library(Seurat)
library(DCATS)

samples <- readRDS("/storage/holab/hcc_data/RDS_files/SCT/HCC.final.SCT_pc10_res0.1.rds")

DefaultAssay(samples) <- "SCT"

SNN_graph_all <- samples@graphs$integrated_snn

Idents(samples) <- "day_type"

for (type in c("Pos", "Neg")) {
	for (x in c("Day3", "Day10", "Day30")) {
		for (y in c("Day10", "Day30")) {
			if (x == y) {
				next
			} else if (x == "Day30" & y == "Day10") {
				next
			}
			
			dayx <- paste0(x,"_",type)
			dayy <- paste0(y,"_",type)
			
			cat(dayx,"vs",dayy,"\n")
			samples.dayx <- subset(samples, idents=dayx)
			samples.dayy <- subset(samples, idents=dayy)
			Idents(samples.dayx) <- "seurat_clusters"
			Idents(samples.dayy) <- "seurat_clusters"
			matx <- t(as.matrix(table(Idents(samples.dayx))))
			colnames(matx) <- NULL
			print(matx)
			maty <- t(as.matrix(table(Idents(samples.dayy))))
			colnames(maty) <- NULL
			print(maty)
			
			
			cells.x.y <- c()
			cells.x.y <- c(cells.x.y, colnames(samples.dayx))
			cells.x.y <- c(cells.x.y, colnames(samples.dayy))
			
			samples.dayx.dayy <- subset(samples, cells = cells.x.y)
			
			SNN_graph <- SNN_graph_all[rownames(samples.dayx.dayy@meta.data),rownames(samples.dayx.dayy@meta.data)]
			
			Idents(samples.dayx.dayy) <- "seurat_clusters"
			cell_labels <- Idents(samples.dayx.dayy)
			
			similarity_mat <- KNN_transition(SNN_graph, cell_labels)
			sorted_similarity_mat <- similarity_mat
			sorted_similarity_mat <- sorted_similarity_mat[sort(rownames(similarity_mat)),sort(rownames(similarity_mat))]

			DCATs_fit_result <- dcats_fit(matx, maty, sorted_similarity_mat)
			rownames(DCATs_fit_result) <- as.numeric(rownames(DCATs_fit_result)) - 1
			
			print(DCATs_fit_result)
			cat("\n\n")
		}
	}
}