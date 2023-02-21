library(ggplot2)
library(Seurat)
library(tidyverse)
library(magrittr)
library(viridis)
library(dplyr)
library(MAST)


#Differential Expression Analysis
DE <- FindMarkers(seurat, ident.1 = 'cluster', assay = 'RNA', test.use = "MAST", logfc.threshold = 0.05, min.pct =0.05) %>% mutate(genes = row.names(.))

de_tne$diffexpressed <- 'NO'
de_tne$diffexpressed[de_tne$avg_log2FC > 1 & -log10(de_tne$p_val) > -log10(10e-10)] <- "LOGandP"
de_tne$diffexpressed[de_tne$avg_log2FC < -1 & -log10(de_tne$p_val) > -log10(10e-10)] <- "LOGandP"
de_tne$diffexpressed[de_tne$avg_log2FC < 1 & de_tne$avg_log2FC > -1 & -log10(de_tne$p_val) > -log10(10e-10)] <- "P"
de_tne$diffexpressed[de_tne$avg_log2FC < -1 & -log10(de_tne$p_val) < -log10(10e-10)] <- "LOG"
de_tne$diffexpressed[de_tne$avg_log2FC > 1 & -log10(de_tne$p_val) < -log10(10e-10)] <- "LOG"

#Subcluster identification
seurat <- FindSubCluster(seurat, "1", graph.name = 'integrated_nn',subcluster.name = "subcluster_X",  resolution = 0.25)

#UMAP Projection
DefaultAssay(seurat1) <- 'integrated'
DefaultAssay(seurat2) <- 'integrated'

seurat.anchors <- FindTransferAnchors(reference = seurat1, query = seurat2, dims = 1:n, reference.reduction = 'pca')
predictions <- TransferData(anchorset = seurat.anchors, refdata = seurat1$MATCHTOTHIS, dims = 1:n)
seurat2 <- AddMetaData(seurat2, metadata = predictions)
seurat1 <- RunUMAP(seurat1, dims = 1:n, reduction = "pca", return.model = TRUE)
seurat2 <- MapQuery(anchorset = seurat.anchors, reference = seurat1, query = seurat2, refdata = list(celltype = "MATCHTOTHIS"), 
               reference.reduction = "pca", reduction.model = "umap")
