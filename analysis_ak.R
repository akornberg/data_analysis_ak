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

#Boxplot generation
sample_data <- function(x) {
  df <- as.data.frame(x)
  
  df <- df %>% mutate(sample_celltype = paste(df$sample, df$celltype),
                      pt_celltype_cluster =  paste(.$sample, .$celltype , .$seurat_clusters),
                      sample_cluster = paste(df$sample, df$seurat_clusters))
  
  cluster_count <- df %>% group_by(seurat_clusters) %>% summarise(cluster_count = n())
  celltype_count <- df %>% group_by(celltype, seurat_clusters) %>% summarise(celltype_count = n())
  sample_count <- df %>% group_by(sample, condition) %>% summarise(total_sample_count = n())
  sample_cluster <- df %>% group_by(sample_cluster) %>% summarise(sample_cluster_count = n())
  
  full <- cluster_count %>% merge(celltype_count) %>% merge(sample_count) %>% 
    mutate(pt_celltype = paste(.$sample, .$celltype), cluster_celltype = paste(.$cluster_label, .$celltype), 
           pt_cluster = paste(.$sample, .$seurat_clusters), pt_celltype_cluster = paste(.$sample, .$celltype , .$seurat_clusters), 
           sample_cluster = paste(.$sample, .$celltype))
  
  patient_celltype <- df %>% group_by(sample_celltype) %>% summarise(total_patient_celltype=n()) 
  full <- full_join(full, patient_celltype, by = c("pt_celltype" = "sample_celltype")) %>% replace(is.na(.), 0)
  pt_celltype_cluster <- df %>% group_by(pt_celltype_cluster) %>% summarise(pt_celltype_cluster_count=n())
  full <- full_join(full, pt_celltype_cluster, by = "pt_celltype_cluster") %>%
    mutate(pt_ct_percentage = 100*(full$total_patient_celltype/full$total_sample_count),
           pt_cl_ct_percentage = 100*(.$pt_celltype_cluster_count/.$total_patient_celltype)) %>% replace(is.na(.), 0) 
}
full <- sample_data(seurat@meta.data)
boxplot <- full %>% filter(celltype == c("CD4"),seurat_clusters =='Th17', total_patient_celltype > 9) 
stat.test <- compare_means(
  pt_cl_ct_percentage ~ condition, data = boxplot,
  method = 'wilcox.test')  %>%  mutate(plot = ifelse((.$p <= 0.05), .$p.signif, ifelse(.$p <= 0.1, .$p.format, NA)), plot2 = ifelse((.$p <= 0.05), .$p.signif, NA)) %>% print(.) %>% na.omit()

#Correlation
cd4_corr <- full %>% filter(celltype == "CD4", seurat_clusters == "Tfh") %>% select(sample, celltype, total_patient_celltype, pt_cl_ct_percentage, condition)
cd4_corr <- full %>% filter(celltype == "CD4", seurat_clusters == "Tregs") %>% select(sample, celltype, total_patient_celltype, pt_cl_ct_percentage, condition) %>%
  set_colnames(c("treg_sample", "treg_celltype", "treg_total_patient_celltype", "treg_pt_cl_ct_percentage", "treg_condition")) %>% cbind(cd4_corr)
cd4_corr$cond <- ifelse((cd4_corr$condition == "CONTROL"), "Control", "Celiac")

cd4_corr <- cd4_corr %>% filter(treg_total_patient_celltype > 9 & total_patient_celltype > 9, condition %in% c("ACD","PCD", "GFD")) %>%
  mutate(log_treg=log10(.$treg_pt_cl_ct_percentage), log_ab = log10(.$pt_cl_ct_percentage))

cd4_corr <- cd4_corr[!is.infinite(cd4_corr$log_treg),]
cd4_corr <- cd4_corr[!is.infinite(cd4_corr$log_ab),]
cor.test(cd4_corr$log_ab, cd4_corr$log_treg, method = 'pearson')

#UMAP Projection
DefaultAssay(seurat1) <- 'integrated'
DefaultAssay(seurat2) <- 'integrated'

seurat.anchors <- FindTransferAnchors(reference = seurat1, query = seurat2, dims = 1:n, reference.reduction = 'pca')
predictions <- TransferData(anchorset = seurat.anchors, refdata = seurat1$MATCHTOTHIS, dims = 1:n)
seurat2 <- AddMetaData(seurat2, metadata = predictions)
seurat1 <- RunUMAP(seurat1, dims = 1:n, reduction = "pca", return.model = TRUE)
seurat2 <- MapQuery(anchorset = seurat.anchors, reference = seurat1, query = seurat2, refdata = list(celltype = "MATCHTOTHIS"), 
               reference.reduction = "pca", reduction.model = "umap")
