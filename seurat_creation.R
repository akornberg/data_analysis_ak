# Data Input

#Required Packages
library(Seurat)
library(tidyverse)
library(ggplot2)
library(magrittr)

#Required data inputs
counts <- Read10X(data.dir = "/PATH/TO/filtered_feature_bc_matrix/")
hto <- read.csv(file="/PATH/TO/ASSIGNED/HASHES.csv", row.names = 1)
metadata <- read.csv(file="/Users/adamkornberg/Documents/Columbia/Han Lab/Data Analysis/celiac/2022_12/metadata.csv")
adt_targets <- read.csv(file="/PATH/TO/ASSIGNED/CITE_SEQ_NAMES.csv")
adt_isotype <- read.csv(file="/PATH/TO/ASSIGNED/ISOTYPE_CONTROL_NAMES.csv", row.names = 1)
tcrab <- read.csv("/PATH/TO/filtered_ab_contig_annotations.csv")
tcrgd <- read.csv("/PATH/TO/filtered_gd_contig_annotations.csv")
experiment <- "_BATCH_X"

# Seurat object creation for combined RNA/HTO/ADT assays 
full_seurat_creation <- function(x, y) {
  
  seurat <- CreateSeuratObject(counts = counts$`Gene Expression`, project = y)
  
  adt_orig <- as.data.frame(counts$`Antibody Capture`)
  row.names(adt_orig) <- adt_targets$x
  seurat[["ADT"]] <- CreateAssayObject(counts = adt_orig)
  
  adt_counts <- (as.data.frame(seurat@assays$ADT@counts)) %>% 
    mutate(adt = row.names(.)) %>% 
    left_join(adt_isotype, by = c("adt"="names")) %>% 
    set_rownames(.$adt) 
  
  adt <-data.frame(row.names = row.names(seurat@meta.data))    
  for(i in levels(as.factor(adt_counts$isotype))){
    df <- subset(adt_counts, isotype == i)
    df <- select(df, 1:(ncol(df)-2))%>% t() %>% as.data.frame()
    ic <- df[[i]]
    finaldf <- (df - ic)
    adt <- cbind(adt, finaldf)
  }
  adt[adt<0] <- 0
  adt <- as.data.frame(t(adt))
  adt <- arrange(adt, row.names(adt))
  seurat[["CITE"]] <- CreateAssayObject(counts = adt)
  
  hto <- hto %>% 
    mutate(barcode = row.names(hto)) %>%
    left_join(metadata, by="sample") %>%
    set_rownames(.$barcode) 
  
  seurat@meta.data <- seurat@meta.data %>% 
    mutate(old.barcode = row.names(.), adj.barcode =  gsub("-1", experiment, old.barcode)) %>% 
    set_rownames(.$adj.barcode)
  
  row.names(hto) <- paste0(row.names(hto), experiment)
  seurat <- AddMetaData(object=seurat, metadata = hto)
  seurat <- RenameCells(seurat, new.names = seurat$adj.barcode)
  
  seurat <- PercentageFeatureSet(seurat, pattern = "^MT-", col.name = "percent.mt")
  seurat <- PercentageFeatureSet(seurat, pattern = "^RPS", col.name = "percent.rps")
  seurat <- PercentageFeatureSet(seurat, pattern = "^RPL", col.name = "percent.rpl")
  seurat@meta.data$percent.rib <- seurat@meta.data$percent.rps+seurat@meta.data$percent.rpl
  
  seurat@meta.data$exclude <- ifelse((seurat@meta.data$percent.rpl < 2), "low_quality",
                                     ifelse(seurat@meta.data$percent.mt < 10, "quality", 
                                            ifelse(seurat@meta.data$percent.mt > 15, "low_quality", 
                                                   ifelse(seurat@meta.data$percent.rpl < 2.5, "low_quality", "quality"))))
  seurat
} 

#TCRAB Input (HUMAN)
tcrab_input <- function(x) {
  tcrab <- as.data.frame(x)
  tcrab <- tcrab %>% mutate(adj.tcr.barcode =  gsub("-1", experiment, barcode)) 
  
  tcra <- tcrab %>% filter(productive == "True", chain == "TRA") %>% 
    mutate(duplicated = duplicated(barcode)) %>% 
    filter(duplicated == FALSE) %>%
    select(adj.tcr.barcode, length, v_gene, d_gene, j_gene, c_gene, cdr3, cdr3_nt, reads, umis, raw_clonotype_id) %>%
    set_colnames(c("adj.tcr.barcode", "length_alpha", "v_gene_alpha", "d_gene_alpha", "j_gene_alpha", "c_gene_alpha", "cdr3a", "cdr3a_nt", "reads_alpha", "umis_alpha", "raw_alpha_clonotype_id"))
  
  tcrb <- tcrab %>% filter(productive == "True", chain == "TRB") %>% 
    mutate(duplicated = duplicated(barcode)) %>% 
    filter(duplicated == FALSE) %>% 
    select(adj.tcr.barcode, length, v_gene, d_gene, j_gene, c_gene, cdr3, cdr3_nt, reads, umis, raw_clonotype_id) %>%
    set_colnames(c("adj.tcr.barcode", "length_beta", "v_gene_beta", "d_gene_beta", "j_gene_beta", "c_gene_beta", "cdr3b", "cdr3b_nt", "reads_beta", "umis_beta", "raw_beta_clonotype_id"))
  
  tcrab.base <- full_join(tcra,tcrb, by="adj.tcr.barcode") %>% set_rownames(.$adj.tcr.barcode)
}

#TCRGD Input (HUMAN)
tcrgd_input <- function(x) {
  tcrgd <- as.data.frame(x)
  tcrgd <- tcrgd %>% mutate(adj.tcr.barcode =  gsub("-1", experiment, barcode)) 
  
  tcrg <- tcrgd %>% filter(chain == "TRG", productive == "True") %>% 
    mutate(duplicated = duplicated(barcode)) %>% 
    filter(duplicated == FALSE) %>% 
    select(adj.tcr.barcode, length, v_gene, d_gene, j_gene, c_gene, cdr3, cdr3_nt, reads, umis, raw_clonotype_id) %>%
    set_colnames(c("adj.tcr.barcode", "length_gamma", "v_gene_gamma", "d_gene_gamma", "j_gene_gamma", "c_gene_gamma", "cdr3g", "cdr3g_nt", "reads_gamma", "umis_gamma", "raw_gamma_clonotype_id"))
  
  tcrd <- tcrgd  %>% filter(chain == "TRD", productive == "True") %>% 
    mutate(duplicated = duplicated(barcode)) %>% 
    filter(duplicated == FALSE) %>% 
    select(adj.tcr.barcode, length, v_gene, d_gene, j_gene, c_gene, cdr3, cdr3_nt, reads, umis, raw_clonotype_id) %>%
    set_colnames(c("adj.tcr.barcode", "length_delta", "v_gene_delta", "d_gene_delta", "j_gene_delta", "c_gene_delta", "cdr3d", "cdr3d_nt", "reads_delta", "umis_delta", "raw_delta_clonotype_id"))
  
  tcrgd.base <- full_join(tcrg,tcrd, by="adj.tcr.barcode") %>% set_rownames(.$adj.tcr.barcode)
}

### TCRAB + TCRGD CALCULATE ###
full_tcr_calculate <- function(x) {
  df <- as.data.frame(x)
  df <- df %>% 
    mutate(cdr3b.sample = paste(.$cdr3b, .$sample), cdr3a.sample = paste(.$cdr3a, .$sample), cdr3b.pt = paste(.$cdr3b, .$patient), cdr3a.pt = paste(.$cdr3a, .$patient),
           cdr3g.sample = paste(.$cdr3g, .$sample), cdr3d.sample = paste(.$cdr3d, .$sample), cdr3g.pt = paste(.$cdr3g, .$patient), cdr3d.pt = paste(.$cdr3d, .$patient),
           id = 1:nrow(x)) %>% set_rownames(.$adj.barcode)
  
  df <- df %>%
    group_by(cdr3b, cdr3b.sample) %>% summarise(cdr3b.sample.count = n()) %>% ungroup() %>% na.omit() %>% select(!cdr3b) %>%
    full_join(df, by = "cdr3b.sample") %>% arrange(id)
  
  df <- df %>%
    group_by(cdr3b, cdr3b.pt) %>% summarise(cdr3b.pt.count = n()) %>% ungroup() %>% na.omit() %>% select(!cdr3b) %>%
    full_join(df, by = "cdr3b.pt") %>% arrange(id)
  
  df <- df %>% 
    mutate(cdr3b.sample.count = df$cdr3b.sample.count, cdr3b.pt.count = df$cdr3b.pt.count) %>% set_rownames(.$adj.barcode) 
  
  tcrs <- df %>% filter(cdr3b.sample.count > 0)
  tcr.by.sample <- tcrs %>% group_by(sample) %>% summarise(sample.tcr.count = n())
  tcr.by.pt <- tcrs %>% group_by(patient) %>% summarise(pt.tcr.count = n())
  
  df <- df %>% full_join(tcr.by.sample, by = "sample") %>%
    mutate(cdr3b.freq = 100*(cdr3b.sample.count/sample.tcr.count)) %>% set_rownames(.$adj.barcode)
  
  df <- df %>% full_join(tcr.by.pt, by = "patient") %>%
    mutate(cdr3b.pt.freq = 100*(cdr3b.pt.count/pt.tcr.count)) %>% set_rownames(.$adj.barcode)
  
  df$cdr3b.sample.count.cut <- cut(df$cdr3b.sample.count, breaks = c(0,1,5,25,100,1000), labels = c("one", "2-5", "5-25", "25-100", "100+"))
  df$cdr3b.freq.adj <- ifelse((df$sample.tcr.count < 50), NA, df$cdr3b.freq)
  
  df <- df %>%
    group_by(cdr3d, cdr3d.sample) %>% summarise(cdr3d.sample.count = n()) %>%  ungroup() %>% na.omit() %>% select(!cdr3d) %>%
    full_join(df, by = "cdr3d.sample") %>% arrange(id)
  
  df <- df %>%
    group_by(cdr3g, cdr3g.sample) %>% summarise(cdr3g.sample.count = n()) %>% ungroup() %>% na.omit() %>% select(!cdr3g) %>%
    full_join(df, by = "cdr3g.sample") %>% arrange(id)
  
  df <- df %>% 
    mutate(cdr3d.sample.count = df$cdr3d.sample.count, cdr3g.sample.count = df$cdr3g.sample.count) %>%  set_rownames(.$adj.barcode)
  
  tcrs <- df %>% filter(cdr3d.sample.count > 0)
  tcr.by.sample <- tcrs %>% group_by(sample) %>% summarise(sample.tcrd.count = n())
  tcr.by.pt <- tcrs %>% group_by(patient) %>% summarise(pt.tcrgd.count = n())
  
  df <- df %>% full_join(tcr.by.sample, by = "sample") %>% mutate(cdr3d.freq = 100*(cdr3d.sample.count/sample.tcrd.count)) %>% set_rownames(.$adj.barcode)
  df$cdr3d.sample.count.cut <- cut(df$cdr3d.sample.count, breaks = c(0,1,5,25,100,1000), labels = c("one", "2-5", "5-25", "25-100", "100+"))
  df$cdr3d.freq.adj <- ifelse((df$sample.tcrd.count < 20), NA, df$cdr3d.freq)
  
  tcrs <- df %>% filter(cdr3g.sample.count > 0)
  tcr.by.sample <- tcrs %>% group_by(sample) %>% summarise(sample.tcrg.count = n())
  tcr.by.pt <- tcrs %>% group_by(patient) %>% summarise(pt.tcrg.count = n())
  
  df <- df %>% full_join(tcr.by.sample, by = "sample") %>% mutate(cdr3g.freq = 100*(cdr3g.sample.count/sample.tcrg.count)) %>% set_rownames(.$adj.barcode)
  df$cdr3g.sample.count.cut <- cut(df$cdr3g.sample.count, breaks = c(0,1,5,25,100,1000), labels = c("one", "2-5", "5-25", "25-100", "100+"))
  df$cdr3g.freq.adj <- ifelse((df$sample.tcrg.count < 20), NA, df$cdr3g.freq)
  
  df$tcr_cut <- ifelse((is.na(df$cdr3b)), "No TCR", 
                       ifelse(df$cdr3b.sample.count == 1, "No Expansion",
                              ifelse(df$cdr3b.freq <= 0.5, "<0.5%",
                                     ifelse(df$cdr3b.freq <= 2, "0.5-2%",
                                            ifelse(df$cdr3b.freq <= 5, "2-5%",
                                                   ifelse(df$cdr3b.freq <= 100, "5%+", "other"))))))
  
  test <- df %>% filter(cdr3b.sample.count > 1)
  
  df$tcr_cut2 <- ifelse((is.na(df$cdr3b)), "No TCR", 
                        ifelse(df$cdr3b.sample.count == 1, "No Expansion",
                               ifelse(df$cdr3b.freq <= summary(test$cdr3b.freq)[[2]], "Q1",
                                      ifelse(df$cdr3b.freq <= summary(test$cdr3b.freq)[[3]], "Q2",
                                             ifelse(df$cdr3b.freq <= summary(test$cdr3b.freq)[[5]], "Q3",
                                                    ifelse(df$cdr3b.freq <= 100, "Q4", "other"))))))
  
  test <- df %>% filter(cdr3d.sample.count > 1)
  
  df$tcrd_cut <- ifelse((is.na(df$cdr3d)), "No TCR", 
                        ifelse(df$cdr3d.sample.count == 1, "No Expansion",
                               ifelse(df$cdr3d.freq <= 0.5, "<0.5%",
                                      ifelse(df$cdr3d.freq <= 2, "0.5-2%",
                                             ifelse(df$cdr3d.freq <= 5, "2-5%",
                                                    ifelse(df$cdr3d.freq <= 100, "5%+", "other"))))))
  
  
  df$tcrd_cut2 <- ifelse((is.na(df$cdr3d)), "No TCR", 
                         ifelse(df$cdr3d.sample.count == 1, "No Expansion",
                                ifelse(df$cdr3d.freq <= summary(test$cdr3d.freq)[[2]], "Q1",
                                       ifelse(df$cdr3d.freq <= summary(test$cdr3d.freq)[[3]], "Q2",
                                              ifelse(df$cdr3d.freq <= summary(test$cdr3d.freq)[[5]], "Q3",
                                                     ifelse(df$cdr3d.freq <= 100, "Q4", "other"))))))
  
  test <- df %>% filter(cdr3d.sample.count > 1)
  
  df$tcrg_cut <- ifelse((is.na(df$cdr3g)), "No TCR", 
                        ifelse(df$cdr3g.sample.count == 1, "No Expansion",
                               ifelse(df$cdr3g.freq <= 0.5, "<0.5%",
                                      ifelse(df$cdr3g.freq <= 2, "0.5-2%",
                                             ifelse(df$cdr3g.freq <= 5, "2-5%",
                                                    ifelse(df$cdr3g.freq <= 100, "5%+", "other"))))))
  
  df$tcrg_cut2 <- ifelse((is.na(df$cdr3g)), "No TCR", 
                         ifelse(df$cdr3g.sample.count == 1, "No Expansion",
                                ifelse(df$cdr3g.freq <= summary(test$cdr3g.freq)[[2]], "Q1",
                                       ifelse(df$cdr3g.freq <= summary(test$cdr3g.freq)[[3]], "Q2",
                                              ifelse(df$cdr3g.freq <= summary(test$cdr3g.freq)[[5]], "Q3",
                                                     ifelse(df$cdr3g.freq <= 100, "Q4", "other"))))))
  df <- as.data.frame(df) %>% set_rownames(.$adj.barcode)
  
}

#Example data setup
seurat <- full_seurat_creation(counts, "Batch X")
tcrab.base <- tcrab_input(tcrab)
tcrgd.base <- tcrgd_input(tcrgd)
seurat <- AddMetaData(object=seurat, metadata = tcrab.base)
seurat <- AddMetaData(object=seurat, metadata = tcrgd.base)
seurat@meta.data <- full_tcr_calculate(seurat@meta.data)
