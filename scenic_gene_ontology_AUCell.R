library(ggplot2)
library(Seurat)
library(tidyverse)
library(magrittr)
library(viridis)
library(dplyr)
library(SCENIC)
library(AUCell)


#SCENIC
## Get data from sce object:
exprMat <- as.matrix(seurat@assays$RNA@data)
cellInfo <- seurat@meta.data

scenicOptions <- initializeScenic(org = "hgnc", dbDir = "PATH/TO/cisTarget_databases/", datasetTitle="celiac") 

hg38_dbs <- list('500bp'= 'hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather', 
                 '10kb' = 'hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather')
db_mcVersion <- 'v9'

scenicOptions@settings$dbs <- hg38_dbs
scenicOptions@settings$db_mcVersion <- db_mcVersion

genesKept <- geneFiltering(exprMat, scenicOptions=scenicOptions, minSamples = 50, minCountsPerGene = 100)
exprMat_filtered <- exprMat[genesKept, ]
runCorrelation(exprMat_filtered, scenicOptions)
exportsForArboreto(exprMat_filtered, scenicOptions, dir = "int")

#Run Arboreto runGRNBoost2 for faster runtime when compared to runGenie3

GRNBoost_output <- read.delim("PATH/TO/int/net1_grn_output.tsv", header=FALSE)
colnames(GRNBoost_output) <- c("TF","Target","weight")
saveRDS(GRNBoost_output, file="PATH/TO/int/1.4_GENIE3_linkList.Rds")

scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
scenicOptions <- runSCENIC_2_createRegulons(scenicOptions)
scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat = exprMat)
regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
scenic <- as.data.frame(t(as.data.frame(getAUC(regulonAUC))))
colnames(scenic) <- gsub("\\s*\\([^\\)]+\\)","",as.character(colnames(scenic)))
colnames(scenic) <- gsub("extended","ext",as.character(colnames(scenic)))

# Gene Ontology AUC
cells_rankings <- AUCell_buildRankings(celiac@assays$RNA@counts)

  ##Input file downloaded from: https://data.broadinstitute.org/gsea-msigdb/msigdb/release/7.2/
c5v7 <- read.gmt("./c5.all.v7.2.symbols.gmt")

colnames(c5v7)[1] <- "ont"
geneSets <- lapply(unique(c5v7$ont), function(x){print(x);c5v7$gene[c5v7$ont == x]})
names(geneSets) <- unique(c5v7$ont)

cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.1)
