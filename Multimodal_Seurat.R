#Load in required packages
library(ggplot2)
library(Seurat)
library(dplyr)
library(ggridges)
library(sctransform)
library(reshape2)
library(tidyr)
library(data.table)

#Create Seurat Object
ibd.data <- Read10X(data.dir = "/Users/akornberg/Documents/Columbia/Han Lab/Data Analysis/2019_08_IBD/GEX/filtered_feature_bc_matrix/")
ibd <- CreateSeuratObject(counts = ibd.data, project = "ibd_IBD_2019_08")

#Load in Assay/Metadata
adt <- read.csv("/Users/akornberg/Documents/Columbia/Han Lab/Data Analysis/2019_08_IBD/AdT/2019_08_IBD_AdT.csv", row.names = 1)
hto <- read.csv("/Users/akornberg/Documents/Columbia/Han Lab/Data Analysis/2019_08_IBD/HTO/2019_08_IBD_hash_general.csv", row.names = 1)
tcr <- read.csv("/Users/akornberg/Documents/Columbia/Han Lab/Data Analysis/2019_08_IBD/TCR/filtered_contig_annotations_IBD_Aug_2019.csv")


#Begin Seurat Analysis
ibd[["percentMT"]]<-PercentageFeatureSet(ibd, pattern = "^MT-")
ibd <- SCTransform(ibd, vars.to.regress = "percentMT", verbose = FALSE)

## EDA plot
VlnPlot(ibd, features = c("nFeature_RNA", "nCount_RNA", "percentMT"), ncol = 3)
##
plot1 <- FeatureScatter(ibd, feature1 = "nCount_RNA", feature2 = "percentMT")
plot2 <- FeatureScatter(ibd, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))

## MV plot
top10<-head(VariableFeatures(ibd),10)
## plot variable features w/top 10 most variable genes labeled
plot1 <- VariableFeaturePlot(ibd)
plot2 <- LabelPoints(plot = plot1, points = top10)
CombinePlots(plots = list(plot1, plot2))

#Run PCA + Choose PC's to include in the analysis
ibd <- RunPCA(ibd, verbose = FALSE)

ElbowPlot(ibd, ndims= 50)
ibd <- RunUMAP(ibd, dims = 1:30, verbose = FALSE)
ibd <- FindNeighbors(ibd, dims = 1:30, verbose = FALSE)
ibd <- FindClusters(ibd, verbose = FALSE)
DimPlot(ibd, label = TRUE) + NoLegend()

#Identify markers of each cluster
clusterMarkers <- FindAllMarkers(ibd, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
topMarkers<-clusterMarkers %>% group_by(cluster) %>% top_n(n = 25, wt = avg_logFC)

shortMarkers<-clusterMarkers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
shortM<-shortMarkers[,c(6:7)]
shortM$rank<-c(1:10)
shortM %>% spread(rank, gene)

#Add in (previously loaded) AdT data as an assay object
ibd[["ADT"]]<-CreateAssayObject(counts= adt)

ibd<- NormalizeData(ibd, assay = "ADT", normalization.method = "CLR")
ibd <- ScaleData(ibd, assay = "ADT")

shortADTnames<-row.names(ibd[["ADT"]])
ADTnames<-paste0("adt_", shortADTnames)

# FEATURE PLOT - Can show AdT or GEX data. Cutoffs tend to make AdT data appear more clean
FeaturePlot(ibd, features = c("CD44","CD69","adt_CD4","adt_CD8"), ncol=2)

RidgePlot(ibd, features = c("adt_CD4", "adt_CD8"), ncol = 2)

# Add (previously loaded) HTO data as metadata
ibd<-AddMetaData(object=ibd, metadata = hto)

# Create new column with dx + activity
ibd@meta.data$dx.activity<-paste0(ibd@meta.data$xx,ibd@meta.data$activity)
ibd@meta.data$new.column <- 

## counts of how many cells from each sample we got
table(hto$sample)
## total cells - should be same as GEX barcodes (?)
sum(table(hto$sample))

# Group by sample
DimPlot(ibd, group.by = "sample", pt=1, cols = c("brown1","darkred",
"darkorange","darkorange3","darkgoldenrod1","darkgoldenrod3",
"darkolivegreen1","darkolivegreen4","chartreuse","chartreuse3",
"blue","blue4","cadetblue1","cadetblue4","cyan", "cyan4",
"darkorchid1","darkorchid4", "gray90", "gray90") ,na.value = "grey90")

# Group by condition
DimPlot(ibd, group.by = "activity", pt=1, cells.highlight = WhichCells(ibd, ibd$activity=="active"))
DimPlot(ibd, group.by = "stim.status", pt=1)

# Group by patient
DimPlot(ibd, group.by = "stim.status", pt=1, split.by = "patient")
DimPlot(ibd, group.by = "activity", pt=1, split.by = "patient")

# Group by region and disease
DimPlot(ibd, group.by = "region", pt=1)
DimPlot(ibd, group.by = "dx", pt=1)

# Activity by region and disease
DimPlot(ibd, group.by = "activity", pt=1, split.by = "region")
DimPlot(ibd, group.by = "activity", pt=1, split.by = "dx")
DimPlot(ibd, group.by = "stim.status", pt=1, split.by = "dx")

# Split by underlying disease and activity
DimPlot(ibd, group.by = "stim.status", pt=1, split.by = "dx.activity")

# Activity by stim
DimPlot(ibd, group.by = "stim.status", pt=1, split.by = "activity")

## TCR Pre-processing. For this we want to:
  # Remove the "-1" from the barcode
  # Separate the chains to plot later
  # Separate productive CDR3s (full length, start codon, nonstop, in-frame, CDR3 exists, and structure)
  # Remove duplicates (to just remove one duplicate: x <- ibd.tcr[!duplicated(ibd.tcr$barcode), ])
tcr$barcode <- gsub("-1", "", tcr$barcode)

TCRproductive<-subset(tcr, productive == "True")

## examine what chains were detected among the productive barcodes (in this case only TRA and TRB which is good given that's what we sorted on)
print(table(TCRproductive$chain))

## seperate into TCA and TRA 
splitTCR<-split(TCRproductive, TCRproductive$chain, drop=TRUE)

# Im sure there is a better way to do this with lapply but I will leave it for another day. seperate the dataframes out and rename their columns to be unique
Multi<-splitTCR[["Multi"]]
colnames(Multi)<-paste0("Multi.",colnames(Multi))

TRA<-splitTCR[["TRA"]]
colnames(TRA)<-paste0("TRA.",colnames(TRA))

TRB<-splitTCR[["TRB"]]
colnames(TRB)<-paste0("TRB.",colnames(TRB))

## Remove duplicated cell IDs
toRemove<-table(TRB$TRB.barcode)>1
removeNames<-names(which(toRemove==TRUE))
TRBfiltered<-subset(TRB, !TRB.barcode %in% removeNames)

View(TRBfiltered)

## add counts of TRB.cdr3
TRBfiltered <- transform(TRBfiltered, count = ave(TRB.barcode, TRB.cdr3, FUN = length)) 
TRBfiltered$TRB.cdr3<-as.numeric(as.character(TRBfiltered$count))
## set rownames as cell barcodes
rownames(TRBfiltered)<-TRBfiltered$TRB.barcode
### add onto our seurat object
ibd<-AddMetaData(object=ibd, metadata = TRBfiltered)
### Featurelot with clonotype frequencies 
FeaturePlot(ibd, "TRB.cdr3")

View(ibd)
#ggplot example - make your plots look nice
plot <- DimPlot(ibd, group.by="dx.activity")

plot + ggtitle("X") + theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values=c("darkblue", "lightblue","firebrick","lightpink"),
  labels=c("A", "B", "C", "D"))



