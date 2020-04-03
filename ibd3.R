#Load in required packages
library(ggplot2)
library(Seurat)
library(dplyr)
library(ggridges)
library(sctransform)
library(reshape2)
library(tidyr)
library(data.table)
library(EnhancedVolcano)
library(MAST)

#Create Seurat Object
ibd.3 <- Read10X(data.dir = "/Users/akornberg/Documents/Columbia/Han Lab/Data Analysis/2019_12_IBD/gex/filtered_feature_bc_matrix/")

#Remove TCR genes
features.3 <- ibd.3@Dimnames[[1]]

features.3.tcr <- grep(c("^TRAV|^TRAJ|^TRAC|^TRBV|^TRBD|^TRBJ|^TRBC"), features.3, value = TRUE)

ibd.3.data <- ibd.3[!rownames(ibd.3) %in% features.3.tcr, ]

ibd3 <- CreateSeuratObject(counts = ibd.3.data, project="ibd3")

#Load in Assay/Metadata
adt.3 <- read.csv("/Users/akornberg/Documents/Columbia/Han Lab/Data Analysis/2019_12_IBD/adt/ibd3_adt.csv", row.names = 1)
hto.3 <- read.csv("/Users/akornberg/Documents/Columbia/Han Lab/Data Analysis/2019_12_IBD/hto/ibd3_hash_labelled.csv", row.names = 1)
tcr.3 <- read.csv("/Users/akornberg/Documents/Columbia/Han Lab/Data Analysis/2019_12_IBD/tcr/filtered_contig_annotations.csv")
cd4.cd8.3 <- read.csv("/Users/akornberg/Documents/Columbia/Han Lab/Data Analysis/2019_12_IBD/adt/ibd3_CD4CD8_assigned.csv", row.names = 1)

# Add in Metadata
ibd3 <- AddMetaData(object=ibd3, metadata = hto.3)
ibd3 <- AddMetaData(object=ibd3, metadata = cd4.cd8.3)

ibd3[["ADT"]]<-CreateAssayObject(counts= adt.3)
ibd3<- NormalizeData(ibd3, assay = "ADT", normalization.method = "CLR")
ibd3 <- ScaleData(ibd3, assay = "ADT")

shortADTnames.3<-row.names(ibd3[["ADT"]])
ADTnames.3<-paste0("adt_", shortADTnames.3)

## TCR Pre-processing. For this we want to:
# Remove the "-1" from the barcode
# Separate the chains to plot later
# Separate productive CDR3s (full length, start codon, nonstop, in-frame, CDR3 exists, and structure)
# Remove duplicates (to just remove one duplicate: x <- ibd.tcr[!duplicated(ibd.tcr$barcode), ])
tcr.3$barcode <- gsub("-1", "", tcr.3$barcode)

TCRproductive.3<-subset(tcr.3, productive == "True")

## examine what chains were detected among the productive barcodes (in this case only TRA and TRB which is good given that's what we sorted on)
print(table(TCRproductive.3$chain))

## seperate into TCA and TRA 
splitTCR.3<-split(TCRproductive.3, TCRproductive.3$chain, drop=TRUE)

# Im sure there is a better way to do this with lapply but I will leave it for another day. seperate the dataframes out and rename their columns to be unique
Multi.3<-splitTCR.3[["Multi"]]
colnames(Multi.3)<-paste0("Multi.",colnames(Multi.3))

TRA.3<-splitTCR.3[["TRA"]]
colnames(TRA.3)<-paste0("TRA.",colnames(TRA.3))

TRB.3<-splitTCR.3[["TRB"]]
colnames(TRB.3)<-paste0("TRB.",colnames(TRB.3))
TRB.3$cdr3b <- TRB.3$TRB.cdr3
## Remove duplicated cell IDs
toRemove.3<-table(TRB.3$TRB.barcode)>1
removeNames.3<-names(which(toRemove.3==TRUE))
TRBfiltered.3<-subset(TRB.3, !TRB.barcode %in% removeNames.3)
## add counts of TRB.cdr3
TRBfiltered.3 <- transform(TRBfiltered.3, TRB.CDR3count = ave(TRB.barcode, TRB.cdr3, FUN = length)) 
TRBfiltered.3$TRB.CDR3count<-as.numeric(as.character(TRBfiltered.3$TRB.CDR3count))
View(TRBfiltered.3)
## create separate factor for most highly expressed CDR3's
highClones<-names(which(table(TRBfiltered.3$TRB.cdr3)>10))
TRBfiltered.3$TRB.topClone<-ifelse(TRBfiltered.3$TRB.cdr3 %in% highClones, as.character(TRBfiltered.3$TRB.cdr3), NA)
## set rownames as cell barcodes
rownames(TRBfiltered.3)<-TRBfiltered.3$TRB.barcode
### add onto our seurat object
ibd3<-AddMetaData(object=ibd3, metadata = TRBfiltered.3)

# Frequency of Clones
ibd3$TRB.repeatedClones<-ifelse(ibd3$TRB.CDR3count > 1, as.numeric(as.character(ibd3$TRB.CDR3count)),NA)

CD3percentage<-data.frame(cell.ID=factor(),TRB.cdr3=factor(),totalDetectedTCR=numeric(),totalDetectedClone =numeric(),percentClone=numeric())

View(CD3sum) 

for(i in levels(as.factor(ibd3$patient))){
  workingDF<-subset(ibd3@meta.data, patient == i)
  resultsDF<-select(workingDF,cell.ID, TRB.cdr3, TRB.barcode, cdr3b)
  resultsDF$TRB.cdr3<-as.factor(as.character(resultsDF$TRB.cdr3))
  resultsDF$cdr3b<-as.factor(as.character(resultsDF$cdr3b))
  CD3countTable<- as.data.frame(table(as.factor(as.character(workingDF$TRB.cdr3))))
  CD3sum<-sum(CD3countTable$Freq)
  resultsDF$totalDetectedTCR<-CD3sum
  resultsDF$TRB.cdr3<-as.factor(as.character(resultsDF$TRB.cdr3))
  resultsDF$cdr3b<-as.factor(as.character(resultsDF$cdr3b))
  resultsDF<-merge(resultsDF,CD3countTable, by.x = "TRB.cdr3", by.y = "Var1")
  resultsDF$percentClone<-resultsDF$Freq/resultsDF$totalDetectedTCR
  CD3percentage<-rbind(CD3percentage,resultsDF)}
row.names(CD3percentage)<-CD3percentage$cell.ID
View(CD3percentage)

ibd3 <- AddMetaData(ibd3, metadata = CD3percentage)

# Calculate CDR3B per sample
ibd3$patient.sample <- substr(ibd3$sample, start = 1, stop = 4)

CD3percentage<-data.frame(cell.ID=factor(),TRB.cdr3=factor(),totalDetectedTCR=numeric(),totalDetectedClone =numeric(),percentClone=numeric())
for(i in levels(as.factor(ibd3$patient.sample))){
  workingDF<-subset(ibd3@meta.data, patient.sample == i)
  resultsDF<-select(workingDF,cell.ID, TRB.cdr3, TRB.barcode, cdr3b)
  resultsDF$TRB.cdr3<-as.factor(as.character(resultsDF$TRB.cdr3))
  resultsDF$cdr3b<-as.factor(as.character(resultsDF$cdr3b))
  CD3countTable<- as.data.frame(table(as.factor(as.character(workingDF$TRB.cdr3))))
  CD3sum<-sum(CD3countTable$Freq)
  resultsDF$totalDetectedTCR<-CD3sum
  resultsDF$TRB.cdr3<-as.factor(as.character(resultsDF$TRB.cdr3))
  resultsDF$cdr3b<-as.factor(as.character(resultsDF$cdr3b))
  resultsDF<-merge(resultsDF,CD3countTable, by.x = "TRB.cdr3", by.y = "Var1")
  resultsDF$percentClone<-resultsDF$Freq/resultsDF$totalDetectedTCR
  CD3percentage<-rbind(CD3percentage,resultsDF)}
row.names(CD3percentage)<-CD3percentage$cell.ID

tcr.by.sample.3 <- select(CD3percentage, "totalDetectedTCR", "Freq", "percentClone")
colnames(tcr.by.sample.3) <- c("TCR.sample", "freq.sample", "pC.sample")
ibd3 <- AddMetaData(ibd3, metadata = tcr.by.sample.3)

#Begin Seurat Analysis
ibd3[["percentMT"]]<-PercentageFeatureSet(ibd3, pattern = "^MT-")

plot1 <- FeatureScatter(ibd3, feature1 = "nCount_RNA", feature2 = "percentMT") + geom_hline(yintercept = 20)
plot2 <- FeatureScatter(ibd3, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_hline(yintercept = 2500)
CombinePlots(plots = list(plot1, plot2))

ibd3 <- subset(ibd3, subset = nFeature_RNA > 5 & nFeature_RNA < 2500 & percentMT < 20)
ibd3$experiment <- "ibd3"

## END HERE ##
ibd3 <- SCTransform(ibd3, vars.to.regress = c("percentMT"),  verbose = FALSE)
ibd3 <- RunPCA(ibd3, verbose = FALSE)
ElbowPlot(ibd3, ndims= 50)
ibd3 <- RunUMAP(ibd3, dims = 1:25, verbose = FALSE)
ibd3 <- FindNeighbors(ibd3, dims = 1:25, verbose = FALSE)
ibd3 <- FindClusters(ibd3, verbose = FALSE)
DimPlot(ibd3, label = TRUE, pt.size = 0.75)

DimPlot(ibd3, group.by = "label", pt.size = 0.75, na.value = "grey90") + 
  ggtitle("CD4/CD8") + theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values=c("darkgreen", "indianred3", "blue", "gray90")) 

DimPlot(ibd3, group.by = "stim.status", pt.size = 0.75, na.value = "grey90") + 
  ggtitle("Stim vs Nonstim") + theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values=c("gray90", "navy", "yellowgreen")) 

DimPlot(ibd3, group.by = "activity", split.by = "label", pt.size =0.75 , na.value = "grey90") + 
  ggtitle("Inflammation") + theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values=c("firebrick", "steelblue2", "steelblue2", "gray90")) 

DimPlot(ibd3, group.by = "experiment" ,pt.size = 0.75, na.value = "grey90") + 
  ggtitle("Experiment") + theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values=c("darkgoldenrod2", "navy")) 

DimPlot(ibd3, group.by = "treatment" ,pt.size = 0.75, na.value = "grey90") + 
  ggtitle("Treatment") + theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values=c("olivedrab3", "indianred", "gray60", "navy", "plum", "firebrick", "darkorange")) 

ibd3$pCadj <- ifelse(ibd3$totalDetectedTCR < 50, NA, ibd3$percentClone)

FeaturePlot(ibd3, features = "Freq", cols = c("lightgoldenrod2", "firebrick")) + 
  ggtitle("TCRB Frequency")

FeaturePlot(ibd3, features = c("GZMB", "PRF1"), ncol=2)

