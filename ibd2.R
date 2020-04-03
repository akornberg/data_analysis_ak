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
ibd.2 <- Read10X(data.dir = "/Users/akornberg/Documents/Columbia/Han Lab/Data Analysis/2019_11_IBD/gex/filtered_feature_bc_matrix/")

#Remove TCR genes
features.2 <- ibd.2@Dimnames[[1]]

features.2.tcr <- grep(c("^TRAV|^TRAJ|^TRAC|^TRBV|^TRBD|^TRBJ|^TRBC"), features.2, value = TRUE)

ibd.2.data <- ibd.2[!rownames(ibd.2) %in% features.2.tcr, ]

ibd2 <- CreateSeuratObject(counts = ibd.2.data, project="ibd2")

#Load in Assay/Metadata
adt.2 <- read.csv("/Users/akornberg/Documents/Columbia/Han Lab/Data Analysis/2019_11_IBD/adt/ibd2_adt.csv", row.names = 1)
hto.2 <- read.csv("/Users/akornberg/Documents/Columbia/Han Lab/Data Analysis/2019_11_IBD/hto/ibd2_hash_labelled.csv", row.names = 1)
tcr.2 <- read.csv("/Users/akornberg/Documents/Columbia/Han Lab/Data Analysis/2019_11_IBD/tcr/filtered_contig_annotations.csv")
cd4.cd8.2 <- read.csv("/Users/akornberg/Documents/Columbia/Han Lab/Data Analysis/2019_11_IBD/adt/ibd2_CD4CD8_assigned.csv", row.names = 1)

# Add in Metadata
ibd2 <- AddMetaData(object=ibd2, metadata = hto.2)
ibd2 <- AddMetaData(object=ibd2, metadata = cd4.cd8.2)

ibd2[["ADT"]]<-CreateAssayObject(counts= adt.2)
ibd2<- NormalizeData(ibd2, assay = "ADT", normalization.method = "CLR")
ibd2 <- ScaleData(ibd2, assay = "ADT")

shortADTnames.2<-row.names(ibd2[["ADT"]])
ADTnames.2<-paste0("adt_", shortADTnames.2)

## TCR Pre-processing. For this we want to:
# Remove the "-1" from the barcode
# Separate the chains to plot later
# Separate productive CDR3s (full length, start codon, nonstop, in-frame, CDR3 exists, and structure)
# Remove duplicates (to just remove one duplicate: x <- ibd.tcr[!duplicated(ibd.tcr$barcode), ])
tcr.2$barcode <- gsub("-1", "", tcr.2$barcode)

TCRproductive.2<-subset(tcr.2, productive == "True")

## examine what chains were detected among the productive barcodes (in this case only TRA and TRB which is good given that's what we sorted on)
print(table(TCRproductive.2$chain))

## seperate into TCA and TRA 
splitTCR.2<-split(TCRproductive.2, TCRproductive.2$chain, drop=TRUE)

# Im sure there is a better way to do this with lapply but I will leave it for another day. seperate the dataframes out and rename their columns to be unique
Multi.2<-splitTCR.2[["Multi"]]
colnames(Multi.2)<-paste0("Multi.",colnames(Multi.2))

TRA.2<-splitTCR.2[["TRA"]]
colnames(TRA.2)<-paste0("TRA.",colnames(TRA.2))

TRB.2<-splitTCR.2[["TRB"]]
colnames(TRB.2)<-paste0("TRB.",colnames(TRB.2))
TRB.2$cdr3b <- TRB.2$TRB.cdr3
## Remove duplicated cell IDs
toRemove.2<-table(TRB.2$TRB.barcode)>1
removeNames.2<-names(which(toRemove.2==TRUE))
TRBfiltered.2<-subset(TRB.2, !TRB.barcode %in% removeNames.2)
View(TRBfiltered.2)
## add counts of TRB.cdr3
TRBfiltered.2 <- transform(TRBfiltered.2, TRB.CDR3count = ave(TRB.barcode, TRB.cdr3, FUN = length)) 
TRBfiltered.2$TRB.CDR3count<-as.numeric(as.character(TRBfiltered.2$TRB.CDR3count))
View(TRBfiltered.2)
## create separate factor for most highly expressed CDR3's
highClones<-names(which(table(TRBfiltered.2$TRB.cdr3)>10))
TRBfiltered.2$TRB.topClone<-ifelse(TRBfiltered.2$TRB.cdr3 %in% highClones, as.character(TRBfiltered.2$TRB.cdr3), NA)
## set rownames as cell barcodes
rownames(TRBfiltered.2)<-TRBfiltered.2$TRB.barcode
### add onto our seurat object
ibd2<-AddMetaData(object=ibd2, metadata = TRBfiltered.2)

# Frequency of Clones
ibd2$TRB.repeatedClones<-ifelse(ibd2$TRB.CDR3count > 1, as.numeric(as.character(ibd2$TRB.CDR3count)),NA)

CD3percentage<-data.frame(cell.ID=factor(),TRB.cdr3=factor(),totalDetectedTCR=numeric(),totalDetectedClone =numeric(),percentClone=numeric())

View(CD3sum) 

for(i in levels(as.factor(ibd2$patient))){
  workingDF<-subset(ibd2@meta.data, patient == i)
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

ibd2 <- AddMetaData(ibd2, metadata = CD3percentage)

# Calculate CDR3B per sample
ibd2$patient.sample <- substr(ibd2$sample, start = 1, stop = 4)

CD3percentage<-data.frame(cell.ID=factor(),TRB.cdr3=factor(),totalDetectedTCR=numeric(),totalDetectedClone =numeric(),percentClone=numeric())
for(i in levels(as.factor(ibd2$patient.sample))){
  if(i=='860A') next
  workingDF<-subset(ibd2@meta.data, patient.sample == i)
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

tcr.by.sample.2 <- select(CD3percentage, "totalDetectedTCR", "Freq", "percentClone")
colnames(tcr.by.sample.2) <- c("TCR.sample", "freq.sample", "pC.sample")
ibd2 <- AddMetaData(ibd2, metadata = tcr.by.sample.2)
#Begin Seurat Analysis
ibd2[["percentMT"]]<-PercentageFeatureSet(ibd2, pattern = "^MT-")

plot1 <- FeatureScatter(ibd2, feature1 = "nCount_RNA", feature2 = "percentMT") + geom_hline(yintercept = 25)
plot2 <- FeatureScatter(ibd2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_hline(yintercept = 3000)
CombinePlots(plots = list(plot1, plot2))

ibd2 <- subset(ibd2, subset = nFeature_RNA > 5 & nFeature_RNA < 3000 & percentMT < 25)
ibd2$experiment <- "ibd2"

ibd2 <- SCTransform(ibd2, vars.to.regress = "percentMT", verbose = FALSE)

# These are now standard steps in the Seurat workflow for visualization and clustering
ibd2 <- RunPCA(ibd2, verbose = FALSE)
ibd2 <- RunUMAP(ibd2, dims = 1:25, verbose = FALSE)
ibd2 <- FindNeighbors(ibd2, dims = 1:25, verbose = FALSE)
ibd2 <- FindClusters(ibd2, verbose = FALSE, resolution = 0.8)
DimPlot(ibd2, label = TRUE, group.by = "label") + NoLegend()


###################################
###################################
#END HERE#
###################################


