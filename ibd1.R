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
ibd.1 <- Read10X(data.dir = "/Users/akornberg/Documents/Columbia/Han Lab/Data Analysis/2019_08_IBD/GEX/filtered_feature_bc_matrix/")

#Remove TCR genes
features.1 <- ibd.1@Dimnames[[1]]

features.1.tcr <- grep(c("^TRAV|^TRAJ|^TRAC|^TRBV|^TRBD|^TRBJ|^TRBC"), features.1, value = TRUE)

ibd.1.data <- ibd.1[!rownames(ibd.1) %in% features.1.tcr, ]

ibd1 <- CreateSeuratObject(counts = ibd.1.data, project="ibd1")

#Load in Assay/Metadata
adt.1 <- read.csv("/Users/akornberg/Documents/Columbia/Han Lab/Data Analysis/2019_08_ibd/AdT/2019_08_ibd_AdT.csv", row.names = 1)
hto.1 <- read.csv("/Users/akornberg/Documents/Columbia/Han Lab/Data Analysis/2019_08_ibd/HTO/2019_08_ibd_hash_general.csv", row.names = 1)
tcr.1 <- read.csv("/Users/akornberg/Documents/Columbia/Han Lab/Data Analysis/2019_08_ibd/TCR/filtered_contig_annotations_ibd_Aug_2019.csv")
cd4.cd8.1 <- read.csv("/Users/akornberg/Documents/Columbia/Han Lab/Data Analysis/2019_08_IBD/AdT/ibd1_CD4CD8_assigned.csv", row.names = 1)

# Add in Metadata
ibd1 <- AddMetaData(object=ibd1, metadata = hto.1)
ibd1 <- AddMetaData(object=ibd1, metadata = cd4.cd8.1)

ibd1[["ADT"]]<-CreateAssayObject(counts= adt.1)
ibd1<- NormalizeData(ibd1, assay = "ADT", normalization.method = "CLR")
ibd1 <- ScaleData(ibd1, assay = "ADT")

shortADTnames.1<-row.names(ibd1[["ADT"]])
ADTnames.1<-paste0("adt_", shortADTnames.1)

## TCR Pre-processing. For this we want to:
# Remove the "-1" from the barcode
# Separate the chains to plot later
# Separate productive CDR3s (full length, start codon, nonstop, in-frame, CDR3 exists, and structure)
# Remove duplicates (to just remove one duplicate: x <- ibd.tcr[!duplicated(ibd.tcr$barcode), ])
tcr.1$barcode <- gsub("-1", "", tcr.1$barcode)

TCRproductive.1<-subset(tcr.1, productive == "True")

## examine what chains were detected among the productive barcodes (in this case only TRA and TRB which is good given that's what we sorted on)
print(table(TCRproductive.1$chain))

## seperate into TCA and TRA 
splitTCR.1<-split(TCRproductive.1, TCRproductive.1$chain, drop=TRUE)

# Im sure there is a better way to do this with lapply but I will leave it for another day. seperate the dataframes out and rename their columns to be unique
Multi.1<-splitTCR.1[["Multi"]]
colnames(Multi.1)<-paste0("Multi.",colnames(Multi.1))

TRA.1<-splitTCR.1[["TRA"]]
colnames(TRA.1)<-paste0("TRA.",colnames(TRA.1))

TRB.1<-splitTCR.1[["TRB"]]
colnames(TRB.1)<-paste0("TRB.",colnames(TRB.1))
TRB.1$cdr3b <- TRB.1$TRB.cdr3
## Remove duplicated cell IDs
toRemove.1<-table(TRB.1$TRB.barcode)>1
removeNames.1<-names(which(toRemove.1==TRUE))
TRBfiltered.1<-subset(TRB.1, !TRB.barcode %in% removeNames.1)
## add counts of TRB.cdr3
TRBfiltered.1 <- transform(TRBfiltered.1, TRB.CDR3count = ave(TRB.barcode, TRB.cdr3, FUN = length)) 
TRBfiltered.1$TRB.CDR3count<-as.numeric(as.character(TRBfiltered.1$TRB.CDR3count))
## create separate factor for most highly expressed CDR3's
highClones<-names(which(table(TRBfiltered.1$TRB.cdr3)>10))
TRBfiltered.1$TRB.topClone<-ifelse(TRBfiltered.1$TRB.cdr3 %in% highClones, as.character(TRBfiltered.1$TRB.cdr3), NA)
## set rownames as cell barcodes
rownames(TRBfiltered.1)<-TRBfiltered.1$TRB.barcode
### add onto our seurat object
ibd1<-AddMetaData(object=ibd1, metadata = TRBfiltered.1)

# Frequency of Clones
ibd1$TRB.repeatedClones<-ifelse(ibd1$TRB.CDR3count > 1, as.numeric(as.character(ibd1$TRB.CDR3count)),NA)

CD3percentage<-data.frame(cell.ID=factor(),TRB.cdr3=factor(),totalDetectedTCR=numeric(),totalDetectedClone =numeric(),percentClone=numeric())

for(i in levels(as.factor(ibd1$patient))){
  workingDF<-subset(ibd1@meta.data, patient == i)
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

ibd1 <- AddMetaData(ibd1, metadata = CD3percentage)

# Calculate CDR3B per sample
ibd1$patient.sample <- substr(ibd1$sample, start = 1, stop = 4)

CD3percentage<-data.frame(cell.ID=factor(),TRB.cdr3=factor(),totalDetectedTCR=numeric(),totalDetectedClone =numeric(),percentClone=numeric())
for(i in levels(as.factor(ibd1$patient.sample))){
  workingDF<-subset(ibd1@meta.data, patient.sample == i)
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

tcr.by.sample.1 <- select(CD3percentage, "totalDetectedTCR", "Freq", "percentClone")
colnames(tcr.by.sample.1) <- c("TCR.sample", "freq.sample", "pC.sample")
ibd1 <- AddMetaData(ibd1, metadata = tcr.by.sample.1)
#Begin Seurat Analysis
ibd1[["percentMT"]]<-PercentageFeatureSet(ibd1, pattern = "^MT-")

plot1 <- FeatureScatter(ibd1, feature1 = "nCount_RNA", feature2 = "percentMT") + geom_hline(yintercept = 25)
plot2 <- FeatureScatter(ibd1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_hline(yintercept = 2500)
CombinePlots(plots = list(plot1, plot2))

ibd1 <- subset(ibd1, subset = nFeature_RNA > 5 & nFeature_RNA < 2500 & percentMT < 25)
ibd1$experiment <- "ibd1"

###################################
#END HERE#
###################################

ibd1 <- SCTransform(ibd.1, vars.to.regress = "percentMT", verbose = FALSE)

## EDA plot
VlnPlot(ibd.1, features = c("nFeature_RNA", "nCount_RNA", "percentMT"), ncol = 3)
##
plot1 <- FeatureScatter(ibd.1, feature1 = "nCount_RNA", feature2 = "percentMT")
plot2 <- FeatureScatter(ibd.1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))

## MV plot
top10<-head(VariableFeatures(ibd.1),10)
## plot variable features w/top 10 most variable genes labeled
plot1 <- VariableFeaturePlot(ibd.1)
plot2 <- LabelPoints(plot = plot1, points = top10)
CombinePlots(plots = list(plot1, plot2))

#Run PCA + Choose PC's to include in the analysis

ibd.1 <- RunPCA(ibd.1, verbose = FALSE)

ElbowPlot(ibd.1, ndims= 50)
ibd.1 <- RunUMAP(ibd.1, dims = 1:30, verbose = FALSE)
ibd.1 <- FindNeighbors(ibd.1, dims = 1:30, verbose = FALSE)
ibd.1 <- FindClusters(ibd.1, verbose = FALSE)
DimPlot(ibd.1, label = TRUE) + NoLegend()

#Identify markers of each cluster
clusterMarkers <- FindAllMarkers(ibd.1, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
topMarkers<-clusterMarkers %>% group_by(cluster) %>% top_n(n = 25, wt = avg_logFC)

shortMarkers<-clusterMarkers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
shortM<-shortMarkers[,c(6:7)]
shortM$rank<-c(1:10)
shortM %>% spread(rank, gene)

#Add in (previously loaded) AdT data as an assay object
ibd.1[["ADT"]]<-CreateAssayObject(counts= adt)

ibd.1<- NormalizeData(ibd.1, assay = "ADT", normalization.method = "CLR")
ibd.1 <- ScaleData(ibd.1, assay = "ADT")

shortADTnames<-row.names(ibd.1[["ADT"]])
ADTnames<-paste0("adt_", shortADTnames)

# FEATURE PLOT - Can show AdT or GEX data. Cutoffs tend to make AdT data appear more clean

FeaturePlot(ibd.1, features = c("CD2", "B2M", "PTPRC", "ATP1B3"), ncol=4, max.cutoff = 10 )
DefaultAssay(ibd.1) <- "RNA"
View(ibd.1@meta.data$stim.status)

View(adt_test2)
cluster0.markers <- FindMarkers(ibd.1, ident.1 = , min.pct = 0.25)
cluster8.markers <- FindMarkers(ibd.1, ident.1 = 0)
class(cluster8.markers)

#Volcano Plot
Idents(object = ibd.1) <- ibd.1@meta.data$seurat_clusters
marker0.10 <- FindMarkers(ibd.1, ident.1 = 0 , ident.2 = 10, logfc.threshold = 0, min.pct = 0.01) 

EnhancedVolcano(marker0.10, lab=rownames(marker0.10), x = "avg_logFC", y = "p_val")

write.csv(marker0.10, file="DEgenes_0.10.csv")

ibd.1.markers <- FindAllMarkers(ibd.1, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
View(ibd.1.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC))

cluster8.markers <- FindMarkers(ibd.1, ident.1 = 8, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
cluster8.markers

ggplot(rna_test2, aes(x=rna_test2$SELL)) +  geom_histogram(binwidth =.01)
ggplot(adt_test2, aes(x=adt_test2$CD25)) +  geom_histogram(binwidth =.01)

scale_x_continuous(limits=c(0,400))

RidgePlot(ibd.1, features = c("adt_CD4", "adt_CD8"), ncol = 2)

# Add (previously loaded) HTO data as metadata
row.names(hto)<-hto$cell.ID
ibd.1<-AddMetaData(object=ibd.1, metadata = hto)

# Create new column with dx + activity
ibd.1@meta.data$dx.activity<-paste(ibd.1@meta.data$dx,ibd.1@meta.data$activity)

## counts of how many cells from each sample we got
test <- table(hto$sample)
View(test)
ggplot(data=test, aes(x=Var1, y=Freq)) + geom_bar(stat = "identity")
test2 <- substr(test$Var1, 1, 4) 
View(test2)
test3 <- gsub("nega", "negative", test2)
View(test3)
samplefreq <- cbind(test, test3)
colnames(samplefreq) <- c("Specific Sample", "Assigned Cells", "Sample" )
View(samplefreq)

write.csv(test, "a_harder_version_frequency_of_cells_practice.csv")
## total cells - should be same as GEX barcodes (?)
sum(table(hto$sample))

# Assign Colors for Plotting
samplecols <- c("brown1","darkred",
                "darkorange","darkorange3","darkgoldenrod1","darkgoldenrod3",
                "darkolivegreen1","darkolivegreen4","chartreuse","chartreuse3",
                "blue","blue4","cadetblue1","cadetblue4","cyan", "cyan4",
                "darkorchid1","darkorchid4", "gray90", "gray90")

activitycolor <- c("#C5050C", "#002D72","#C99E10","#DADFE1")
activitylabel <- c("Active", "Control", "Inactive", "Unknown")

stimcolor <- c("#DADFE1", "#DADFE1","#002D72","#FF5910")
stimlabel <- c("Multiplet", "Negative", "Nonstim", "Stim")

regioncolor <- c("#003D73", "#0878A4","#C05640","#D59B2D","#93A806", "#DADFE1")
regionlabel <- c("Ascending Colon","Ileum", "Rectum", "Sigmoid Colon", "Splenic Flexure", "Unknown")

dxcolor <- c("#C5050C", "#002D72", "#C99E10", "#DADFE1")
dxlabel <- c("Chron's Disease", "Healthy", "Ulcerative Colitis","Unknown")
# Group by sample
DimPlot(ibd.1, group.by = "sample", pt=1, cols = samplecols, na.value = "grey90")

# Group by condition
DimPlot(ibd.1, group.by = "activity", pt=1) + 
  ggtitle("Inflammation") + theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values=activitycolor, labels=activitylabel)


DimPlot(ibd.1, group.by = "stim.status", pt=1) +
  ggtitle("Stim Status") + theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values=stimcolor, labels=stimlabel)

# Group by patient
DimPlot(ibd.1, group.by = "stim.status", pt=1, split.by = "patient") +
  ggtitle("Stim Status by Patient") + theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values=stimcolor, labels=stimlabel)


DimPlot(ibd.1, group.by = "activity", pt=1, split.by = "patient") + 
  ggtitle("Inflammation Activity by Patient") + theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values=activitycolor, labels=activitylabel)

# Group by region and disease
DimPlot(ibd.1, group.by = "region", pt=1) +
  ggtitle("Region Plot") + theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values=regioncolor, labels=regionlabel)


DimPlot(ibd.1, group.by = "dx", pt=1) +
  ggtitle("Diagnosis ") + theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values=dxcolor, labels=dxlabel)

# Activity by region and disease
DimPlot(ibd.1, group.by = "activity", pt=1, split.by = "region") +
  ggtitle("Activity by Region") + theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values=activitycolor, labels=activitylabel)

DimPlot(ibd.1, group.by = "activity", pt=1, split.by = "dx") +
  ggtitle("Activity by Dx") + theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values=activitycolor, labels=activitylabel)


DimPlot(ibd.1, group.by = "stim.status", pt=1, split.by = "dx") +
  ggtitle("Stim status by Dx") + theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values=stimcolor, labels=stimlabel)

# Split by underlying disease and activity
DimPlot(ibd.1, group.by = "stim.status", pt=1, split.by = "dx.activity") +
  ggtitle("Stim status by Dx Activity") + theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values=stimcolor, labels=stimlabel)

# Activity by stim
DimPlot(ibd.1, group.by = "stim.status", pt=1, split.by = "activity") +
  ggtitle("Stim by Activity") + theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values=stimcolor, labels=stimlabel)

## TCR Pre-processing. For this we want to:
# Remove the "-1" from the barcode
# Separate the chains to plot later
# Separate productive CDR3s (full length, start codon, nonstop, in-frame, CDR3 exists, and structure)
# Remove duplicates (to just remove one duplicate: x <- ibd.1.tcr[!duplicated(ibd.1.tcr$barcode), ])
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
TRB$cdr3b <- TRB$TRB.cdr3
## Remove duplicated cell IDs
toRemove<-table(TRB$TRB.barcode)>1
removeNames<-names(which(toRemove==TRUE))
TRBfiltered<-subset(TRB, !TRB.barcode %in% removeNames)
View(TRBfiltered)
## add counts of TRB.cdr3
TRBfiltered <- transform(TRBfiltered, TRB.CDR3count = ave(TRB.barcode, TRB.cdr3, FUN = length)) 
TRBfiltered$TRB.CDR3count<-as.numeric(as.character(TRBfiltered$TRB.CDR3count))
View(TRBfiltered)
## create separate factor for most highly expressed CDR3's
highClones<-names(which(table(TRBfiltered$TRB.cdr3)>10))
TRBfiltered$TRB.topClone<-ifelse(TRBfiltered$TRB.cdr3 %in% highClones, as.character(TRBfiltered$TRB.cdr3), NA)
## set rownames as cell barcodes
rownames(TRBfiltered)<-TRBfiltered$TRB.barcode
### add onto our seurat object
ibd.1<-AddMetaData(object=ibd.1, metadata = TRBfiltered)
### Feature plot with clonotype frequencies 
FeaturePlot(ibd.1, "TRB.CDR3count", cols = c("lightgoldenrod2", "firebrick")) + 
  ggtitle("TCRB Frequency")

FeaturePlot(ibd.1, "TRB.CDR3count", split.by = "patient", cols= c("lightgoldenrod2", "firebrick")) +
  ggtitle("TCRB Frequency by Patient") 

DimPlot(ibd.1, group.by = "TRB.topClone") + 
  ggtitle("CDR3B Appearing > 10") + theme(plot.title = element_text(hjust = 0.5))

# Add metadata column combining stim and pt factors
ibd.1$cloneBySample<-ifelse(ibd.1$TRB.topClone != "NA", paste0(ibd.1$TRB.topClone, 
                                                           "_",as.character(ibd.1$patient),as.character(ibd.1$stim.status)), NA)

DimPlot(ibd.1, group.by = "cloneBySample") + 
  ggtitle("CDR3B Expanded > 10x by Sample") + theme(plot.title = element_text(hjust = 0.5))

# Add metadata column combining activity and pt factors
ibd.1$cloneBySampleActivity<-ifelse(ibd.1$TRB.topClone != "NA", paste0(ibd.1$TRB.topClone, "_",as.character(ibd.1$patient),as.character(ibd.1$activity)), NA)

DimPlot(ibd.1, group.by = "cloneBySampleActivity", pt.size = 1) + 
  ggtitle("CDR3B Expanded > 10x by Sample") + theme(plot.title = element_text(hjust = 0.5))

table(ibd.1$cloneBySampleActivity)

# Frequency of Clones
ibd.1$TRB.repeatedClones<-ifelse(ibd.1$TRB.CDR3count > 1, as.numeric(as.character(ibd.1$TRB.CDR3count)),NA)

CD3percentage<-data.frame(cell.ID=factor(),TRB.cdr3=factor(),totalDetectedTCR=numeric(),totalDetectedClone =numeric(),percentClone=numeric())

View(CD3sum) 

for(i in levels(as.factor(ibd.1$patient))){
  workingDF<-subset(ibd.1@meta.data, patient == i)
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

ibd.1 <- AddMetaData(ibd.1, metadata = CD3percentage)

##
FeaturePlot(ibd.1, features = "TRB.repeatedClones", split.by = "patient", cols =c("lightgoldenrod2", "firebrick")) + 
  ggtitle("CDR3B Frequency by Patient") +
  theme(plot.title = element_text(hjust = 0.5))

FeaturePlot(ibd.1, features = "percentClone") +
  scale_color_gradient(low="lightgoldenrod2", high="firebrick", na.value = "grey90") +
  ggtitle("CDR3B Frequency") + theme(plot.title = element_text(hjust = 0.5))

TRBfiltered$new_id <- ifelse(grepl("CAISEESQETQYF", TRBfiltered$TRB.topClone), "clone_100", "other")
ibd.1$
  View(TRBfiltered)

Idents(object = ibd.1) <- ibd.1@meta.data$seurat_clusters
high.clone.marker <- FindMarkers(ibd.1, ident.1 = 'clone_100', only.pos = TRUE)
View(high.clone.marker)
write.csv(high.clone.marker, "clone_100.csv")


View(ibd.1$orig.ident)
View(ibd.1$TRB.CDR3count)
Idents(ibd.1)

View(hto$cell.ID)
hto$cellTag <- hto$cell.ID

View(cellTag)







# Proportion of CDR3s detected in individual sample
CD3percentage2<-data.frame(cell.ID=factor(),TRB.cdr3=factor(),totalDetectedTCR=numeric(),totalDetectedClone =numeric(),percentClone=numeric())
View(ibd.1$patient.activity)

View(levels(as.factor(ibd.1$patient)))
View(CD3sum)

for(i in levels(as.factor(ibd.1$patient.activity))){
  workingDF<-subset(ibd.1@meta.data, patient.activity == i)
  resultsDF<-select(workingDF,cell.ID, TRB.cdr3)
  resultsDF$TRB.cdr3<-as.factor(as.character(resultsDF$TRB.cdr3))}
CD3countTable<- as.data.frame(table(as.factor(as.character(workingDF$TRB.cdr3))))
CD3sum<-sum(CD3countTable$Freq)
resultsDF$totalDetectedTCR<-CD3sum
resultsDF$TRB.cdr3<-as.factor(as.character(resultsDF$TRB.cdr3))
resultsDF<-merge(resultsDF,CD3countTable, by.x = "TRB.cdr3", by.y = "Var1")
resultsDF<-resultsDF%>%rename(totalDetectedClone = Freq)
resultsDF$percentCloneperactivity<-resultsDF$totalDetectedClone/resultsDF$totalDetectedTCR
CD3percentage2<-rbind(CD3percentage2,resultsDF)
}

row.names(CD3percentage2)<-CD3percentage2$cell.ID
View(CD3percentage)

ibd.1 <- AddMetaData(ibd.1, metadata = CD3percentage2)
View(ibd.1@meta.data)

cdr3b.freq <- ibd.1@meta.data
test <- subset(cdr3b.freq, cdr3b.freq$patient == 821)
test2 <- test[!is.na(test$percentCloneperactivity),]
active <- subset(test2, test2$activity == "active")
inactive <- subset(test2, test2$activity == "inactive")

test5 <- as.data.frame(table(inactive$cdr3b))
View(test3)
test5$merge <- test5$Var1
inactive$merge <- inactive$cdr3b
test4 <- merge(active, test3, by="merge")
test6 <- merge(inactive, test5, by = "merge")
test6$inactive <- test6$percentCloneperactivity
test7 <- merge(test4, test6, by = "merge")
View(test7)
ggplot(test, aes(x=test$activity=="active", y=test$activity=="inactive")) + geom_point()
#Cluster on AdT's
DefaultAssay(ibd.1) <- "RNA"
ibd.1 <- RunPCA(ibd.1, features = rownames(ibd.1), reduction.name = "pca_adt", reduction.key = "pca_adt_", 
              verbose = FALSE)
DimPlot(ibd.1, reduction = "pca_adt")

adt.data <- GetAssayData(ibd.1, slot = "data")
adt.dist <- dist(t(adt.data))

ibd.1[["rnaClusterID"]] <- Idents(ibd.1)

ibd.1[["tsne_adt"]] <- RunTSNE(adt.dist, assay = "ADT", reduction.key = "adtTSNE_")
ibd.1[["adt_snn"]] <- FindNeighbors(adt.dist)$snn
ibd.1 <- FindClusters(ibd.1, resolution = 0.2, graph.name = "adt_snn")
tsne_adtClusters <- DimPlot(ibd.1, reduction = "tsne_adt", pt.size = 0.5) + NoLegend()
tsne_adtClusters    
