#Load in required packages
library(ggplot2)
library(Seurat)
library(dplyr)
library(ggridges)
library(sctransform)
library(reshape2)
library(tidyr)
library(Matrix)
library(data.table)
library(EnhancedVolcano)
library(rmarkdown)
library(MAST)

#set Working Directory
setwd("~/OneDrive - cumc.columbia.edu/Samhita Rao/2020")

#exp4, you have bead a and bead b
#gex
#add in gex raw data. this is a matrix
sr.beada <- Read10X(data.dir = "~/OneDrive - cumc.columbia.edu/Samhita Rao/2019/May2019toDec2019/10X/EXP 4/Data/GEX/AS001/filtered_feature_bc_matrix/")
sr.beadb <- Read10X(data.dir = "~/OneDrive - cumc.columbia.edu/Samhita Rao/2019/May2019toDec2019/10X/EXP 4/Data/GEX/AS002/filtered_feature_bc_matrix/")

#row name of all 'features'
features.a <- sr.beada@Dimnames[[1]]
features.b <- sr.beadb@Dimnames[[1]]

#pulls out a list of all genes that start with trav... all tcr genes (makes a vector)
features.a.tcr <- grep(c("^Trav|^Traj|^Trac|^Trbv|^Trbd|^Trbj|^Trbc"), features.a, value = TRUE)
features.b.tcr <- grep(c("^Trav|^Traj|^Trac|^Trbv|^Trbd|^Trbj|^Trbc"), features.b, value = TRUE)

#from the data, we want to remove every tcr gene form this matrix 
sr.beada.data <- sr.beada[!rownames(sr.beada) %in% features.a.tcr, ]
sr.beadb.data <- sr.beadb[!rownames(sr.beada) %in% features.b.tcr, ]

#this is where we create the object for bead a and bead b 
sr.a <- CreateSeuratObject(counts = sr.beada.data, project = "sr_2019_11_a")
sr.b <- CreateSeuratObject(counts = sr.beadb.data, project="sr_2019_11_b")


#load hto files and add as metadata
hto.a <- read.csv(file.choose(), row.names = 1)
hto.b <- read.csv(file.choose(), row.names = 1)
sr.a <- AddMetaData(object=sr.a, metadata = hto.a)
sr.b <- AddMetaData(object=sr.b, metadata = hto.b)

#load adt files and add as metadata
adt.a <- read.csv(file.choose(), row.names = 1)
adt.b <- read.csv(file.choose(), row.names = 1)

sr.a[["ADT"]]<-CreateAssayObject(counts= adt.a)
sr.b[["ADT"]]<-CreateAssayObject(counts= adt.b)

sr.a<- NormalizeData(sr.a, assay = "ADT", normalization.method = "CLR")
sr.a <- ScaleData(sr.a, assay = "ADT")

shortADTnames.a<-row.names(sr.a[["ADT"]])
ADTnames.a<-paste0("adt_", shortADTnames.a)

sr.b<- NormalizeData(sr.b, assay = "ADT", normalization.method = "CLR")
sr.b <- ScaleData(sr.b, assay = "ADT")

shortADTnames.b<-row.names(sr.b[["ADT"]])
ADTnames.b<-paste0("adt_", shortADTnames.b)

#load in cd4cd8_beada/b_assigned files
cd4.cd8.a <- read.csv(file.choose(), row.names = 1)
cd4.cd8.b <- read.csv(file.choose(), row.names = 1)

sr.a <- AddMetaData(object=sr.a, metadata = cd4.cd8.a)
sr.b <- AddMetaData(object=sr.b, metadata = cd4.cd8.b)


#merge objects
sr.combined <- merge(sr.a, y=sr.b, project = "sr_2019_11")

#TCR Preprocessing
#load in bead a tcr
tcr.a <- read.csv("~/OneDrive - cumc.columbia.edu/Samhita Rao/2019/May2019toDec2019/10X/EXP 4/Data/TCR/filtered_contig_annotations_beada_b6.csv")
#load in bead b tcr
tcr.b <- read.csv("~/OneDrive - cumc.columbia.edu/Samhita Rao/2019/May2019toDec2019/10X/EXP 4/Data/TCR/filtered_contig_annotations_beadb_b7.csv")

#stash bead a and bead b 
tcr.a$old.barcode <- tcr.a$barcode
tcr.b$old.barcode <- tcr.b$barcode

#in 10x, tcr names comes with "-1" and we want to remove it and replace it with nothing
tcr.a$barcode <- gsub("-1", "", tcr.a$barcode)
tcr.b$barcode <- gsub("-1", "", tcr.b$barcode)

#taking just the productive beta chains for bead a and bead b so here 'a'/ 'b' refer to beads and not alpha/ beta chain
TCRproductive.a<-subset(tcr.a, productive == "True")
TCRproductive.b<-subset(tcr.b, productive == "True")

#print table for bead a 
print(table(TCRproductive.a$chain))

#print table for bead b
print(table(TCRproductive.b$chain))

#Seperate into TCRa and TCRb

#this is bead a 
splitTCR.a<-split(TCRproductive.a, TCRproductive.a$chain, drop=TRUE)

#this is bead b
splitTCR.b<-split(TCRproductive.b, TCRproductive.b$chain, drop=TRUE)

#alpha chain of bead a
TRA.a<-splitTCR.a[["TRA"]]
colnames(TRA.a)<-paste0("TRA.",colnames(TRA.a))

#alpha of bead b
TRA.b<-splitTCR.b[["TRA"]]
colnames(TRA.b)<-paste0("TRA.",colnames(TRA.b))

#beta of bead a 
TRB.a<-splitTCR.a[["TRB"]]
colnames(TRB.a)<-paste0("TRB.a",colnames(TRB.a))
TRB.a$cdr3b <- TRB.a$TRB.acdr3

#beta of bead b
TRB.b<-splitTCR.b[["TRB"]]
colnames(TRB.b)<-paste0("TRB.b",colnames(TRB.b))
TRB.b$cdr3b <- TRB.b$TRB.bcdr3

#TCR pre processing
#remove duplicates for bead a 
toRemove.a<-table(TRB.a$TRB.abarcode)>1
removeNames.a<-names(which(toRemove.a==TRUE))
TRBfiltered.a<-subset(TRB.a, !TRB.abarcode %in% removeNames.a)

TRBfiltered.a <- transform(TRBfiltered.a, TRB.aCDR3count = ave(TRB.abarcode, TRB.acdr3, FUN = length)) 
TRBfiltered.a$TRB.aCDR3count<-as.numeric(as.character(TRBfiltered.a$TRB.aCDR3count))

#remove duplicates for bead b 
toRemove.b<-table(TRB.b$TRB.bbarcode)>1
removeNames.b<-names(which(toRemove.b==TRUE))
TRBfiltered.b<-subset(TRB.b, !TRB.bbarcode %in% removeNames.b)

TRBfiltered.b <- transform(TRBfiltered.b, TRB.bCDR3count = ave(TRB.bbarcode, TRB.bcdr3, FUN = length)) 
TRBfiltered.b$TRB.bCDR3count<-as.numeric(as.character(TRBfiltered.b$TRB.bCDR3count))

#the old barcode is now a column barcode 
#for bead a
TRBfiltered.a$barcode <- TRBfiltered.a$TRB.aold.barcode
#for bead b
TRBfiltered.b$barcode <- TRBfiltered.b$TRB.bold.barcode

#takes all the barcodes and adds a "_1" if from bead a and "_2" if from bead b 
barcoder <- function(df, prefix, trim="\\-1"){
  
  df$barcode <- gsub(trim, "", df$barcode)
  df$barcode <- paste0(df$barcode, prefix)
  
  df
}

TRBfiltered.a <- barcoder(TRBfiltered.a, prefix = "_1")
TRBfiltered.b <- barcoder(TRBfiltered.b, prefix = "_2")

rownames(TRBfiltered.a)<-TRBfiltered.a$barcode
rownames(TRBfiltered.b)<-TRBfiltered.b$barcode


#these are the column names, you can view TRBfiltered.a or b to confirm this
trb.colnames <- c("trb.barcode", "is.cell", "contig.id", "high.confidence", "length", "chain", "v.gene", "d.gene",
                  "j.gene", "c.gene", "full.length", "productive", "cdr3", "cdr3.nt", "reads", "umis", 
                  "raw.colonotype.id", "raw.consensus.id", "old.barcode", "cdr3b", "cdr3.count", "barcode.rowname") 

#assign column names to bead a and b 
colnames(TRBfiltered.a) <- trb.colnames
colnames(TRBfiltered.b) <- trb.colnames

#bind bead a and bead b into one 
trb.bound <- rbind(TRBfiltered.a, TRBfiltered.b)

#add this metadata into main object 
sr.combined<-AddMetaData(object=sr.combined, metadata = trb.bound)
#this just creates a column called cell.id
sr.combined$cell.id <- rownames(sr.combined@meta.data)

#gex: subset
#remove low quality cells
sr.combined[["percentMT"]]<-PercentageFeatureSet(sr.combined, pattern = "^mt-")

plot1 <- FeatureScatter(sr.combined, feature1 = "nCount_RNA", feature2 = "percentMT") + geom_hline(yintercept = 42)
plot2 <- FeatureScatter(sr.combined, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_hline(yintercept = 5500)
CombinePlots(plots = list(plot1, plot2))

sr.combined <- subset(sr.combined, subset = nFeature_RNA > 5 & nFeature_RNA < 5500 & percentMT < 42)

#split by experiment
sr.turnover <- subset(sr.combined, subset = experiment == "turnover")
sr.icb <- subset(sr.combined, subset = experiment == "icb")

#for turnover, this counts cdr3b frequency per mouse 
#create a dataframe with the following column names
CD3percentage<-data.frame(cell.ID=factor(),TRB.cdr3=factor(),totalDetectedTCR=numeric(),totalDetectedClone =numeric(),percentClone=numeric())

sr.turnover$day<-as.character(sr.turnover$day)

sr.turnover$day<- ifelse((sr.turnover$day=="na"), NA, sr.turnover$day)

for(i in levels(as.factor(sr.turnover$day))){
  workingDF<-subset(sr.turnover@meta.data, day == i)
  resultsDF<-select(workingDF,cell.id, cdr3, trb.barcode, cdr3b)
  resultsDF$cdr3<-as.factor(as.character(resultsDF$cdr3))
  resultsDF$cdr3b<-as.factor(as.character(resultsDF$cdr3b))
  CD3countTable<- as.data.frame(table(as.factor(as.character(workingDF$cdr3))))
  CD3sum<-sum(CD3countTable$Freq)
  resultsDF$totalDetectedTCR<-CD3sum
  resultsDF$cdr3<-as.factor(as.character(resultsDF$cdr3))
  resultsDF$cdr3b<-as.factor(as.character(resultsDF$cdr3b))
  resultsDF<-merge(resultsDF,CD3countTable, by.x = "cdr3", by.y = "Var1")
  resultsDF$percentClone<-resultsDF$Freq/resultsDF$totalDetectedTCR
  CD3percentage<-rbind(CD3percentage,resultsDF)}
row.names(CD3percentage)<-CD3percentage$cell.id

#add as metadata into main object
sr.turnover <- AddMetaData(sr.turnover, metadata = CD3percentage)

#SC Transform --> Normalization turnover
sr.turnover <- SCTransform(sr.turnover, vars.to.regress = c("percentMT"),  verbose = FALSE)
sr.turnover <- RunPCA(sr.turnover, verbose = FALSE)
ElbowPlot(sr.turnover, ndims= 50)
sr.turnover <- RunUMAP(sr.turnover, dims = 1:30, verbose = FALSE)
sr.turnover <- FindNeighbors(sr.turnover, dims = 1:30, verbose = FALSE)
sr.turnover <- FindClusters(sr.turnover, verbose = FALSE)

DimPlot(sr.turnover, pt=0.267, label = TRUE) + NoLegend() 
DimPlot(sr.turnover, group.by = "day", pt=0.267, cols = c("gold3","gold3",
                                                          "gold3","navy","gold3", "firebrick",
                                                          "gold3", "salmon", "gold3",
                                                          "deepskyblue1", "grey90") ,na.value = "grey90")

#cluster biomarkers 
clusterMarkers <- FindAllMarkers(sr.turnover, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
topMarkers<-clusterMarkers %>% group_by(cluster) %>% top_n(n = 25, wt = avg_logFC)

shortMarkers<-clusterMarkers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
shortM<-shortMarkers[,c(6:7)]
shortM$rank<-c(1:10)
View(shortM %>% spread(rank, gene))

cluster.7.markers <- FindMarkers(sr.turnover, ident.1 = "7", ident.2 = , min.diff.pct = 0)
View(cluster.3.markers)

library(EnhancedVolcano)

EnhancedVolcano(cluster.3.markers, lab=rownames(cluster.3.markers), x = "avg_logFC", y = "p_val_adj", 
                title = "Cluster 3 gene expression")

EnhancedVolcano(cluster.3.markers,
                lab = rownames(cluster.3.markers),
                x = 'avg_logFC',
                y = 'p_val_adj',
                xlim = c(-2, 2),
                title = 'Cluster 3 gene expression',
                FCcutoff = 1,
                pointSize = 3.0,
                labSize = 3.0)


#plot cd8s/cd4s/adts on sr.turnover
FeaturePlot(sr.icb, features = c("Cd4"), pt=0.267, min.cutoff =0, max.cutoff =1)

#plot tcr percent clone for turnover
FeaturePlot(sr.turnover, features = "percentClone") +
  scale_color_gradient(low="lightgoldenrod2", high="firebrick", na.value = "grey90") +
  ggtitle("CDR3B percentClone") + theme(plot.title = element_text(hjust = 0.75))


#divide cd4 and cd8 for turnover
#Separate cd4 and cd8
sr.turnover.cd4 <- subset(sr.turnover, subset = label == "CD4")
sr.turnover.cd8 <- subset(sr.turnover, subset = label == "CD8")

#cd8 for turnover
sr.turnover.cd8 <- SCTransform(sr.turnover.cd8, vars.to.regress = c("percentMT"),  verbose = FALSE)
sr.turnover.cd8 <- RunPCA(sr.turnover.cd8, verbose = FALSE)
ElbowPlot(sr.turnover.cd8, ndims= 50)
sr.turnover.cd8 <- RunUMAP(sr.turnover.cd8, dims = 1:30, verbose = FALSE)
sr.turnover.cd8 <- FindNeighbors(sr.turnover.cd8, dims = 1:30, verbose = FALSE)
sr.turnover.cd8 <- FindClusters(sr.turnover.cd8, verbose = FALSE)

DimPlot(sr.turnover.cd8, pt=0.5, label = TRUE) + NoLegend()
DimPlot(sr.turnover.cd8, group.by = "pre_post", pt=0.3, cols = c("grey90","firebrick",
                                                                 "grey90","grey90","gold3", "orange3",
                                                                 "gold3", "red4", "gold3",
                                                                 "deepskyblue1", "grey90") ,na.value = "grey90")

clusterMarkers.cd8 <- FindAllMarkers(sr.turnover.cd8, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
topMarkers<-clusterMarkers.cd8 %>% group_by(cluster) %>% top_n(n = 25, wt = avg_logFC)

shortMarkers<-clusterMarkers.cd8 %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
shortM<-shortMarkers[,c(6:7)]
shortM$rank<-c(1:10)
View(shortM %>% spread(rank, gene))

cluster.7.markers <- FindMarkers(sr.turnover.cd8, ident.1 = "7", ident.2 = , min.diff.pct = 0.05,logfc.threshold = 0.05)
View(cluster.1.markers)

EnhancedVolcano(cluster.7.markers,
                lab = rownames(cluster.7.markers),
                x = 'avg_logFC',
                y = 'p_val_adj',
                xlim = c(-2, 2),
                title = 'CD8+ T cell cluster 7 gene expression',
                FCcutoff = 0.5,
                pointSize = 3.0,
                labSize = 3.0)

#exhaustion dot plot
ex_genes <- c("Ctla4", "Pdcd1",
              "Lag3", "Havcr2", "Tigit", "Gzmb", "Prf1")
DotPlot(object = sr.turnover.cd8, features = ex_genes)
qu_genes <- c("Cd69","Klf2", "Id3",
              "Tcf7", "Itgae")
DotPlot(object = sr.turnover.cd8, features = qu_genes)


#tcr plots for sr.turnover.cd8
FeaturePlot(sr.turnover.cd8, features = "percentClone") +
  scale_color_gradient(low="lightgoldenrod2", high="firebrick", na.value = "grey90") +
  ggtitle("CDR3B percentClone") + theme(plot.title = element_text(hjust = 0.5))


# diversity plot to compare tcrs between mice for turnover
# Cdr3b split by sample --> For diversity, expansion, and contraction analysis FOR TURNOVER:

tcr.general <- select(sr.turnover.cd8@meta.data, day, cdr3, percentClone)
tcr.general$barcode <- row.names(tcr.general)

cdr3 <- levels(as.factor(tcr.general$cdr3))

for(i in unique(as.factor(sr.turnover.cd8$day))){
  tcr.general[[i]]<-tcr.general[i,2] 
  tcr.general[ncol(tcr.general)] <- ifelse((tcr.general$day == i), tcr.general$percentClone, NA)
}

tcr.sample <- select(tcr.general, cdr3, 'D0_1', 'D0_2', 'D3_1', 'D3_2', 'D6_1', 
                     'D6_2',"D9_1", 'D9_2', 'D12_1', 'D12_2')

cdr3.sample <- as.data.frame(setDT(tcr.sample)[, lapply(.SD, function(x) sort(x)[1L]), by = .(cdr3)])

to.merge <- select(tcr.general, cdr3, barcode)

tcr.merged <- merge(cdr3.sample, to.merge, by = 'cdr3', all = TRUE)

rownames(tcr.merged) <- tcr.merged$barcode
tcr.merged$barcode <- NULL
colnames(tcr.merged) <- c("cdr3", 'D0_1', 'D0_2', 'D3_1', 'D3_2', 'D6_1', 
                          'D6_2',"D9_1", 'D9_2', 'D12_1', 'D12_2')

tcr.merged$d3<- ifelse((tcr.merged$pc.d3 == tcr.merged$pc.d0.d3), "Other",
                       ifelse(tcr.merged$pc.d3>tcr.merged$pc.d0.d3, "Expand", "Contract"))

tcr.merged[is.na(tcr.merged)]


sr.turnover.cd8 <- AddMetaData(sr.turnover.cd8, metadata = tcr.merged)

#To make 0
sr.turnover.cd8@meta.data[49:58][is.na(sr.turnover.cd8@meta.data[49:58])] <- 0

#To make NA
sr.turnover.cd8@meta.data[49:58][sr.turnover.cd8@meta.data[49:58] == 0] <- NA

ggplot(sr.turnover.cd8@meta.data, aes(x=D0_1, y=D0_2)) + geom_point(position="jitter") + 
  ggtitle("D0_1 v D0_2, scaled, with 0s") + theme(plot.title = element_text(hjust = 0.5)) + xlim(0, 0.05) + ylim(0, 0.05)


ggplot(sr.turnover@meta.data, aes(x=D0_1, y=D0_2)) + geom_point() + 
  ggtitle("D0_1 v D0_2, scaled, with 0s") + theme(plot.title = element_text(hjust = 0.5)) + scale_x_continuous(name = "D0_1", limits = c(0,0.1))
+ scale_y_continuous(name= "D0_2", limits=c(0,0.1))  

#heatmaps for days and clusters for cd8 only
tcr.cluster <- select(sr.turnover.cd8@meta.data, day, cdr3, percentClone)
tcr.cluster$barcode <- row.names(tcr.cluster)

cdr3 <- levels(as.factor(tcr.cluster$cdr3))

for(i in levels(as.factor(sr.turnover.cd8$day))){
  tcr.cluster[[i]]<-tcr.cluster[i,2]
  tcr.cluster[,ncol(tcr.cluster)] <- ifelse((tcr.cluster$day == i), tcr.cluster$percentClone, NA)
}

colnames(tcr.cluster) <- c("seurat clusters", "cdr3", "percentClone", "barcode","cluster.0","cluster.1", "cluster.2", "cluster.3", "cluster.4", "cluster.5", "cluster.6", "cluster.7", "cluster.8", "cluster.9", "cluster.10", "cluster.11", "cluster.12", "cluster.13", "cluster.14")
cluster.sample <- select(tcr.cluster, cdr3,D0_1, D0_2, D3_1,D3_2,D6_1,D6_2,D9_1,D9_2,D12_1,D12_2)

cdr3.cluster <- as.data.frame(setDT(cluster.sample)[, lapply(.SD, function(x) sort(x)[1L]), by = .(cdr3)])

to.merge <- select(tcr.cluster, cdr3, barcode)

tcr.merged.clusters <- merge(cdr3.cluster, to.merge, by = 'cdr3', all = TRUE)

rownames(tcr.merged.clusters) <- tcr.merged.clusters$barcode
tcr.merged.clusters$barcode <- NULL

tcr.merged.clusters[,2:ncol(tcr.merged.clusters)][is.na(tcr.merged.clusters[,2:ncol(tcr.merged.clusters)])] <- 0

sr.turnover.cd8 <- AddMetaData(sr.turnover.cd8, metadata = tcr.merged.clusters)

sr.turnover.cd8$shared <- ifelse((sr.turnover.cd8$D0_1 > 0 & sr.turnover.cd8$D0_2 > 0), "Shared",
                                 ifelse(sr.turnover.cd8$D3_1 > 0 & sr.turnover.cd8$D3_2 > 0, "Shared",
                                        ifelse(sr.turnover.cd8$D6_1 > 0 & sr.turnover.cd8$D6_2 > 0, "Shared",
                                               ifelse(sr.turnover.cd8$D9_1 > 0 & sr.turnover.cd8$D9_2 > 0, "Shared",
                                                      ifelse(sr.turnover.cd8$D12_1 > 0 & sr.turnover.cd8$D12_2 > 0, "Shared", NA)))))
sr.turnover.cd8$turnover.day <- factor(sr.turnover.cd8$day, levels = unique(sr.turnover.cd8$day))
sr.turnover.cd8$turnover.day <- factor(sr.turnover.cd8$turnover.day, levels =c("D0_1", "D0_2", "D3_1", "D3_2", "D6_1","D6_2", "D9_1","D9_2","D12_1","D12_2"))
#Take Total barcodes with CDR3B's - Then only unique
#test is all
#df is shared
test <- select(sr.turnover.cd8@meta.data, cdr3, percentClone, turnover.day, day, seurat_clusters)
test2 <- na.omit(test)
test3 <- unique(test2)
test3$pC <- test3$percentClone*100

test3$Y1 <- cut(test3$pC, breaks=c(0, 0.1488, 0.5633, 1.1904, 5, 20))
test3$Y2 <- cut(test3$pC, breaks=c(0,0.01, 0.25, 0.5,1, 2, 20))
test3$Y3 <- cut(test3$pC, breaks=6)

ggplot(test3, aes(x = turnover.day, y = cdr3)) +
  geom_tile(aes(fill = Y2), colour= "white") +
  scale_fill_brewer(palette = "OrRd") +
  theme(panel.grid.minor=element_blank(), axis.text.y =element_text(size = 1))

df <- select(sr.turnover.cd8@meta.data, cdr3, percentClone,day, shared, turnover.day,seurat_clusters)
df <- na.omit(df)
df <- unique(df)
df$pC <- df$percentClone*100

df$Y1 <- cut(df$pC, breaks=c(0, 0.1488, 0.5633, 1.1904, 5, 20))
df$Y2 <- cut(df$pC, breaks=c(0, 0.5,1, 2, 20))
df$Y3 <- cut(df$pC, breaks=6)

ggplot(df, aes(x = turnover.day, y = cdr3)) +
  geom_tile(aes(fill = Y2), colour= "white") +
  scale_fill_brewer(palette = "OrRd") +
  theme(panel.grid.minor=element_blank(), axis.text.y =element_text(size = 4))

#Specific Day 
df$d0 <- ifelse((df$day == "D0_1"), "D0_1", ifelse(df$day == "D0_2", "D0_2", NA))
df$d3 <- ifelse((df$day == "D3_1"), "D3_1", ifelse(df$day == "D3_2", "D3_2", NA))
df$d6 <- ifelse((df$day == "D6_1"), "D6_1", ifelse(df$day == "D6_2", "D6_2", NA))
df$d9 <- ifelse((df$day == "D9_1"), "D9_1", ifelse(df$day == "D9_2", "D9_2", NA))
df$d12 <- ifelse((df$day == "D12_1"), "D12_1", ifelse(df$day == "D12_2", "D12_2", NA))

ggplot(df, aes(x = d3, y = cdr3)) +
  geom_tile(aes(fill = Y2), colour= "white") +
  scale_fill_brewer(palette = "OrRd") +
  theme(panel.grid.minor=element_blank(), axis.text.y =element_text(size = 4))

#Clusters
ggplot(test3, aes(x = turnover.day, y = cdr3)) +
  geom_tile(aes(fill = seurat_clusters))+
  theme(panel.grid.minor=element_blank(), axis.text.y =element_text(size = 4))

sr.turnover$day.adj <- ifelse((sr.turnover$day == "D0_1"), "D0",
                       ifelse(sr.turnover$day == "D0_2", "D0", 
                       ifelse(sr.turnover$day == "D3_1", "D0",
                       ifelse(sr.turnover$day == "D6_1", "D0",
                       ifelse(sr.turnover$day == "D9_1", "D0",
                       ifelse(sr.turnover$day == "D12_1", "D0", sr.turnover$day))))))

sr.turnover$day.general <- gsub("_2", "", sr.turnover$day.adj)

#cd4 for turnover
sr.turnover.cd4 <- subset(x=sr.turnover, subset = label == "CD4")

sr.turnover.cd4 <- RunPCA(sr.turnover.cd4, verbose = FALSE)
ElbowPlot(sr.turnover.cd4, ndims= 50)
sr.turnover.cd4 <- RunUMAP(sr.turnover.cd4, dims = 1:30, verbose = FALSE)
sr.turnover.cd4 <- FindNeighbors(sr.turnover.cd4, dims = 1:30, verbose = FALSE)
sr.turnover.cd4 <- FindClusters(sr.turnover.cd4, verbose = FALSE)

DimPlot(sr.turnover.cd4, pt=0.267, label = TRUE) + NoLegend()
DimPlot(sr.turnover.cd4, group.by = "day.general", pt=0.267) + 
  scale_color_manual(values = c("gray90", "deepskyblue2", "lightgoldenrod2", "firebrick", "navy"))


#tcr plots for sr.turnover.cd4
FeaturePlot(sr.turnover.cd4, features = "percentClone") +
  scale_color_gradient(low="lightgoldenrod2", high="firebrick", na.value = "grey90") +
  ggtitle("CDR3B percentClone") + theme(plot.title = element_text(hjust = 0.5))

#cluster markers for turnover cd4
sr.turnover.cd4 <- FindAllMarkers(sr.turnover.cd4, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
topMarkers<-sr.turnover.cd4 %>% group_by(cluster) %>% top_n(n = 25, wt = avg_logFC)

shortMarkers<-sr.turnover.cd4 %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
shortM<-shortMarkers[,c(6:7)]
shortM$rank<-c(1:10)
View(shortM %>% spread(rank, gene))

turnover.cd4.c0 <- FindMarkers(sr.turnover.cd4, ident.1 = 0, ident.2 = , min.diff.pct = ,logfc.threshold = , test.use = "MAST")
View(cluster.1.markers)

FeaturePlot(sr.turnover.cd4 , features = c("Foxp3"), pt=0.267, min.cutoff =0, max.cutoff =)

#Shared
DimPlot(sr.turnover.cd4, group.by = "shared.day", pt=0.267, order = c("Shared D12_1", "Shared D12_2")) + 
  scale_color_manual(values = c("gray90", "gray90", "gray90", "gray90", "gray90", "gray90", "gray90", "gray90", "gray90", "purple4", "darkorange"))+
  ggtitle("Shared D12")















#exp4 icb 
# FOR ICB
CD3percentage2<-data.frame(cell.ID=factor(),TRB.cdr3=factor(),totalDetectedTCR=numeric(),totalDetectedClone =numeric(),percentClone=numeric())

for(i in levels(as.factor(sr.icb$mouse))){
  workingDF<-subset(sr.icb@meta.data, mouse == i)
  resultsDF<-select(workingDF,cell.id, cdr3, trb.barcode, cdr3b)
  resultsDF$cdr3<-as.factor(as.character(resultsDF$cdr3))
  resultsDF$cdr3b<-as.factor(as.character(resultsDF$cdr3b))
  CD3countTable<- as.data.frame(table(as.factor(as.character(workingDF$cdr3))))
  CD3sum<-sum(CD3countTable$Freq)
  resultsDF$totalDetectedTCR<-CD3sum
  resultsDF$cdr3<-as.factor(as.character(resultsDF$cdr3))
  resultsDF$cdr3b<-as.factor(as.character(resultsDF$cdr3b))
  resultsDF<-merge(resultsDF,CD3countTable, by.x = "cdr3", by.y = "Var1")
  resultsDF$percentClone<-resultsDF$Freq/resultsDF$totalDetectedTCR
  CD3percentage2<-rbind(CD3percentage2,resultsDF)}
row.names(CD3percentage2)<-CD3percentage2$cell.id

sr.icb <- AddMetaData(sr.icb, metadata = CD3percentage2)





#exp5
#gex: load the data which should be in the filtered_feature_bc_matrix format
exp5 <- Read10X(data.dir = "~/OneDrive - cumc.columbia.edu/Samhita Rao/2019/May2019toDec2019/10X/EXP 5/GEX/filtered_feature_bc_matrix/")
features <- sr.beada@Dimnames[[1]]
features.tcr <- grep(c("^Trav|^Traj|^Trac|^Trbv|^Trbd|^Trbj|^Trbc"), features, value = TRUE)
exp5.data <- exp5[!rownames(exp5) %in% features.tcr, ]

#gex: create seurat object 
exp5 <- CreateSeuratObject(counts = exp5.data, project = "exp5")
#hto: load data
hto <- read.csv(file.choose(), row.names = 1)
#hto: create seurat object 
exp5 <- AddMetaData(object=exp5, metadata = hto)
#adt: load data
adt <- read.csv(file.choose(), row.names = 1)
#adt: create seurat object, normalize and name each ADT
exp5[["ADT"]]<-CreateAssayObject(counts= adt)
exp5<- NormalizeData(exp5, assay = "ADT", normalization.method = "CLR")
exp5 <- ScaleData(exp5, assay = "ADT")
shortADTnames<-row.names(exp5[["ADT"]])
ADTnames<-paste0("adt_", shortADTnames)
#adt: last step is to add in the manually assigned file for cd4cd8
exp5.cd4.cd8<- read.csv(file.choose(), row.names = 1)
exp5 <- AddMetaData(object=exp5, metadata = exp5.cd4.cd8)

#TCR Preprocessing
#load in tcr file
tcr <- read.csv(file.choose())

#stash barcode
tcr$old.barcode <- tcr$barcode

#in 10x, tcr names comes with "-1" and we want to remove it and replace it with nothing
tcr$barcode <- gsub("-1", "", tcr$barcode)

#taking just the productive beta chains for bead a and bead b so here 'a'/ 'b' refer to beads and not alpha/ beta chain
TCRproductive<-subset(tcr, productive == "TRUE")

#print table 
print(table(TCRproductive$chain))

#Seperate into TCRalpha and TCRbeta
splitTCR<-split(TCRproductive, TCRproductive$chain, drop=TRUE)

#alpha chain 
TRA<-splitTCR[["TRA"]]
colnames(TRA)<-paste0("TRA",colnames(TRA))

#beta chain
TRB<-splitTCR[["TRB"]]
colnames(TRB)<-paste0("TRB",colnames(TRB))
TRB$cdr3b <- TRB$TRBcdr3

#TCR pre processing
#remove duplicates 
toRemove<-table(TRB$TRBbarcode)>1
removeNames<-names(which(toRemove==TRUE))
TRBfiltered<-subset(TRB, !TRBbarcode %in% removeNames)

TRBfiltered <- transform(TRBfiltered, TRBCDR3count = ave(TRBbarcode, TRBcdr3, FUN = length)) 
TRBfiltered$TRBCDR3count<-as.numeric(as.character(TRBfiltered$TRBCDR3count))

#the old barcode is now a column barcode 

TRBfiltered$barcode <- TRBfiltered$TRBold.barcode

#had to add in an extra column for cdr3b
TRBfiltered$cdr3b<- TRBfiltered$TRBcdr3

#these are the column names, you can view TRBfiltered.a or b to confirm this
trb.colnames <- c("trb.barcode", "is.cell", "contig.id", "high.confidence", "length", "chain", "v.gene", "d.gene",
                  "j.gene", "c.gene", "full.length", "productive", "cdr3", "cdr3.nt", "reads", "umis", 
                  "raw.colonotype.id", "raw.consensus.id", "old.barcode", "cdr3b", "cdr3.count", "barcode.rowname") 


#assign column names 
colnames(TRBfiltered) <- trb.colnames
row.names(TRBfiltered) <- TRBfiltered$trb.barcode

#add this metadata into main object 
exp5<-AddMetaData(object=exp5, metadata = TRBfiltered)

#this just creates a column called cell.id
exp5$cell.id <- rownames(exp5@meta.data)




CD3percentage2<-data.frame(cell.ID=factor(),TRB.cdr3=factor(),totalDetectedTCR=numeric(),totalDetectedClone =numeric(),percentClone=numeric())

for(i in levels(as.factor(exp5$mouse))){
  workingDF<-subset(exp5@meta.data, mouse == i)
  resultsDF<-select(workingDF,cell.id, cdr3, trb.barcode, cdr3b)
  resultsDF$cdr3<-as.factor(as.character(resultsDF$cdr3))
  CD3countTable<- as.data.frame(table(as.factor(as.character(workingDF$cdr3))))
  CD3sum<-sum(CD3countTable$Freq)
  resultsDF$totalDetectedTCR<-CD3sum
  resultsDF$cdr3<-as.factor(as.character(resultsDF$cdr3))
  resultsDF<-merge(resultsDF,CD3countTable, by.x = "cdr3", by.y = "Var1")
  resultsDF$percentClone<-resultsDF$Freq/resultsDF$totalDetectedTCR
  CD3percentage2<-rbind(CD3percentage2,resultsDF)}
row.names(CD3percentage2)<-CD3percentage2$cell.id

exp5 <- AddMetaData(exp5, metadata = CD3percentage2)

#combine exp4 icb and exp5
#Make a new column for each experiment that will allow the merge program to call the right cells
sr.icb$experiment <- "exp4"
exp5$experiment <- "exp5"
#merge 
exp4exp5 <- merge(exp5, y = sr.icb, add.cell.ids = c("exp5", "exp4"), project = "exp4exp5")
#post merge processing 
exp4exp5.list <- SplitObject(exp4exp5, split.by = "experiment")
exp4exp5.list <- exp4exp5.list[c("exp5", "exp4")]
for (i in 1:length(exp4exp5.list)) {
  exp4exp5.list[[i]] <- SCTransform(exp4exp5.list[[i]], verbose = FALSE)
}
#prepSCTintegration, ensuring all pearson residuals have been calculated
exp4exp5.features <- SelectIntegrationFeatures(object.list = exp4exp5.list, nfeatures = 3000)
options(future.globals.maxSize = 4000 * 1024^2)
exp4exp5.list <- PrepSCTIntegration(object.list = exp4exp5.list, anchor.features = exp4exp5.features, 
                                    verbose = FALSE)
#Next, identify anchors and integrate the datasets. Commands are identical to the standard workflow, 
#but make sure to set normalization.method = 'SCT':
exp4exp5.anchors <- FindIntegrationAnchors(object.list = exp4exp5.list, normalization.method = "SCT", 
                                           anchor.features = exp4exp5.features, verbose = FALSE)
sr.icb <- IntegrateData(anchorset = exp4exp5.anchors, normalization.method = "SCT", 
                              verbose = FALSE)
#run PCA and umap, find clusters
sr.icb <- RunPCA(sr.icb, verbose = FALSE)
ElbowPlot(sr.icb, ndims= 50)
sr.icb <- RunUMAP(sr.icb, dims = 1:40, verbose = FALSE)
sr.icb <- FindNeighbors(sr.icb, dims = 1:40, verbose = FALSE)
sr.icb <- FindClusters(sr.icb, verbose = FALSE)

#make plots
DimPlot(sr.icb, pt=0.3, label = FALSE, group.by = "experiment") 

sr.icb$icb.treatment <- ifelse((sr.icb$treatment == "C"), "C",
                               ifelse(sr.icb$treatment == "P", "P",
                                      ifelse(sr.icb$treatment == "PBS", "PBS",
                                             ifelse(sr.icb$treatment == "PC", "PC", NA))))

sr.icb$pre_postadj <- ifelse((sr.icb$pre_post == "POST"), "POST", 
                           ifelse(sr.icb$pre_post == "PRE", "PRE", NA))

sr.icb$treat_pre_post <- paste(sr.icb$icb.treatment, sr.icb$pre_postadj)

DimPlot(sr.icb, group.by  = "treat_pre_post", na.value = "grey90", pt.size = .5, order = c("PBS POST", "PBS PRE")) +
  scale_color_manual(values=c("gray90", "gray90","gray90", "gray90", "gray90", "gray90", "goldenrod2", "gray90")) +
  ggtitle("ICB - PBS")

sr.icb$expansions <- sr.icb$Freq - 1
sr.icb$expansionfreq <- sr.icb$expansions / sr.icb$totalDetectedTCR
sr.icb$expansionfreq[sr.icb$expansionfreq == 0] <- NA

FeaturePlot(sr.icb, features = "expansionfreq")+ ggtitle("Expansion Freq") + 
  scale_color_gradient(low="lightgoldenrod2",high="firebrick", na.value = "grey90") 

sr.icb$expansionfreq.pc.pre <- ifelse((sr.icb$treat_pre_post == "PC PRE"), sr.icb$expansionfreq, NA)
sr.icb$expansionfreq.pc.post <- ifelse((sr.icb$treat_pre_post == "PC POST"), sr.icb$expansionfreq, NA)
sr.icb$expansionfreq.p.pre <- ifelse((sr.icb$treat_pre_post == "P PRE"), sr.icb$expansionfreq, NA)
sr.icb$expansionfreq.p.post <- ifelse((sr.icb$treat_pre_post == "P POST"), sr.icb$expansionfreq, NA)
sr.icb$expansionfreq.c.pre <- ifelse((sr.icb$treat_pre_post == "C PRE"), sr.icb$expansionfreq, NA)
sr.icb$expansionfreq.c.post <- ifelse((sr.icb$treat_pre_post == "C POST"), sr.icb$expansionfreq, NA)
sr.icb$expansionfreq.pbs.pre <- ifelse((sr.icb$treat_pre_post == "PBS PRE"), sr.icb$expansionfreq, NA)
sr.icb$expansionfreq.pbs.post <- ifelse((sr.icb$treat_pre_post == "PBS POST"), sr.icb$expansionfreq, NA)

FeaturePlot(sr.icb, features = "expansionfreq.pbs.post")+ ggtitle("Expansion Freq PBS Post") + 
  scale_color_gradient(low="lightgoldenrod2",high="firebrick", na.value = "grey90") 

## SHARED
tcr.cluster <- select(sr.icb@meta.data, mouse, cdr3, percentClone)
tcr.cluster$barcode <- row.names(tcr.cluster)

cdr3 <- levels(as.factor(tcr.cluster$cdr3))

for(i in levels(as.factor(sr.icb$mouse))){
  tcr.cluster[[i]]<-tcr.cluster[i,2] 
  tcr.cluster[,ncol(tcr.cluster)] <- ifelse((tcr.cluster$mouse == i), tcr.cluster$percentClone, NA)
}

cluster.sample <- select(tcr.cluster, cdr3,C13M1POST,C13M1PRE,C13M2POST,C13M2PRE,C13M3POST,C13M3PRE,C14M1POST,C14M1PRE,C14M2POST,C14M2PRE,    
C14M3POST,C14M3PRE,C14M5POST,C14M5PRE,C15M1POST,C15M1PRE,C15M4POST,C15M4PRE,C15M5POST,C15M5PRE,C16M3POST,C16M3PRE,C16M4POST,C16M4PRE)

cdr3.cluster <- as.data.frame(setDT(cluster.sample)[, lapply(.SD, function(x) sort(x)[1L]), by = .(cdr3)])

to.merge <- select(tcr.cluster, cdr3, barcode)

tcr.merged.clusters <- merge(cdr3.cluster, to.merge, by = 'cdr3', all = TRUE)

rownames(tcr.merged.clusters) <- tcr.merged.clusters$barcode
tcr.merged.clusters$barcode <- NULL

tcr.merged.clusters[,2:ncol(tcr.merged.clusters)][is.na(tcr.merged.clusters[,2:ncol(tcr.merged.clusters)])] <- 0

sr.icb <- AddMetaData(sr.icb, metadata = tcr.merged.clusters)

sr.icb$shared <- ifelse((sr.icb$C13M1POST > 0 & sr.icb$C13M1PRE > 0), "Shared",
                      ifelse(sr.icb$C13M2POST > 0 & sr.icb$C13M2PRE > 0, "Shared",
                      ifelse(sr.icb$C13M3POST > 0 & sr.icb$C13M3PRE > 0, "Shared",
                      ifelse(sr.icb$C14M1POST > 0 & sr.icb$C14M1PRE > 0, "Shared",
                      ifelse(sr.icb$C14M2POST > 0 & sr.icb$C14M2PRE > 0, "Shared",
                      ifelse(sr.icb$C14M3POST > 0 & sr.icb$C14M3PRE > 0, "Shared",
                      ifelse(sr.icb$C14M5POST > 0 & sr.icb$C14M5PRE > 0, "Shared",
                      ifelse(sr.icb$C15M1POST > 0 & sr.icb$C15M1PRE > 0, "Shared",
                      ifelse(sr.icb$C15M4POST > 0 & sr.icb$C15M4PRE > 0, "Shared",
                      ifelse(sr.icb$C15M5POST > 0 & sr.icb$C15M5PRE > 0, "Shared",
                      ifelse(sr.icb$C16M3POST > 0 & sr.icb$C16M3PRE > 0, "Shared",
                      ifelse(sr.icb$C16M4POST > 0 & sr.icb$C16M4PRE > 0, "Shared", NA))))))))))))       

# Shared Expanded  
sr.icb$shared.day <- ifelse((sr.icb$C13M1POST > 0 & sr.icb$C13M1PRE > 0 & sr.icb$mouse == "C13M1PRE"), "Shared PC Pre",
                 ifelse(sr.icb$C13M1POST > 0 & sr.icb$C13M1PRE > 0 & sr.icb$mouse == "C13M1POST", "Shared PC Post", 
                 ifelse(sr.icb$C13M2POST > 0 & sr.icb$C13M2PRE > 0 & sr.icb$mouse == "C13M2PRE", "Shared PC Pre",
                 ifelse(sr.icb$C13M2POST > 0 & sr.icb$C13M2PRE > 0 & sr.icb$mouse == "C13M2POST", "Shared PC Post",
                 ifelse(sr.icb$C13M3POST > 0 & sr.icb$C13M3PRE > 0 & sr.icb$mouse == "C13M3PRE", "Shared PBS Pre",
                 ifelse(sr.icb$C13M3POST > 0 & sr.icb$C13M3PRE > 0 & sr.icb$mouse == "C13M3POST", "Shared PBS Post",
                 ifelse(sr.icb$C14M1POST > 0 & sr.icb$C14M1PRE > 0 & sr.icb$mouse == "C14M1PRE", "Shared C Pre",
                 ifelse(sr.icb$C14M1POST > 0 & sr.icb$C14M1PRE > 0 & sr.icb$mouse == "C14M1POST", "Shared C Post",
                 ifelse(sr.icb$C14M2POST > 0 & sr.icb$C14M2PRE > 0 & sr.icb$mouse == "C14M2PRE", "Shared PC Pre",
                 ifelse(sr.icb$C14M2POST > 0 & sr.icb$C14M2PRE > 0 & sr.icb$mouse == "C14M2POST", "Shared PC Post", 
                 ifelse(sr.icb$C14M3POST > 0 & sr.icb$C14M3PRE > 0 & sr.icb$mouse == "C14M3PRE", "Shared PC Pre",
                 ifelse(sr.icb$C14M3POST > 0 & sr.icb$C14M3PRE > 0 & sr.icb$mouse == "C14M3POST", "Shared PC Post",
                 ifelse(sr.icb$C14M5POST > 0 & sr.icb$C14M5PRE > 0 & sr.icb$mouse == "C14M5PRE", "Shared PBS Pre",
                 ifelse(sr.icb$C14M5POST > 0 & sr.icb$C14M5PRE > 0 & sr.icb$mouse == "C14M5POST", "Shared PBS Post",
                 ifelse(sr.icb$C15M1POST > 0 & sr.icb$C15M1PRE > 0 & sr.icb$mouse == "C15M1PRE", "Shared C Pre",
                 ifelse(sr.icb$C15M1POST > 0 & sr.icb$C15M1PRE > 0 & sr.icb$mouse == "C15M1POST", "Shared C Post",
                 ifelse(sr.icb$C15M4POST > 0 & sr.icb$C15M4PRE > 0 & sr.icb$mouse == "C15M4PRE", "Shared PBS Pre",
                 ifelse(sr.icb$C15M4POST > 0 & sr.icb$C15M4PRE > 0 & sr.icb$mouse == "C15M4POST", "Shared PBS Post",
                 ifelse(sr.icb$C15M5POST > 0 & sr.icb$C15M5PRE > 0 & sr.icb$mouse == "C15M5PRE", "Shared P Pre",
                 ifelse(sr.icb$C15M5POST > 0 & sr.icb$C15M5PRE > 0 & sr.icb$mouse == "C15M5POST", "Shared P Post",
                 ifelse(sr.icb$C16M3POST > 0 & sr.icb$C16M3PRE > 0 & sr.icb$mouse == "C16M3PRE", "Shared P Pre",
                 ifelse(sr.icb$C16M3POST > 0 & sr.icb$C16M3PRE > 0 & sr.icb$mouse == "C16M3POST", "Shared P Post",
                 ifelse(sr.icb$C16M4POST > 0 & sr.icb$C16M4PRE > 0 & sr.icb$mouse == "C16M4PRE", "Shared PC Pre",
                 ifelse(sr.icb$C16M4POST > 0 & sr.icb$C16M4PRE > 0 & sr.icb$mouse == "C16M4POST", "Shared PC Post", NA))))))))))))))))))))))))

DimPlot(sr.icb, group.by = "shared.day",pt.size = 0.5, na.value = "grey90", order = c("Shared PBS Pre", "Shared PBS Post")) + 
  ggtitle("Shared PBS") + theme(plot.title = element_text(hjust = 0.5)) + 
  scale_color_manual(values=c("gray90", "gray90", "gray90", "gray90", "gray90", "gray90", "darkorange", "purple4"))

# Shared by Mouse
tcr.cluster <- select(sr.icb@meta.data, mouse, cdr3, percentClone)
tcr.cluster$barcode <- row.names(tcr.cluster)

cdr3 <- levels(as.factor(tcr.cluster$cdr3))

for(i in levels(as.factor(sr.icb$mouse))){
  tcr.cluster[[i]]<-tcr.cluster[i,2] 
  tcr.cluster[,ncol(tcr.cluster)] <- ifelse((tcr.cluster$mouse == i), tcr.cluster$percentClone, NA)
}

cluster.sample <- select(tcr.cluster, cdr3,"C13M1POST","C13M1PRE","C13M2POST","C13M2PRE","C13M3POST","C13M3PRE","C14M2POST","C14M2PRE")

cdr3.cluster <- as.data.frame(setDT(cluster.sample)[, lapply(.SD, function(x) sort(x)[1L]), by = .(cdr3)])

to.merge <- select(tcr.cluster, cdr3, barcode)

tcr.merged.clusters <- merge(cdr3.cluster, to.merge, by = 'cdr3', all = TRUE)

rownames(tcr.merged.clusters) <- tcr.merged.clusters$barcode
tcr.merged.clusters$barcode <- NULL

tcr.merged.clusters[,2:ncol(tcr.merged.clusters)][is.na(tcr.merged.clusters[,2:ncol(tcr.merged.clusters)])] <- 0

sr.icb <- AddMetaData(sr.icb, metadata = tcr.merged.clusters)

View(sr.icb@meta.data)

sr.icb$shared <- ifelse((sr.icb$C13M1PRE > 0 & sr.icb$C13M1POST > 0), "Shared",
                        ifelse(sr.icb$C13M2PRE > 0 & sr.icb$C13M2POST > 0, "Shared",
                               ifelse(sr.icb$C13M3PRE > 0 & sr.icb$C13M3POST > 0, "Shared",
                                      ifelse(sr.icb$C14M2PRE > 0 & sr.icb$C14M2POST > 0, "Shared", NA))))
DimPlot(sr.icb, group.by = "shared", pt.size = 0.5, na.value = "grey90") + 
  ggtitle("Shared") + theme(plot.title = element_text(hjust = 0.5)) + 
  scale_color_manual(values=c("gray90", "gray90", "gray90", "gray90", "gray90", "gray90",
                              "gray90", "gray90", "gray90", "darkorange", "purple4"))

sr.icb$mouse.day <- as.factor(sr.icb$mouse)

sr.icb$shared.day <- ifelse((sr.icb$C13M1PRE > 0 & sr.icb$C13M1POST > 0 & sr.icb$mouse.day == "C13M1PRE"), "Shared PC PRE",
                            ifelse(sr.icb$C13M1PRE > 0 & sr.icb$C13M1POST > 0 & sr.icb$mouse.day == "C13M1POST", "Shared PC POST", 
                                   ifelse(sr.icb$C13M2PRE > 0 & sr.icb$C13M2POST > 0 & sr.icb$mouse.day == "C13M2PRE", "Shared PC PRE",
                                          ifelse(sr.icb$C13M2PRE > 0 & sr.icb$C13M2POST > 0 & sr.icb$mouse.day == "C13M2POST", "Shared PC POST",
                                                 ifelse(sr.icb$C13M3PRE > 0 & sr.icb$C13M3POST > 0 & sr.icb$mouse.day == "C13M3PRE", "Shared PBS PRE",
                                                        ifelse(sr.icb$C13M3PRE > 0 & sr.icb$C13M3POST > 0 & sr.icb$mouse.day == "C13M3POST", "Shared PBS POST", 
                                                               ifelse(sr.icb$C14M2PRE > 0 & sr.icb$C14M2POST > 0 & sr.icb$mouse.day == "C14M2PRE", "Shared PC PRE",
                                                                      ifelse(sr.icb$C14M2PRE > 0 & sr.icb$C14M2POST > 0 & sr.icb$mouse.day == "C14M2POST", "Shared PC POST", "Unshared"))))))))

DimPlot(sr.icb, group.by = "shared.day", pt.size = 0.5, na.value = "grey90") + 
  ggtitle("Shared") + theme(plot.title = element_text(hjust = 0.5)) + 
  scale_color_manual(values=c("gray90", "gray90", "gray90", "gray90", "gray90", "gray90",
                              "gray90", "gray90", "gray90", "darkorange", "purple4"))                     
# EXPANSION or CONTRACTION
expand.threshold <- 0.05
contract.threshold <- 0.05

sr.icb$expand <- ifelse((sr.icb$C13M1POST > sr.icb$C13M1PRE + expand.threshold), "Expand",
                        ifelse(sr.icb$C13M2POST > sr.icb$C13M2PRE + expand.threshold, "Expand",
                               ifelse(sr.icb$C13M3POST > sr.icb$C13M3PRE + expand.threshold, "Expand", 
                                      ifelse(sr.icb$C14M2POST > sr.icb$C14M2PRE + expand.threshold, "Expand", "Other"))))

sr.icb$contract <-  ifelse((sr.icb$C13M1PRE > sr.icb$C13M1POST + contract.threshold), "Contract",
                           ifelse(sr.icb$C13M2PRE > sr.icb$C13M2POST + contract.threshold, "Contact",
                                  ifelse(sr.icb$C13M3PRE > sr.icb$C13M3POST + contract.threshold, "Contract", 
                                         ifelse(sr.icb$C14M2PRE > sr.icb$C14M2POST + contract.threshold, "Contract", "Other"))))

DimPlot(sr.icb, group.by = "expand", pt.size = 0.5, na.value = "grey90") + 
  ggtitle("Expansion") + theme(plot.title = element_text(hjust = 0.5)) + 
  scale_color_manual(values=c("gray90", "gray90", "gray90", "gray90", "gray90", "gray90",
                              "gray90", "gray90", "gray90", "darkorange", "purple4"))  

# Merge our two experiment types - Create "Experiment.merge" column to correctly put pieces back together
sr.turnover$experiment.merge <- "exp4"
sr.icb$experiment.merge <- sr.icb$experiment

options(future.globals.maxSize = 4000 * 1024^2)

#merge 

exp5[["percentMT"]]<-PercentageFeatureSet(exp5, pattern = "^mt-")
FeatureScatter(exp5, feature1 = "nCount_RNA", feature2 = "percentMT") + geom_hline(yintercept = 25) + geom_vline(xintercept = 38000)


sr.all <- merge(sr.turnover, y = sr.icb, add.cell.ids = c("turnover", "icb"), project = "sr.all")

all.list <- SplitObject(sr.all, split.by = "experiment.merge")
all.list <- all.list[c("exp4", "exp5")]
for (i in 1:length(all.list)) {
  all.list[[i]] <- SCTransform(all.list[[i]], verbose = FALSE)
}

all.features <- SelectIntegrationFeatures(object.list = all.list, nfeatures = 3000)
all.list <- PrepSCTIntegration(object.list = all.list, anchor.features = all.features, 
                                    verbose = FALSE)
all.anchors <- FindIntegrationAnchors(object.list = all.list, normalization.method = "SCT", 
                                           anchor.features = all.features, verbose = FALSE)
sr.combined <- IntegrateData(anchorset = all.anchors, normalization.method = "SCT", 
                                     verbose = FALSE)

sr.combined <- RunPCA(sr.combined, verbose = FALSE)
ElbowPlot(sr.combined, ndims= 50)
sr.combined <- RunUMAP(sr.combined, dims = 1:30)
sr.combined <- FindNeighbors(sr.combined, dims = 1:30, verbose = FALSE)
sr.combined <- FindClusters(sr.combined, verbose = FALSE, resolution = )
DimPlot(sr.combined, label = TRUE, pt.size = 0.35)

DimPlot(sr.combined, group.by = "experiment", pt=, order = c("exp5")) + 
  scale_color_manual(values = c("gray90", "gray90", "indianred3"))

DimPlot(sr.combined, group.by = "treat_pre_post", pt=, order = c("PC PRE", "PC POST")) + 
  scale_color_manual(values = c("gray90", "gray90", "gray90","gray90", "gray90", "gray90","navy", "darkgoldenrod2", "gray90"))

DimPlot(sr.combined, group.by = "day.general", pt=, order = c("D12")) + 
  scale_color_manual(values = c("gray90", "gray90", "gray90","gray90", "purple4", "gray90","navy", "darkgoldenrod2", "gray90", "gray90"))

FeaturePlot(sr.combined , features = c("Cd4"), pt=, min.cutoff =0, max.cutoff =5)
FeaturePlot(sr.combined, features = c("adt_CD4"), ncol= 1, min.cutoff = "q10", max.cutoff = "q95")

DimPlot(sr.combined, group.by = "shared.day", pt.size = , na.value = "grey90", order=c("Shared D12_1", "Shared D12_2")) + 
  ggtitle("Shared D12") + theme(plot.title = element_text(hjust = 0.5)) + 
  scale_color_manual(values=c("gray90", "gray90", "gray90", "gray90", "gray90", "gray90",
                              "gray90", "gray90", "gray90", "gray90", "gray90","gray90", "gray90", "gray90", "gray90", "gray90", "gray90",
                              "purple4", "darkorange2", "gray90", "gray90", "gray90")) 

FeaturePlot(sr.combined, features = "percentClone") + ggtitle("Clonality") + 
  scale_color_gradient(low="lightgoldenrod2",high="firebrick", na.value = "grey90") 

sr.combined$expansions <- sr.combined$Freq - 1
sr.combined$expansionfreq <- sr.combined$expansions / sr.combined$totalDetectedTCR
sr.combined$expansionfreq[sr.combined$expansionfreq == 0] <- NA

FeaturePlot(sr.combined, features = "expansionfreq") + ggtitle("Expansion Freq") + 
  scale_color_gradient(low="lightgoldenrod2",high="firebrick", na.value = "grey90") 

sr.combined$expansionfreq.cd4 <- ifelse((sr.combined$label == "CD4"), sr.combined$expansionfreq, NA)
sr.combined$expansionfreq.cd8 <- ifelse((sr.combined$label == "CD8"), sr.combined$expansionfreq, NA)

FeaturePlot(sr.combined, features = "expansionfreq.cd4") + ggtitle("Expansion Freq CD4") + 
  scale_color_gradient(low="lightgoldenrod2",high="firebrick", na.value = "grey90") 

FeaturePlot(sr.combined, features = "expansionfreq.cd8") + ggtitle("Expansion Freq CD8") + 
  scale_color_gradient(low="lightgoldenrod2",high="firebrick", na.value = "grey90") 

#Look for CD4/CD8
#Retry integration (Normalize by experiment)
options(future.globals.maxSize = 4000 * 1024^2)

all.list <- SplitObject(sr.combined, split.by = "integrate")
all.list <- all.list[c("exp4 icb", "exp5 icb", "turnover")]
for (i in 1:length(all.list)) {
  all.list[[i]] <- SCTransform(all.list[[i]], verbose = FALSE)
}

all.features <- SelectIntegrationFeatures(object.list = all.list, nfeatures = 3000)
all.list <- PrepSCTIntegration(object.list = all.list, anchor.features = all.features, 
                               verbose = FALSE)
all.anchors <- FindIntegrationAnchors(object.list = all.list, normalization.method = "SCT", 
                                      anchor.features = all.features, verbose = FALSE)
sr.merged <- IntegrateData(anchorset = all.anchors, normalization.method = "SCT", 
                             verbose = FALSE)

sr.merged.cd4 <- subset(sr.merged, subset = label == "CD4")
sr.merged.cd8 <- subset(sr.merged, subset = label == "CD8")

sr.merged.cd4 <- RunPCA(sr.merged.cd4, verbose = FALSE)
ElbowPlot(sr.merged.cd4, ndims= 50)
sr.merged.cd4 <- RunUMAP(sr.merged.cd4, dims = 1:25)
sr.merged.cd4 <- FindNeighbors(sr.merged.cd4, dims = 1:25, verbose = FALSE)
sr.merged.cd4 <- FindClusters(sr.merged.cd4, verbose = FALSE, resolution = )
DimPlot(sr.merged.cd4, label = TRUE, pt.size = ) + ggtitle("CD4 Clusters")

sr.merged.cd8 <- RunPCA(sr.merged.cd8, verbose = FALSE)
ElbowPlot(sr.merged.cd8, ndims= 50)
sr.merged.cd8 <- RunUMAP(sr.merged.cd8, dims = 1:25)
sr.merged.cd8 <- FindNeighbors(sr.merged.cd8, dims = 1:25, verbose = FALSE)
sr.merged.cd8 <- FindClusters(sr.merged.cd8, verbose = FALSE, resolution = )
DimPlot(sr.merged.cd8, label = TRUE, pt.size = ) + ggtitle("CD8 Clusters")

plot2 <- DimPlot(sr.merged.cd8, group.by = "experiment", pt=, order = c("exp4")) + 
  scale_color_manual(values = c("deepskyblue2", "gray90", "indianred3")) + ggtitle("CD8 Separate Merge")

plot1 <- FeatureScatter(sr.merged, feature1 = "nCount_RNA", feature2 = "percentMT") + geom_hline(yintercept = 25) + geom_vline(xintercept = 38000)
plot2 <- FeatureScatter(sr.merged, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_hline(yintercept = 200) + geom_hline(yintercept = 0) + geom_vline(xintercept = 37000)
CombinePlots(plots = list(plot1, plot2))

sr.combined$remove <- ifelse(sr.combined$experiment.merge == "exp5", 0, sr.combined$percentMT)
  
sr.combined <- subset(x=sr.combined, label != "negative" & stimvnonstim != "negative" & nFeature_RNA < 5250 & remove < 20 & nCount_RNA<37000 & nFeature_RNA > 200)

sr.combined$test <- ifelse((sr.combined$day.general == "D0"), "turnover", 
                           ifelse(sr.combined$day.general == "D3", "turnover", 
                                  ifelse(sr.combined$day.general == "D6", "turnover", 
                                         ifelse(sr.combined$day.general == "D9", "turnover", 
                                                ifelse(sr.combined$day.general == "D12", "turnover", "icb")))))

sr.combined.cd4 <- subset(sr.combined, subset = label == "CD4")
sr.combined.cd8 <- subset(sr.combined, subset = label == "CD8")

sr.combined.cd4 <- RunPCA(sr.combined.cd4, verbose = FALSE)
ElbowPlot(sr.combined.cd4, ndims= 50)
sr.combined.cd4 <- RunUMAP(sr.combined.cd4, dims = 1:25)
sr.combined.cd4 <- FindNeighbors(sr.combined.cd4, dims = 1:25, verbose = FALSE)
sr.combined.cd4 <- FindClusters(sr.combined.cd4, verbose = FALSE, resolution = )
DimPlot(sr.merged.cd4, label = TRUE, pt.size = ) + ggtitle("CD4 Clusters")

sr.combined.cd8 <- RunPCA(sr.combined.cd8, verbose = FALSE)
ElbowPlot(sr.combined.cd8, ndims= 50)
sr.combined.cd8 <- RunUMAP(sr.combined.cd8, dims = 1:25)
sr.combined.cd8 <- FindNeighbors(sr.combined.cd8, dims = 1:25, verbose = FALSE)
sr.combined.cd8 <- FindClusters(sr.combined.cd8, verbose = FALSE, resolution = )
DimPlot(sr.merged.cd8, label = TRUE, pt.size = ) + ggtitle("CD8 Clusters")

DimPlot(sr.merged.cd4, group.by = "experiment", pt=, order = c("exp4")) + 
  scale_color_manual(values = c("deepskyblue2", "gray90", "indianred3")) + ggtitle("CD8 Combined Merge")

DimPlot(sr.merged.cd8, group.by = "treat.full", pt=, order = c("PBS PRE", "PBS POST")) + 
  scale_color_manual(values = c("gray90", "gray90", "gray90","gray90", "gray90", "gray90","navy", "darkgoldenrod2", "gray90")) +
  ggtitle("CD8 - PBS Pre vs. Post")

DimPlot(sr.merged.cd8, group.by = "day.general", pt=, order = c("D0")) + 
  scale_color_manual(values = c("gray90", "gray90", "gray90","gray90", "purple4")) +
  ggtitle("CD8 - D0")

DimPlot(sr.merged.cd8, group.by = "stimvnonstim", pt=, order = ) + 
  scale_color_manual(values = c("gray90", "firebrick", "lightgreen","gray90", "purple4", "gray90","navy", "darkgoldenrod2", "gray90", "gray90"))


DefaultAssay(sr.merged.cd4) <- "RNA"
DefaultAssay(sr.merged.cd8) <- "RNA"
sr.merged.cd4 <- NormalizeData(sr.merged.cd4, verbose = FALSE)
sr.merged.cd8 <- NormalizeData(sr.merged.cd8, verbose = FALSE)

FeaturePlot(sr.merged.cd4 , features = c("Cd40lg"), pt=, min.cutoff =, max.cutoff =)
FeaturePlot(sr.merged.cd4, features = c("adt_CD4"), ncol= 1, min.cutoff = "q10", max.cutoff = "q95")

FeaturePlot(sr.merged.cd4, features = "expansionfreq") + ggtitle("CD4 Expansion Freq") + 
  scale_color_gradient(low="lightgoldenrod2",high="firebrick", na.value = "grey90") 

DimPlot(sr.merged.cd8, group.by = "shared.day.ns", pt.size = , na.value = "grey90", order=c("Shared D12_1", "Shared D12_2")) + 
  ggtitle("CD8 Shared D12") + theme(plot.title = element_text(hjust = 0.5)) + 
  scale_color_manual(values=c("gray90", "gray90", "gray90", "gray90", "gray90", "gray90",
                              "gray90", "gray90", "gray90", "gray90", "gray90","gray90", "gray90", "gray90", "gray90", "gray90", "gray90","gray90",
                              "purple4", "darkorange2", "gray90", "gray90", "gray90")) 

clusterMarkers.cd8 <- FindAllMarkers(sr.merged.cd8, only.pos = TRUE, min.pct = 0.15, logfc.threshold = 0.15, test.use = "MAST")
topMarkers.cd8<-clusterMarkers.cd8 %>% group_by(cluster) %>% top_n(n = 25, wt = avg_logFC)

shortMarkers.cd8<-clusterMarkers.cd8 %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
shortM.cd4<-shortMarkers.cd4[,c(6:7)]
shortM.cd4$rank<-c(1:20)
head(shortM %>% spread(rank, gene))

cluster.11.markers.cd4 <- FindMarkers(sr.merged.cd4, ident.1 =11, ident.2 = , logfc.threshold = 0.1,min.diff.pct = 0.1, test.use = "MAST")

write.csv(cluster.1v2.markers.cd4, file="TregC1C2.cd4.csv")

FeaturePlot(sr.merged.cd4 , features = c("Sell"), pt=, min.cutoff =0, max.cutoff =5)
FeaturePlot(sr.merged.cd4, features = c("adt_CD62L"), ncol= 1, min.cutoff = "q10", max.cutoff = "q95")

# END #


