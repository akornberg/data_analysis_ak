#Expected TCR Sharing Between Clusters - Adam Kornberg 
library(ggplot2)
library(Seurat)
library(tidyverse)
library(data.table)
library(corrplot)

DimPlot(sr.merged.cd4, label = TRUE, pt.size =0.5) 

#Improve CDR3B Count
tcr.count <- sr.merged.cd4@meta.data %>% group_by(mouse,cdr3b) %>% summarise(count = n())
sr.merged.cd4@meta.data <- merge(sr.merged.cd4@meta.data, tcr.count)

tcr.count <- sr.merged.cd8@meta.data %>% group_by(mouse,cdr3b) %>% summarise(count = n())
sr.merged.cd8@meta.data <- merge(sr.merged.cd8@meta.data, tcr.count)

# tcrs is relavant TCR data frame from Samhita's merged CD4 dataset. All NA's have been removed.
tcrs <- subset(sr.merged.cd4@meta.data, subset = cdr3.count >= 0)

# count is the frequency of CDR3Bs (NOT BY SAMPLE)
count <- select(tcrs,cdr3b, count, mouse)

# Total # of TCRs by each sample or each cluster
tcr.by.sample <- tcrs %>% group_by(sample) %>% summarise(sample.count = n())
tcr.by.cluster <- tcrs %>% group_by(seurat_clusters) %>% summarise(cluster.count = n())
tcr.by.both <- tcrs %>% group_by(seurat_clusters, sample) %>% summarise(cluster.count = n())

# CDR3B counted by cluster or by sample
cdr3.by.cluster <- tcrs %>% group_by(cdr3b, seurat_clusters) %>% summarise(cdr3.cluster.count=n())
cdr3.by.sample <- tcrs %>% group_by(cdr3b, sample) %>% summarise(cdr3.sample.count=n())
cdr3.by.both <- tcrs %>% group_by(cdr3b, sample, seurat_clusters) %>% summarise(cdr3.sample.count=n())

# TCR "Cluster Frequency", is the sample-dependent frequency of a CDR3B in a given cluster. Here, each CDR3B is given it's own ID (sample dependent == id, sample-independent == cdr3b.id, row.name == ID of the row to re-bind back later)
tcr.full <- merge(tcr.by.both, cdr3.by.both)
tcr.full <- merge(tcr.full, tcr.by.sample, by="sample")
tcr.full$row.name <- row.names(tcr.full)
tcr.full$cluster.freq <- (tcr.full$cdr3.sample.count/tcr.full$cluster.count)*100
tcr.full$cdr3b.sample <- paste(tcr.full$sample, tcr.full$cdr3b)
tcr.full <- transform(tcr.full, cdr3b.id=match(cdr3b, unique(cdr3b)))
tcr.full <- transform(tcr.full, id=match(cdr3b.sample, unique(cdr3b.sample)))

# Use this function to determine a few important variables, all of this is sample-dependent:
## max => Most expanded CDR3B of any cluster (sample-dependent)
## ID.cluster.sum => Total cells from a given sample in clusters featuring shared CDR3B of interest
## Denominator => ID.cluster.sum - actual # of cells in a cluster
## TotalCDR3B => Total # of shared CDR3B 
## Numerator => TotalCDR3B - cdr3.sample.count
## Expected.sharing => Numerator/Denominator . Calculated for each clone individually
sharing <- data.frame()
for(i in levels(as.factor(tcr.full$id))){
  workingDF <- subset(tcr.full, id == i)
  workingDF$max <- max(workingDF$cdr3.sample.count)
  workingDF$max.id <- ifelse((workingDF$max == workingDF$cdr3.sample.count), workingDF$cluster.count, 0)
  workingDF$max.cluster.count <- ifelse((workingDF$max.id == 0), max(workingDF$max.id), workingDF$max.id)
  workingDF$ID.cluster.sum <- sum(workingDF$cluster.count)
  workingDF$denominator <- workingDF$sample.count - workingDF$max.cluster.count
  workingDF$totalCDR3B <- sum(workingDF$cdr3.sample.count)
  workingDF$numerator <- workingDF$totalCDR3B - workingDF$max
  workingDF$expected.sharing <- (workingDF$numerator/workingDF$denominator)*100
  workingDF$difference <- workingDF$cluster.freq - workingDF$expected.sharing
  resultsDF<-select(workingDF,row.name,id, cdr3.sample.count, max,max.id,max.cluster.count, ID.cluster.sum,denominator, totalCDR3B, numerator, expected.sharing, difference)
  sharing <- rbind(sharing, resultsDF)
  sharing$row.name <- as.numeric(sharing$row.name)
  sharing <- sharing[order(sharing[,1]),]
}
sharing.bind <- select(sharing, max,max.id, max.cluster.count, ID.cluster.sum, totalCDR3B, numerator, denominator, expected.sharing, difference)
tcr.full <- cbind(tcr.full, sharing.bind)

## Correlation Analysis
ggplot(tcr.full,aes(x=cluster.count)) + geom_histogram(binwidth=0.05) +
  scale_x_continuous(trans='log10',limits=c(1,500)) 
ggplot(tcr.full,aes(x=cluster.freq)) + geom_histogram(binwidth=0.05) +
  scale_x_continuous(trans='log10',limits=c(1,1000)) 

#Store TCR.Full for later
tcr.store <- tcr.full
#Choose Cutoffs for Dataset (Based of Cluster Size or Appearances in a cluster) SPLIT UP FILES
cutoff <- 10
tcr.full <- subset(tcr.full, subset = cluster.count >= cutoff)
#Create Matrix for Correlation Analysis
tcr.general <- select(tcr.full, seurat_clusters, cdr3b,cdr3b.sample, cluster.freq, cluster.count, id)
tcr.general$row.name <- row.names(tcr.general)

cdr3 <- levels(as.factor(tcr.general$cdr3b.sample))

for(i in unique(as.factor(tcr.store$seurat_clusters))){
  tcr.general[[i]]<-tcr.general[i,2] 
  tcr.general[ncol(tcr.general)] <- ifelse((tcr.general$seurat_clusters == i), tcr.general$cluster.freq, NA)
}

#ADJUST THE AMOUNT OF CLUSTERS!
tcr.sample <- select(tcr.general, cdr3b.sample, '0', '1', '2', '3', '4', '5', '6', '7', '8', '9','10','11','12','13','14','15')

cdr3.sample <- as.data.frame(setDT(tcr.sample)[, lapply(.SD, function(x) sort(x)[1L]), by = .(cdr3b.sample)])

to.merge <- select(tcr.general, cdr3b.sample, id)

tcr.merged <- merge(cdr3.sample, to.merge, by = 'cdr3b.sample', all = TRUE)

tcr.merged[is.na(tcr.merged)] <- 0
colnames(tcr.merged) <- c("cdr3b.sample", 'Cluster 0', 'Cluster 1','Cluster 2','Cluster 3','Cluster 4','Cluster 5',
                          'Cluster 6','Cluster 7','Cluster 8','Cluster 9','Cluster 10','Cluster 11','Cluster 12',
                          'Cluster 13','Cluster 14','Cluster 15',"id")
# Preparation for Correlation Analysis. Remove all extra columns + columns with all 0's.
tcr.merged$id <- NULL
tcr.merged <- unique(tcr.merged)
row.names(tcr.merged) <- tcr.merged$cdr3b.sample
tcr.merged$cdr3b.sample <- NULL
tcr.merged <- as.matrix(tcr.merged)

# Correlation Analysis. Pearson or spearman can be used. 
c <- cor(tcr.merged, method = c("pearson"))

corrplot(c, method = c("circle"), type = "upper")

#Options for Corrplot are below:
# "circle","square", "ellipse", "number", "pie", "shade" and "color"

## EXAMPLE
CTCSAGTGGYEQYF
test <- subset(tcr.full, subset = id == 296)

## NOT NECESSARY BUT TO MERGE TCR.MERGED/TCR.FULL FOLLOW BELOW. MUST PERFORM WITH NO CUTOFFS.
tcr.merged <- tcr.merged[order(tcr.merged$id),]
tcr.full <- tcr.full[order(tcr.full$id),]
tcr.full <- cbind(tcr.full, tcr.merged)
tcr.full$row.name <- as.numeric(tcr.full$row.name)
tcr.full <- tcr.full[order(tcr.full$row.name),]
