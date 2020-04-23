# Expected TCR Sharing Between Clusters - Adam Kornberg 
library(ggplot2)
library(Seurat)
library(tidyverse)
library(data.table)
library(corrplot)

# Import Data 
import <- read.csv(file.choose())

# Remove all cells without a TCR
tcrs <- subset(import, subset = TRBCDR3count >= 0)

# Improve CDR3B Count (No need for second line since I don't have the Seurat Object but refer to old script to merge this back to your Seurat Object if needed)
tcr.count <- tcrs %>% group_by(sample,TRBcdr3) %>% summarise(count = n())
SeuratObject.CD8@meta.data <- merge(SeuratObject.CD8@meta.data, tcr.count)

# Remove the old count, add corrected 'count' variable. 
tcrs$TRBCDR3count <- NULL
tcrs <- merge(tcrs, tcr.count)

# Change column names to match script
colnames(tcrs) <- c("cdr3b", "sample", "seurat_clusters", "count")

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

# Store TCR.Full for later
tcr.store <- tcr.full
# Choose Cutoffs for Dataset (Based of Cluster Size or Appearances in a cluster) SPLIT UP FILES
cutoff <- #INSERT NUMBER HERE
  tcr.full <- subset(tcr.full, subset = cluster.count >= cutoff)

##IF YOU WANT TO SUBSET A SPECIFIC PATIENT DO IT HERE#########

################# #########   ########################3
# Create Matrix for Correlation Analysis
tcr.general <- select(tcr.full, seurat_clusters, cdr3b,cdr3b.sample, cluster.freq, cluster.count, id)
tcr.general$row.name <- row.names(tcr.general)

cdr3 <- levels(as.factor(tcr.general$cdr3b.sample))

for(i in unique(as.factor(tcr.store$seurat_clusters))){
  tcr.general[[i]]<-tcr.general[i,2] 
  tcr.general[ncol(tcr.general)] <- ifelse((tcr.general$seurat_clusters == i), tcr.general$cluster.freq, NA)
}

# ADJUST THE AMOUNT OF CLUSTERS!
tcr.sample <- select(tcr.general, cdr3b.sample, '0', '1', '2', '3', '4', '5', '6', '7', '8')

cdr3.sample <- as.data.frame(setDT(tcr.sample)[, lapply(.SD, function(x) sort(x)[1L]), by = .(cdr3b.sample)])

to.merge <- select(tcr.general, cdr3b.sample, id)

tcr.merged <- merge(cdr3.sample, to.merge, by = 'cdr3b.sample', all = TRUE)

tcr.merged[is.na(tcr.merged)] <- 0
colnames(tcr.merged) <- c("cdr3b.sample", 'Cluster 0', 'Cluster 1','Cluster 2','Cluster 3','Cluster 4','Cluster 5',
                          'Cluster 6','Cluster 7','Cluster 8',"id")
# Preparation for Correlation Analysis. Remove all extra columns + columns with all 0's.
tcr.merged$id <- NULL
tcr.merged <- unique(tcr.merged)
row.names(tcr.merged) <- tcr.merged$cdr3b.sample
tcr.merged$cdr3b.sample <- NULL
tcr.merged <- as.matrix(tcr.merged)

# Correlation Analysis. Pearson or spearman can be used. 
c <- cor(tcr.merged, method = c("pearson"))

corrplot(c, method = c("circle"), type = "upper")

# Options for Corrplot are below:
# "circle","square", "ellipse", "number", "pie", "shade" and "color"

## EXAMPLE
CAISTSGSYYNEQFF
test <- subset(tcr.full, subset = id == 296)

## NOT NECESSARY BUT TO MERGE TCR.MERGED/TCR.FULL FOLLOW BELOW. MUST PERFORM WITH NO CUTOFFS.
tcr.merged <- tcr.merged[order(tcr.merged$id),]
tcr.full <- tcr.full[order(tcr.full$id),]
tcr.full <- cbind(tcr.full, tcr.merged)
tcr.full$row.name <- as.numeric(tcr.full$row.name)
tcr.full <- tcr.full[order(tcr.full$row.name),]
