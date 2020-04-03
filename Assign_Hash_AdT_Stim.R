#Load packages
library(tidyverse)
library(ggplot2)

#Set working directory
setwd("~/Documents/Columbia/Han Lab/Data Analysis/2019_12_IBD/hto/")

#Import and format data
hash <- read.csv(file.choose(), row.names = 1)
hash <- t(hash)
class(hash) <- "numeric"
hash <- as.data.frame(hash)
colnames(hash) <- c("C13M1 PRE","C13M1 POST", "C13M2 PRE","C13M2 POST", "C13M3 PRE",
"C13M3 POST","C14M2 PRE","C14M2 POST", "C21M5 1", "C21M5 2", "C22M5 1", "C22M5 2","C18M4 1",
"C18M4 2", "C19M3 1", "C19M3 2", "C22M1 1","C22M1 2", "stim", "nonstim", "unmapped")

View(hash)
#Create UMItotal column
hash$UMItotal <- rowSums(hash)
ggplot(hash,aes(x=UMItotal)) + geom_histogram(binwidth=0.01) +
  scale_x_continuous(trans='log10',limits=c(1,10000)) + 
  geom_vline(xintercept = 8000, color = "red") + geom_vline(xintercept = 8, color = "red") + ggtitle("Removal of Low/High UMIs")
#normalize data to total number of UMIs per barcode. Hashnormalized values are a % of total UMI counts
hashnormalized <- (hash/hash$UMItotal)*100
View(hashnormalized)
#create subset dataframes (HTO's use stim or sample - AdT's use CD4CD8)
hashstim <- hashnormalized[,c("stim", "nonstim"
  )]
hashsample <- hashnormalized[,c("816A", "816B", "811A", "811B", "794A", "794B", "756A", "756B",
                                "745AB", "742A", "742B", "741A", "741B", "740A", "740B", 
                                "739A", "739B", "730A", "730B", "705A", "705B")]

#Stim vs Nonstim
ggplot(hashstim) + geom_point(mapping=aes(x=nonstim,y=stim), size=1)+
  scale_x_continuous(trans='log10') + scale_y_continuous(trans='log10')
ggplot(hashstim,aes(x=nonstim)) + geom_histogram(binwidth=0.01) +
  scale_x_continuous(trans='log10',limits=c(1,1000))
ggplot(hashstim,aes(x=stim)) + geom_histogram(binwidth=0.01) +
  scale_x_continuous(trans='log10',limits=c(1,1000))

ggplot(hashsample,aes(x=hashsample$`860A`)) + geom_histogram(binwidth=0.05) +
  scale_x_continuous(trans='log10',limits=c(1,1000))

hashstim$stimvnonstim <- ifelse((hashstim$stim == 0 & hashstim$nonstim == 0),"negative",
                          ifelse(hashstim$stim > 0.1 & hashstim$nonstim < 1, "stim",
                          ifelse(hashstim$nonstim > 0.1 & hashstim$stim < 0.1, "nonstim", 
                            ifelse(hashstim$stim > 5, "stim",
                         ifelse(hashstim$stim > hashstim$nonstim - 10, "stim", "nonstim")))))

ggplot(hashstim) + geom_point(mapping=aes(x=nonstim,y=stim,color=hashstim$stimvnonstim), size=1) +
  scale_x_continuous(trans='log10') + scale_y_continuous(trans='log10')
View(hashstim)
#HTO - Assignhash function. If threshold should be set, adjust length(hash)
assignhash <- function(x) {
  max = max(x)
  hash <- names(x)[which(x==max)]
  if(length(hash)>1) {
    if(max==0) {
      return("negative")}
    else {return("multiplet")}
  }
  else {return(hash)}
}
hashsample$label <- apply(hashsample, MARGIN=1, FUN=assignhash)
View(hashsample$label)
table(hashsample$label)

#EXPORT DATA. Use only labels that are needed
hashlabels <- data.frame("stimvnonstim" = hashstim$stimvnonstim,
                        "sample" = hashsample$label, row.names=rownames(hash))
View(hashlabels)

write.csv(hashlabels, file = "ibd3_hash_labelled.csv")


