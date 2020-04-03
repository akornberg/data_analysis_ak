library(circlize)
library(dplyr)
library(reshape2)

## Load in the data 
vdj <- read.csv(file.choose())

## Remove any columns without relevant data
vdjordered <- vdj[order(vdj$trb_v_genes),]
vdjcleaning <- vdjordered[!(is.na(vdjordered$trb_v_genes) | vdjordered$trb_v_genes==""), ]

View(vdjfinal)
## Remove any columns with multiple V or J regions
vdjcleaned <- vdjcleaning[!grepl(';'), vdjcleaning$trb_v_genes]

vdjcleaned <- vdjcleaning[!(is.na(vdjcleaning$trb_v_genes) | vdjcleaning$trb_v_genes==";"), ]
View(vdjcleaning)






View(vdj)
chordDiagram(vdj) 
title("KM Circular Plot - No Threshold")



Vregion <- unique(vdj$from)
Jregion <- unique(vdj$to)



vdjordered <- vdj[order(vdj$trb_v_genes),]
vdjfinal <- vdjordered[!(is.na(vdjordered$trb_v_genes) | vdjordered$trb_v_genes==""), ]
View(vdj)



View(Jregion)

Jtranspose <- t(Jregion)
View(J)


df = data.frame(from = vdj$from, to = vdj$to, value = vdj$value)
df
chordDiagram(vdj)

mat2 = matrix(vdjfinal$trb_j_genes_genes, 4512)
colnames(mat) = vdj$trb_v_genes
mat

View(mat)

cbind(mat,mat2)

chordDiagram(vdj)
circos.clear()

unique(vdj$from)


vdjmatrix <- matrix(x=vdj, nrow=nrow(vdj), ncol=ncol(vdj))
View(vdjmatrix)
