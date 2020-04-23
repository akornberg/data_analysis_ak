#Run scTransform individually for each object. 
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
library(clustree)
library(cluster, quietly = TRUE)
library(monocle3)

plot1 <- FeatureScatter(ibd.merged, feature1 = "nCount_RNA", feature2 = "percentMT") + geom_hline(yintercept = 12) + geom_vline(xintercept = 8000) + geom_vline(xintercept = 400)
plot2 <- FeatureScatter(ibd.merged, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_hline(yintercept = 2500) + geom_hline(yintercept = 400) + geom_vline(xintercept = 8000) + geom_vline(xintercept = 400)
CombinePlots(plots = list(p1,p2,p3,p4))

ibd.object <- merge(x=ibd2, y =ibd3)
ibd.object <- subset(x=ibd.object, subset = activity != "unknown" & activity != "control" & dx != "cd" & patient != 730 & label != "negative" & nCount_RNA > 400 & nFeature_RNA < 2500 & percentMT < 12 & nCount_RNA<8000)

ibd.list <- SplitObject(ibd.object, split.by = "experiment")
ibd.list <- ibd.list[c("ibd2", "ibd3")]
for (i in 1:length(ibd.list)) {
  ibd.list[[i]] <- SCTransform(ibd.list[[i]], vars.to.regress = "percentMT", verbose = FALSE)
}

ibd.features <- SelectIntegrationFeatures(object.list = ibd.list, nfeatures = 3000)
ibd.list <- PrepSCTIntegration(object.list = ibd.list, anchor.features = ibd.features, 
                                    verbose = FALSE)
ibd.anchors <- FindIntegrationAnchors(object.list = ibd.list, normalization.method = "SCT", 
                                           anchor.features = ibd.features, verbose = FALSE)
ibd.merged <- IntegrateData(anchorset = ibd.anchors, normalization.method = "SCT", 
                                     verbose = FALSE)

#Anchor positions have been found, we will proceed with normal analysis from here! 
ibd.merged <- RunPCA(ibd.merged, verbose = FALSE)
ElbowPlot(ibd.merged, ndims= 50)
ibd.merged <- RunUMAP(ibd.merged, dims = 1:25)
ibd.merged <- FindNeighbors(ibd.merged, dims = 1:25, verbose = FALSE)
ibd.merged <- FindClusters(ibd.merged, verbose = FALSE, resolution = 1)
DimPlot(ibd.merged, label = TRUE, pt.size = 0.75) + ggtitle("UMAP Clusters") + theme(plot.title = element_text(hjust = 0.5)) 

ibd.merged$test <- ifelse((ibd.merged$nCount_RNA > 12000), "High", "Low")
Idents(ibd.merged) <- ibd.merged$test
#Optional clustertree
clustree(ibd.merged, prefix="integrated_snn_res.")
#We have our general plot, now we can plot out interesting metadata
DimPlot(ibd.merged, group.by = "label", pt.size = 0.75, na.value = "grey90") + 
  ggtitle("CD4/CD8") + theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values=c("darkgreen", "indianred3", "blue", "gray90")) 

DimPlot(ibd.merged, group.by = "label", pt.size = 0.75, na.value = "grey90") + 
  ggtitle("Sample A vs Sample B") + theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values=c("gray90", "darkorange2", "blue", "gray90"), label = c("Sample A", "Sample B")) 


DimPlot(ibd.merged, group.by = "stim.status", pt.size = 0.75, na.value = "grey90") + 
  ggtitle("Stim vs Nonstim") + theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values=c("gray90", "navy", "yellowgreen")) 

DimPlot(ibd.merged, group.by = "activity", pt.size =0.75 , na.value = "grey90") + 
  ggtitle("Inflammation") + theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values=c("firebrick", "steelblue2", "steelblue2", "gray90")) 

DimPlot(ibd.merged, group.by = "experiment" ,pt.size = 0.75, na.value = "grey90") + 
  ggtitle("Experiment") + theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values=c("darkgoldenrod2", "navy")) 

DimPlot(ibd.merged, group.by = "cells" ,pt.size = 0.75, na.value = "grey90") + 
  ggtitle("Treatment") + theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values=c("gray50", "indianred", "gray50", "navy", "plum", "firebrick", "darkorange")) 

DimPlot(ibd.merged, group.by = "activity.and.stim" ,pt.size = 0.75, na.value = "grey90") + 
  ggtitle("Activity & Stim") + theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values=c("gray90", "gray90", "gray90", "gray90", "navy", "deepskyblue2", "gray90")) 

DimPlot(ibd.merged, group.by = "region.general", split.by = "activity",pt.size = 0.5, na.value = "grey90") + 
  ggtitle("Region") + theme(plot.title = element_text(hjust = 0.5)) + 
  scale_color_manual(values=c("gray90", "gray90", "firebrick", "gray90", "navy", "gray90",
  "gray90", "firebrick", "gray90", "gray90", "gray90", "gray90", "gray90", "gray90"))

ibd.merged$region.general <- ifelse((ibd.merged$region == "Rectum"), "Rectum",
                              ifelse(ibd.merged$region == "Ileum (terminal)", "Ileum",
                              ifelse(ibd.merged$region == "Ileum (Nonspecific)", "Ileum",
                              ifelse(ibd.merged$region == "Sigmoid Colon (nonspecific)", "Sigmoid Colon",
                              ifelse(ibd.merged$region == "Ascending Colon (nonspecific)", "Ascending Colon",
                              ifelse(ibd.merged$region == "Cecum", "Cecum", "Other"))))))
# Create new "Treated Categories"
ibd.merged$treated <- ifelse((ibd.merged$treatment == "none"), "None", "Treated")

DimPlot(ibd.merged, group.by = "treated" ,pt.size = 0.75, na.value = "grey90") + 
  ggtitle("Treated") + theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values=c("darkgoldenrod2", "forestgreen")) 

  ## Count aminosalicylates as no treatment
ibd.merged$treatedadj <- ifelse((ibd.merged$treatment == "none"), "None",
                         ifelse(ibd.merged$treatment == "Aminosalicylate", "None", "Treated"))

DimPlot(ibd.merged, group.by = "treatedadj" ,pt.size = 0.75, na.value = "grey90") + 
  ggtitle("Treated Adj") + theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values=c("darkgoldenrod2", "forestgreen")) 
# Group by patient sample -- CAN PLOT INDIVIDUAL PATIENT
DimPlot(ibd.merged, group.by = "patient.sample", pt.size = 0.5, na.value = "grey90") + 
  ggtitle("Patient 705 - F 38 - Pancolitis - Untreated") + theme(plot.title = element_text(hjust = 0.5)) + 
  scale_color_manual(values=c("navy", "firebrick", "gray90", "gray90", "gray90", "gray90",
                              "gray90", "gray90", "gray90", "gray90", "gray90", "gray90",
                              "gray90", "gray90", "gray90", "gray90"))

DimPlot(ibd.merged, group.by = "patient", pt.size = 0.5, na.value = "grey90", order=c("867")) + 
  ggtitle("Sample B") + theme(plot.title = element_text(hjust = 0.5)) + 
  scale_color_manual(values=c("gray90", "gray90", "gray90", "gray90", "gray90", "gray90",
                              "gray90", "black", "gray90", "gray90", "gray90", "gray90",
                              "gray90", "gray90", "gray90", "gray90"))

# Group by condition
DimPlot(ibd.merged, group.by = "dx", pt.size = 0.5) + 
  ggtitle("Diagnosis") + theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values=c("navy", "peru", "navy", "gray90")) 

##GEX PLOT##
FeaturePlot(ibd.merged, features = c("TCF7"),ncol= 1, min.cutoff = 0, max.cutoff = ) 
scale_color_gradient(low="gray90",high="forestgreen", na.value = "gray90") 

##AdT PLOT##
FeaturePlot(ibd.merged, features = c("adt_CD69"), ncol= 1, min.cutoff = "q60", max.cutoff = "q95")

#compare AdT normalization strategies

ggplot(test4,aes(x=CD3)) + geom_histogram(binwidth=0.05) +
  scale_x_continuous(trans='log10',limits=c(-0.1,5)) + 
  geom_vline(xintercept = 2.1, color = "red") + geom_vline(xintercept = 0.3, color = "yellow")

#DE
# Differential Expression - MAST 
cluster.0<- FindMarkers(ibd.merged, ident.1 = 0 ,ident.2= , logfc.threshold = 0.15, test.use = "MAST")

ibd.merged$cluster53 <- ifelse((ibd.merged$seurat_clusters == 5 | ibd.merged$seurat_clusters == 3), "c5.c3", ibd.merged$seurat_clusters)
# Compare cluster ID to only other CD4 or CD8 cells
cluster.id <- select(ibd.merged@meta.data, "seurat_clusters", "label")
for(i in levels(as.factor(ibd.merged@meta.data$seurat_clusters))){
  ibd.merged@meta.data[[i]]<-ibd.merged@meta.data[i,2] 
  ibd.merged@meta.data[,ncol(ibd.merged@meta.data)] <-ifelse((ibd.merged@meta.data$seurat_clusters == i), i, ibd.merged@meta.data$label)
  
}
colnames(ibd.merged@meta.data[58:66]) <- c("cluster.0", "cluster.1", "cluster.2", "cluster.3",
                   "cluster.4", "cluster.5", "cluster.6", "cluster.7",
                   "cluster.8")


ibd.merged@meta.data[58:66] <- as.factor(ibd.merged@meta.data[58:66])

# Adjusted CDR3B frequency values
ibd.merged$sfadj <- ifelse((ibd.merged$TCR.sample < 25), NA, ibd.merged$pC.sample)

#Compare Total v Sample Freq for CD4 & CD8
ibd.merged$label <- as.factor(ibd.merged$label)
ibd.merged$sfadj.cd4 <- ifelse((ibd.merged$label == "CD4"), ibd.merged$sfadj, NA)
ibd.merged$sfadj.cd8 <- ifelse((ibd.merged$label == "CD8"), ibd.merged$sfadj, NA)


FeaturePlot(ibd.merged, features = "sfadj", pt.size = 0.75) + ggtitle("Sample TCRB Frequency") + 
  scale_color_gradient(low="lightgoldenrod2",high="firebrick", na.value = "grey90") 

FeaturePlot(ibd.merged, features = "sfadj.cd4", pt.size = 0.75) + ggtitle("Sample CD4 TCRB Frequency") + 
  scale_color_gradient(low="lightgoldenrod2",high="firebrick", na.value = "grey90") 

FeaturePlot(ibd.merged, features = "sfadj.cd8", pt.size = 0.75) + ggtitle("Sample CD8 TCRB Frequency") + 
  scale_color_gradient(low="lightgoldenrod2",high="firebrick", na.value = "grey90") 
#Compare clonality in inflamed and noninflamed tissues
ibd.merged$activity <- as.factor(ibd.merged$activity)

ibd.merged$sfadj.active <- ifelse((ibd.merged$activity == "active"), ibd.merged$sfadj, NA)
ibd.merged$sfadj.active.cd4 <- ifelse((ibd.merged$activity == "active"), ibd.merged$sfadj.cd4, NA)
ibd.merged$sfadj.active.cd8 <- ifelse((ibd.merged$activity == "active"), ibd.merged$sfadj.cd8, NA)

FeaturePlot(ibd.merged, features = "sfadj.active") + ggtitle("Active Inflammation") + 
  scale_color_gradient(low="lightgoldenrod2",high="firebrick", na.value = "grey90") 

ibd.merged$sfadj.inactive <- ifelse((ibd.merged$activity == "inactive"), ibd.merged$sfadj, NA)
ibd.merged$sfadj.inactive.cd4 <- ifelse((ibd.merged$activity == "inactive"), ibd.merged$sfadj.cd4, NA)
ibd.merged$sfadj.inactive.cd8 <- ifelse((ibd.merged$activity == "inactive"), ibd.merged$sfadj.cd8, NA)

FeaturePlot(ibd.merged, features = "sfadj.cd4") + ggtitle("Inactive Inflammation") + 
  scale_color_gradient(low="gray90",high="firebrick", na.value = "grey90") 

# Calculate Frequency of TCR per patient sample
tcr.general <- select(ibd.merged@meta.data, patient.sample, cdr3b, pC.sample)
tcr.general$barcode <- row.names(tcr.general)

cdr3b <- levels(as.factor(tcr.general$cdr3b))

for(i in levels(as.factor(ibd.merged$patient.sample))){
  tcr.general[[i]]<-tcr.general[i,2] 
  tcr.general[,ncol(tcr.general)] <- ifelse((tcr.general$patient.sample == i), tcr.general$pC.sample, NA)
}


tcr.sample <- select(tcr.general, cdr3b,"705A","705B","756B","794A","794B","811A","811B","816A","816B","860A","860B",         
"867A","867B","867C","878A","878B")

cdr3.sample <- as.data.frame(setDT(tcr.sample)[, lapply(.SD, function(x) sort(x)[1L]), by = .(cdr3b)])

to.merge <- select(tcr.general, cdr3b, barcode)

tcr.merged <- merge(cdr3.sample, to.merge, by = 'cdr3b', all = TRUE)

rownames(tcr.merged) <- tcr.merged$barcode
tcr.merged$barcode <- NULL
colnames(tcr.merged) <- c("cdr3", "patient_705A","patient_705B","patient_756B","patient_794A","patient_794B","patient_811A",
                          "patient_811B","patient_816A","patient_816B","patient_860A","patient_860B","patient_867A",
                          "patient_867B","patient_867C","patient_878A","patient_878B")

tcr.merged[,2:ncol(tcr.merged)][is.na(tcr.merged[,2:ncol(tcr.merged)])] <- 0

ibd.merged <- AddMetaData(object=ibd.merged, metadata = tcr.merged)
##Shared CDR3B's
ibd.merged$shared <- ifelse((ibd.merged$patient_705A > 0 & ibd.merged$patient_705B > 0), "Shared",
                     ifelse(ibd.merged$patient_794A > 0 & ibd.merged$patient_794B > 0, "Shared",
                     ifelse(ibd.merged$patient_811A > 0 & ibd.merged$patient_811B > 0, "Shared",
                     ifelse(ibd.merged$patient_816A > 0 & ibd.merged$patient_816B > 0, "Shared",
                     ifelse(ibd.merged$patient_860A > 0 & ibd.merged$patient_860B > 0, "Shared",
                     ifelse(ibd.merged$patient_867A > 0 & ibd.merged$patient_867B > 0, "Shared",
                     ifelse(ibd.merged$patient_867A > 0 & ibd.merged$patient_867C > 0, "Shared",
                     ifelse(ibd.merged$patient_878A > 0 & ibd.merged$patient_878B > 0, "Shared", "Other"))))))))
        
DimPlot(ibd.merged, group.by = "shared", split.by = "patient",pt.size = 0.5, order = c("Shared")) + 
  ggtitle("Diagnosis") + theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values=c("gray90", "navy", "gray90", "gray90"))               
                  
ibd.merged$shared <- as.factor(ibd.merged$shared)
ibd.merged$shared.activity <- ifelse((ibd.merged$activity == "active" & ibd.merged$shared == "Shared"), "Shared CDR3B Active",
                              ifelse(ibd.merged$activity == "inactive" & ibd.merged$shared == "Shared", "Shared CDR3B Inactive", "Other"))

DimPlot(ibd.merged, group.by = "shared.activity", pt.size = 0.5, order = c("Shared CDR3B Inactive", "Shared CDR3B Active")) + 
  ggtitle("Shared TCR's") + theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values=c("gray90", "firebrick","navy"))

ibd.merged$shared.activity <- as.factor(ibd.merged$shared.activity)
ibd.merged$shared.activity.patient <- ifelse((ibd.merged$shared.activity == "Shared CDR3B Active"), paste(ibd.merged$patient, ibd.merged$shared.activity), 
                                      ifelse(ibd.merged$shared.activity == "Shared CDR3B Inactive",paste(ibd.merged$patient, ibd.merged$shared.activity), "Other"))

DimPlot(ibd.merged, group.by = "shared.activity.patient", pt.size = 0.5, order = c("867 Shared CDR3B Inactive","867 Shared CDR3B Active")) + 
  ggtitle("Shared TCR's - Patient 867") + theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values=c("gray90", "gray90", "gray90", "gray90", "gray90", "gray90", "gray90", "firebrick", "navy"))

ibd.merged$cdr3 <- as.character(ibd.merged$cdr3)
ibd.merged$test <- ifelse((ibd.merged$shared.activity.patient == "867 Shared CDR3B Active"), ibd.merged$cdr3,
                   ifelse(ibd.merged$shared.activity.patient == "867 Shared CDR3B Inactive", ibd.merged$cdr3, "Other"))

#Expansion // Contraction
##ADJUSTABLE THRESHOLDS
expand.threshold <- 0.03
contract.threshold <- 0.03 

ibd.merged$expand.active <- ifelse((ibd.merged$patient_705B > ibd.merged$patient_705A + expand.threshold), "Inflamed Expansion",
                            ifelse(ibd.merged$patient_794B > ibd.merged$patient_794A + expand.threshold, "Inflamed Expansion",
                            ifelse(ibd.merged$patient_811A > ibd.merged$patient_811B + expand.threshold, "Inflamed Expansion", 
                            ifelse(ibd.merged$patient_816A > ibd.merged$patient_816B + expand.threshold, "Inflamed Expansion", 
                            ifelse(ibd.merged$patient_811A > ibd.merged$patient_811B + expand.threshold, "Inflamed Expansion", 
                            ifelse(ibd.merged$patient_816A > ibd.merged$patient_816B + expand.threshold, "Inflamed Expansion", 
                            ifelse(ibd.merged$patient_860B > ibd.merged$patient_860A + expand.threshold, "Inflamed Expansion", 
                            ifelse(ibd.merged$patient_867A > ibd.merged$patient_867B + expand.threshold, "Inflamed Expansion",
                            ifelse(ibd.merged$patient_867A > ibd.merged$patient_867C + expand.threshold, "Inflamed Expansion", 
                            ifelse(ibd.merged$patient_876A > ibd.merged$patient_876B + expand.threshold, "Inflamed Expansion", 
                            ifelse(ibd.merged$patient_878A > ibd.merged$patient_878B + expand.threshold, "Inflamed Expansion", "Other")))))))))))
                                                 
                                                        
ibd.merged$inflamed.contraction <-  ifelse((ibd.merged$patient_705B < ibd.merged$patient_705A - contract.threshold), "Inflamed Contraction",
                               ifelse(ibd.merged$patient_794B < ibd.merged$patient_794A - contract.threshold, "Inflamed Contraction",
                               ifelse(ibd.merged$patient_811A < ibd.merged$patient_811B - contract.threshold, "Inflamed Contraction", 
                               ifelse(ibd.merged$patient_816A < ibd.merged$patient_816B - contract.threshold, "Inflamed Contraction", 
                               ifelse(ibd.merged$patient_811A < ibd.merged$patient_811B - contract.threshold, "Inflamed Contraction", 
                               ifelse(ibd.merged$patient_816A < ibd.merged$patient_816B - contract.threshold, "Inflamed Contraction", 
                               ifelse(ibd.merged$patient_860B < ibd.merged$patient_860A - contract.threshold, "Inflamed Contraction", 
                               ifelse(ibd.merged$patient_867A < ibd.merged$patient_867B - contract.threshold, "Inflamed Contraction",
                               ifelse(ibd.merged$patient_867A < ibd.merged$patient_867C - contract.threshold, "Inflamed Contraction", 
                               ifelse(ibd.merged$patient_876A < ibd.merged$patient_876B - contract.threshold, "Inflamed Contraction", 
                               ifelse(ibd.merged$patient_878A < ibd.merged$patient_878B - contract.threshold, "Inflamed Contraction", "Other")))))))))))

ibd.merged$activity <- as.factor(ibd.merged$activity)
ibd.merged$expand.active <- as.factor(ibd.merged$expand.active)
ibd.merged$inflamed.contraction <- as.factor(ibd.merged$inflamed.contraction)


ibd.merged$active.expansion <- ifelse((ibd.merged$activity == "active" & ibd.merged$expand.active == "Inflamed Expansion"), "Inflamed Expansion",
                         ifelse(ibd.merged$activity == "inactive" & ibd.merged$expand.active == "Inflamed Expansion", "Noninflamed Expansion", "Other"))

ibd.merged$active.contraction <- ifelse((ibd.merged$activity == "active" & ibd.merged$inflamed.contraction == "Inflamed Contraction"), "Inflamed Contraction",
                                      ifelse(ibd.merged$activity == "inactive" & ibd.merged$inflamed.contraction == "Inflamed Contraction", "Noninflamed Contraction", "Other"))

ibd.merged$active.expansion <- as.factor(ibd.merged$active.expansion)
ibd.merged$active.contraction <- as.factor(ibd.merged$active.contraction)

ibd.merged$active.expansionadj <- ifelse((ibd.merged$TCR.sample < 30), "Other", 
                                  ifelse(ibd.merged$TCR.sample > 30, ibd.merged$active.expansion, "Other"))
ibd.merged$active.contractionadj <- ifelse((ibd.merged$TCR.sample < 30), "Other",
                                  ifelse(ibd.merged$TCR.sample > 30, ibd.merged$active.contraction, "Other"))

DimPlot(ibd.merged, group.by = "active.expansion",split.by = "patient", pt.size = 0.5, order=c("Noninflamed Expansion", "Inflamed Expansion")) + 
  ggtitle("Active Inflammation TCR > Inactive") + theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values=c("gray90", "navy", "deepskyblue2", "gray90")) 

DimPlot(ibd.merged, group.by = "active.contractionadj",split.by = "patient", pt.size = 0.5, order=c("Noninflamed Contraction", "Inflamed Contraction")) + 
  ggtitle("Inactive Inflammation TCR > Active") + theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values=c("gray90", "navy", "deepskyblue2", "gray90")) 

## Calculate TCR per cluster
 

# TCR Expansion
ibd.merged$expansions <- ibd.merged$freq.sample - 1
ibd.merged$expansionfreq <- ibd.merged$expansions / ibd.merged$TCR.sample
ibd.merged$expansionfreq[ibd.merged$expansionfreq == 0] <- NA

FeaturePlot(ibd.merged, features = "CD4", split.by = "label")+ ggtitle("Expansion Freq") + 
  scale_color_gradient(low="gray90",high="firebrick", na.value = "grey90") 

ibd.merged$expansionfreq.cd4 <- ifelse((ibd.merged$label == "CD4"), ibd.merged$expansionfreq, NA)
ibd.merged$expansionfreq.cd8 <- ifelse((ibd.merged$label == "CD8"), ibd.merged$expansionfreq, NA)

FeaturePlot(ibd.merged, features = "expansionfreq.cd4") + ggtitle("Expansion Freq CD4") + 
  scale_color_gradient(low="lightgoldenrod2",high="firebrick", na.value = "grey90") 

FeaturePlot(ibd.merged, features = "expansionfreq.cd8") + ggtitle("CD8 Expanded T Cells") + 
  scale_color_gradient(low="lightgoldenrod2",high="firebrick", na.value = "grey90") 


### FIGURES FOR QUALIFYING EXAM
DimPlot(ibd.merged, label = TRUE, pt.size = 0.75) + ggtitle("UMAP Clusters") + theme(plot.title = element_text(hjust = 0.5),legend.text = element_text(size=20, face ="bold")) + guides(colour = guide_legend(override.aes = list(size=7.5)))

DimPlot(ibd.merged, label = FALSE, pt.size = 0.75, group.by = "cluster.label") + ggtitle("UMAP Clusters") + theme(plot.title = element_text(hjust = 0.5),legend.text = element_text(size=10, face =)) + guides(colour = guide_legend(override.aes = list(size=7.5)))

DimPlot(ibd.merged, group.by = "activity", pt.size =0.75 , na.value = "grey90") + 
  ggtitle("Inflammation") + theme(plot.title = element_text(hjust = 0.5), legend.text = element_text(size=20, face ="bold")) + 
  scale_color_manual(values=c("firebrick", "steelblue2", "steelblue2", "gray90"), labels= c("Active", "Inactive")) + guides(colour = guide_legend(override.aes = list(size=7.5)))

DimPlot(ibd.merged, group.by = "label", pt.size = 0.75, na.value = "grey90") + 
  ggtitle("CD4 vs CD8") + theme(plot.title = element_text(hjust = 0.5), legend.text = element_text(size=20, face = "bold")) +
  scale_color_manual(values=c("darkorange3","darkgreen")) +  guides(colour = guide_legend(override.aes = list(size=7.5)))

CombinePlots(plots = list(cd28,klrg1,tcr ), ncol = 3)

EnhancedVolcano(cluster.5.3, lab=rownames(cluster.5.3), x = "avg_logFC", y = "p_val", FCcutoff = 1,pCutoff = 10e-7,
                title = "Differential Gene Expression", subtitle= "Non-inflamed (Left) vs. Inflamed (Right)", ylim = c(0,40),
                col = c("gray50", "gray50", "gray50", "firebrick"))
###Cluster ID List###
# 0 - TIGIT high - Treg-like?
# 1 - Highly activated, low cytotoxicity
# 2 - IL7R high 
# 3 - Tcm (CCR7+ CD62L+)
# 4 - Tissue resident / Regulatory CD8's
# 5 - Low proliferation CD4
# 6 - Cytotoxic CD8+ T cells
# 7 - FOXP3 Tregs (IL32/Fas/Lag3/CD74)
# 8 - Activated CD4 Th1/Th17
# 9 - TNF producing, activated Th17
# 10 - Pro-inflammatory Th17 

ibd.merged$cluster.label <- ifelse((ibd.merged$seurat_clusters == 0), "0 - TIGIT High",
                            ifelse(ibd.merged$seurat_clusters == 1, "1 - Stimulated CD4/CD8",
                            ifelse(ibd.merged$seurat_clusters == 2, "2 - IL7R+", 
                            ifelse(ibd.merged$seurat_clusters == 3, "3 - Tcm",
                            ifelse(ibd.merged$seurat_clusters == 4, "4 - Regulatory CD8", 
                            ifelse(ibd.merged$seurat_clusters == 5, "5 - Resting CD4",
                            ifelse(ibd.merged$seurat_clusters == 6, "6 - Cytotoxic CD8", 
                            ifelse(ibd.merged$seurat_clusters == 7, "7 - Tregs",
                            ifelse(ibd.merged$seurat_clusters == 8, "8 - Stimulated CD4/CD8", 
                            ifelse(ibd.merged$seurat_clusters == 9, "9 - Stimulated Th1/Th17",
                            ifelse(ibd.merged$seurat_clusters == 10, "l0 - IL17 Producing", NA)))))))))))

DimPlot(ibd.merged, label = FALSE, pt.size = 0.75, group.by = "cluster.label") + ggtitle("UMAP Clusters") + theme(plot.title = element_text(hjust = 0.5),legend.text = element_text(size=15, face ="bold")) + guides(colour = guide_legend(override.aes = list(size=7.5)))

FeaturePlot(ibd.merged, features = c("IL17F"),ncol= 1, min.cutoff = 0, max.cutoff = 1) 
scale_color_gradient(low="gray90",high="forestgreen", na.value = "gray90") 

FeaturePlot(ibd.merged, features = c("adt_ICOS"), ncol= 1, min.cutoff = "q20", max.cutoff = "q95")

ibd.merged$stim.c6.c4 <- ifelse((ibd.merged$stim.status == "nonstim"), "nonstim",
                         ifelse(ibd.merged$seurat_clusters == 4, 4,
                         ifelse(ibd.merged$seurat_clusters == 6, 6, "other")))

ibd.merged$nonstim.c6.c4 <- ifelse((ibd.merged$stim.status == "stim"), "stim",
                                ifelse(ibd.merged$seurat_clusters == 4, 4,
                                       ifelse(ibd.merged$seurat_clusters == 6, 6, "other")))

cluster.stim.6.4 <- FindMarkers(ibd.merged, ident.1 = 6 ,ident.2= 4, logfc.threshold = 0.15, test.use = "MAST")

## OTHER ##
ibd.merged$c2expansion <- ifelse((ibd.merged$expansionfreq > 0 & ibd.merged$seurat_clusters == 2), "C2_Expand", 
                          ifelse(ibd.merged$seurat_clusters == 5, "C5", "Other"))
cluster.stim.split<- FindMarkers(ibd.merged, ident.1 = 1 ,ident.2= 9 , logfc.threshold = 0.05,min.diff.pct = ,  test.use = "MAST")
# CD4 CD8 GEX vs AdT
test <- as.data.frame(ibd.merged@assays$SCT@scale.data)
test2 <- as.data.frame(t(as.data.frame(test)))
test3 <- select(test2, CD4, CD8B)
ibd.merged <- AddMetaData(object=ibd.merged, metadata = test3)

ggplot(ibd.merged@meta.data) + geom_point(mapping=aes(x=CD4,y=CD8B, color = ibd.merged$label), size=1)

ggplot(test3,aes(x=CD4)) + geom_histogram(binwidth=0.01)

ggplot(test3,aes(x=CD8B)) + geom_histogram(binwidth=0.01) 

EnhancedVolcano(cluster.nonstim.6.4, lab=rownames(cluster.nonstim.6.4), x = "avg_logFC", y = "p_val", FCcutoff = 1,pCutoff = 10e-7,
                title = "Differential Gene Expression - Unstimulated Cells", subtitle= "Non-inflamed (Left) vs. Inflamed (Right)", ylim = c(0,40),
                col = c("gray50", "gray50", "gray50", "firebrick"), drawConnectors = TRUE,widthConnectors = 0.5,    labCol = 'black',
                labFace = 'bold',
                boxedLabels = TRUE,
                legendPosition = "top",
                labSize = 7.5,
                selectLab = c("EOMES", "GZMK", "TMIGD2", "NKG7", "IL7R", "CCR9", "KLRG1", "CXCR3"))
if (!requireNamespace("BiocManager", quietly = TRUE)) +
  install.packages("BiocManager")
BiocManager::install()

BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
                       'limma', 'S4Vectors', 'SingleCellExperiment',
                       'SummarizedExperiment', 'batchelor'))
BiocManager::install(c('monocle3'))

install.packages("devtools")
devtools::install_github('cole-trapnell-lab/leidenbase')
devtools::install_github('cole-trapnell-lab/monocle3')

expression_matrix <- as.matrix(GetAssayData(ibd3,assay = "RNA", slot = "counts"))

expression_matrix <- eight@assays$RNA@counts
cells <- eight@meta.data
genes <- eight@assays$RNA@counts@Dimnames[[1]]
gene_df <- data.frame(ids = genes, row.names = genes)


cds <- new_cell_data_set(expression_matrix,
                         cell_metadata = cells,
                         gene_metadata = gene_df)

cds.8 <- as.CellDataSet(eight, assay = "SCT", reduction = "umap")

eight <- subset(x=ibd.merged, subset = monocle == "8114" | monocle == "8116")
  
cds <- preprocess_cds(cds, num_dim = 100)
cds <- align_cds(cds, alignment_group = "batch")
cds <- reduce_dimension(cds)
cds <- cluster_cells(cds)
cds <- learn_graph(cds)
cds <- order_cells(cds)
plot_cells(cds, genes = c("KLRG1"))

##END
