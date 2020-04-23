library(Seurat)
library(ggplot2)
library(MAST)
# Colon atlas
atlas <- ReadH5AD(file="/Users/akornberg/scanpy/data/Colon_cell_atlas.h5ad")
test <- atlas@meta.data
atlas.metadata <- read.csv(file="/Users/akornberg/Documents/Columbia/Han Lab/Data Analysis/2019_12_IBD/atlas/gex/Colon_immune_metadata.csv")
row.names(atlas.metadata) <- atlas.metadata$index
atlas@meta.data <- atlas.metadata

atlas <- AddMetaData(object=atlas, metadata = atlas.metadata)

# Initial Processing of Data
atlas$donor.region <- paste(atlas$donor, atlas$region)
atlas$experiment <- atlas$donor
atlas$percentMT <- atlas$percent_mito

# Subset out only T cells from Colon (Remove other cell types and mLN)
atlas$test <- ifelse((atlas$cell_type == "Activated CD4 T"), 1, 
              ifelse(atlas$cell_type == "CD8 T", 1,
              ifelse(atlas$cell_type == "Tcm", 1,
              ifelse(atlas$cell_type == "Tfh", 1,
              ifelse(atlas$cell_type == "Th1", 1,
              ifelse(atlas$cell_type == "Th17", 1,
              ifelse(atlas$cell_type == "Treg", 1, NA)))))))
colon <- subset(x=atlas, subset = test == 1)
colon <- subset(x=colon, subset = region != "mLN")

#Remove BCR metadata
colon@meta.data[,12:28] <- NULL
# Integrate Datasets
options(future.globals.maxSize = 4000 * 1024^2)
ibd.short <- merge(x=ibd2, y =ibd3)
ibd.object <- merge(x=colon, y=ibd.short) 
ibd.object <- merge(x=ibd.merged, y=colon.merged)

ibd.list <- SplitObject(ibd.object, split.by = "experiment")
ibd.list <- ibd.list[c("290b", "298c", "302c", "390c", "417c", "ibd2", "ibd3")]
for (i in 1:length(ibd.list)) {
  ibd.list[[i]] <- SCTransform(ibd.list[[i]], vars.to.regress = "percentMT", verbose = FALSE)
}

c("4861STDY7135911","4861STDY7135912","4861STDY7135913","4861STDY7135914","4861STDY7135915","4861STDY7135916","4861STDY7135917", "4861STDY7135918",
"4861STDY7208410","4861STDY7208411","4861STDY7208412","4861STDY7208413","4861STDY7208414","4861STDY7208415","4861STDY7208416","4861STDY7208417", 
"Human_colon_16S7255675","Human_colon_16S7255676","Human_colon_16S7255677","Human_colon_16S7255678", 
"Human_colon_16S7255679","Human_colon_16S7255680","Human_colon_16S7255681","Human_colon_16S7255682", 
"Human_colon_16S8000484","Human_colon_16S8001863","Human_colon_16S8001867","Human_colon_16S8001871","Pan_T7935487","Pan_T7935488","Pan_T7935489","Pan_T7935494") 


ibd.features <- SelectIntegrationFeatures(object.list = ibd.list, nfeatures = 3000)
ibd.list <- PrepSCTIntegration(object.list = ibd.list, anchor.features = ibd.features, 
                               verbose = FALSE)
ibd.anchors <- FindIntegrationAnchors(object.list = ibd.list, normalization.method = "SCT", 
                                      anchor.features = ibd.features, verbose = FALSE)
ibd.full <- IntegrateData(anchorset = ibd.anchors, normalization.method = "SCT", 
                            verbose = FALSE)

ibd.full <- RunPCA(ibd.full, verbose = FALSE)
ElbowPlot(ibd.full, ndims= 50)
ibd.full <- RunUMAP(ibd.full, dims = 1:25)
ibd.full <- FindNeighbors(ibd.full, dims = 1:25, verbose = FALSE)
ibd.full <- FindClusters(ibd.full, verbose = FALSE, resolution = )
DimPlot(ibd.full, label = TRUE, pt.size = 0.25) + ggtitle("UMAP Clusters") + theme(plot.title = element_text(hjust = 0.5)) 
DefaultAssay(colon.merged) <- "RNA"
colon.merged <- NormalizeData(colon.merged, verbose = FALSE)


colon.merged <- RunPCA(colon.merged, verbose = FALSE)
ElbowPlot(colon.merged, ndims= 50)
colon.merged <- RunUMAP(colon.merged, dims = 1:25)
colon.merged <- FindNeighbors(colon.merged, dims = 1:25, verbose = FALSE)
colon.merged <- FindClusters(colon.merged, verbose = FALSE, resolution = )
DimPlot(colon.merged, label = TRUE, pt.size = 0.25) + ggtitle("UMAP Clusters") + theme(plot.title = element_text(hjust = 0.5)) 
DefaultAssay(colon.merged) <- "RNA"
colon.merged <- NormalizeData(colon.merged, verbose = FALSE)

DimPlot(ibd.full, group.by = "patient", pt.size = 0.25, na.value = "grey90")  + 
  ggtitle("Cell Type") + theme(plot.title = element_text(hjust = 0.5)) 

DimPlot(colon.merged, group.by = "region", pt.size = 0.25, na.value = "grey90", order=c("Sigmoid")) + 
  ggtitle("Sigmoid") + theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values=c("gray90", "gray90", "gray90", "firebrick", "gray90")) 

DimPlot(ibd.full, group.by = "activity.full", pt.size = 0.25, na.value = "grey90", order=c("active", "inactive") ) + 
  ggtitle("Merged") + theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values=c("gray90", "gray90", "olivedrab", "deepskyblue2", "firebrick", "gray90", "firebrick")) 

   FeaturePlot(colon.merged, features = c("HLA-E"),ncol= 1, min.cutoff = , max.cutoff = ) 
FeaturePlot(ibd.merged, features = c("CD103"),ncol= 1, min.cutoff = , max.cutoff = ) 

cluster.atlas.5v13 <- FindMarkers(colon.merged, ident.1 = 5 ,ident.2= 13, logfc.threshold = 0.15, test.use = "MAST")

DotPlot(colon.merged, features = c("EOMES", "KLRG1", "GZMA", "CD28", "TMIGD2")) + RotatedAxis()
DotPlot(colon.merged, features = c("GZMK","GZMA","CCL5","CD74","PRF1", "EOMES", "KLRG1", "CD27", "CD28")) + RotatedAxis() + ggtitle("Inflamed Cluster Relevent Genes") + theme(plot.title = element_text(hjust = 0.5)) 
#REMOVE TCR GENES FROM ATLAS
features <- atlas@assays$RNA@counts@Dimnames[[1]]

features.tcr <- grep(c("^TRAV|^TRAJ|^TRAC|^TRBV|^TRBD|^TRBJ|^TRBC"), features, value = TRUE)

atlas.test <- atlas@assays$RNA@counts@Dimnames[[1]][!rownames(ibd.3) %in% features.3.tcr, ]

# Add Corrected TCR GENES 
tcr.390.cae <- read.csv("/Users/akornberg/Documents/Columbia/Han Lab/Data Analysis/2019_12_IBD/atlas/tcr/390c/filtered_contig_annotations_390c_498_TCR_CAE.csv")
tcr.390.tcl <- read.csv("/Users/akornberg/Documents/Columbia/Han Lab/Data Analysis/2019_12_IBD/atlas/tcr/390c/filtered_contig_annotations_390c_499_TCR_TCL.csv")
tcr.390.scl <- read.csv("/Users/akornberg/Documents/Columbia/Han Lab/Data Analysis/2019_12_IBD/atlas/tcr/390c/filtered_contig_annotations_390c_500_TCR_SCL.csv")

tcr.417.cae <- read.csv("/Users/akornberg/Documents/Columbia/Han Lab/Data Analysis/2019_12_IBD/atlas/tcr/417c/filtered_contig_annotations_417c_831_TCR_CAE.csv")
tcr.417.tcl <- read.csv("/Users/akornberg/Documents/Columbia/Han Lab/Data Analysis/2019_12_IBD/atlas/tcr/417c/filtered_contig_annotations_417c_833_TCR_TCL.csv")
tcr.417.scl <- read.csv("/Users/akornberg/Documents/Columbia/Han Lab/Data Analysis/2019_12_IBD/atlas/tcr/417c/filtered_contig_annotations_417c_835_TCR_SCL.csv")
tcr.417.mln <- read.csv("/Users/akornberg/Documents/Columbia/Han Lab/Data Analysis/2019_12_IBD/atlas/tcr/417c/filtered_contig_annotations_417c_838_TCR_MLN.csv")






atlas <- SCTransform(atlas, vars.to.regress = "percent.mito", verbose = FALSE)
atlas <- RunPCA(atlas, verbose = FALSE)
atlas <- RunUMAP(atlas, dims = 1:30, verbose = FALSE)
pbmc <- FindNeighbors(pbmc, dims = 1:30, verbose = FALSE)
pbmc <- FindClusters(pbmc, verbose = FALSE)
DimPlot(pbmc, label = TRUE) + NoLegend()

atlas.metadata <- read.csv(file="/Users/akornberg/Documents/Columbia/Han Lab/Data Analysis/2019_12_IBD/atlas/gex/Colon_immune_metadata.csv")
atlas <- AddMetaData(object=atlas, metadata = atlas.metadata)
save.image(file="atlas.RData")

# Take T cells
# Integrate w/SCTransform protocol
## Donor Info ## 
# 0 - 290b
# 1 - 298c
# 2 - 302c
# 3 - 390c
# 4 - 417c

## Region Info ## 
# 0 - Cecum
# 1 - Sigmoid Colon
# 2 - Transverse Colon
# 3 - MLN

## Cell Type Info ##
# 0 - Activated CD4
# 1 - B IgA
# 2 - B IgG
# 3 - B cycling
# 4 - B follicular
# 5 - B cell memory
# 6 - CD8
# 7 - ILC
# 8 - Lymphoid DC
# 9 - Monocyte
# 10 - Mast
# 11 - Macrophages
# 12 - LYVE1 Macrophage
# 13 - NK
# 14 - Tcm
# 15 - Tfh
# 16 - Th1
# 17 - Th17
# 18 - Treg
# 19 - CDC 1
# 20 - CDC 2
# 21 - Cycling DC
# 22 - PDC
# 23 - GD T cells
# 24 - Cycling GD T cells