library(Seurat)

# Colon atlas
atlas <- ReadH5AD(file="/Users/akornberg/scanpy/data/Colon_cell_atlas.h5ad")
test <- atlas@meta.data
atlas.metadata <- read.csv(file="/Users/akornberg/Documents/Columbia/Han Lab/Data Analysis/2019_12_IBD/atlas/gex/Colon_immune_metadata.csv")
row.names(atlas.metadata) <- atlas.metadata$index
atlas@meta.data <- atlas.metadata

atlas <- AddMetaData(object=atlas, metadata = atlas.metadata)

# Subset out only T cells
atlas$test <- ifelse((atlas$cell_type == "Activated CD4 T"), 1, 
              ifelse(atlas$cell_type == "CD8 T", 1,
              ifelse(atlas$cell_type == "Tcm", 1,
              ifelse(atlas$cell_type == "Tfh", 1,
              ifelse(atlas$cell_type == "Th1", 1,
              ifelse(atlas$cell_type == "Th17", 1,
              ifelse(atlas$cell_type == "Treg", 1, NA)))))))
colon <- subset(x=atlas, subset = test == 1)

#Remove BCR metadata
colon@meta.data[,12:28] <- NULL

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