library(Seurat)
library(Matrix)

hto <- Read10X("../../Volumes/hanlab_1/Han_Lab/Data/Raw_Data/sr_exp6/umi_count/", gene.column = 1)
matrix_dir = "../../Volumes/hanlab_1/Han_Lab/Data/Raw_Data/sr_exp6/umi_count/"
barcode.path <- paste0(matrix_dir, "barcodes.tsv.gz")
features.path <- paste0(matrix_dir, "features.tsv.gz")
matrix.path <- paste0(matrix_dir, "matrix.mtx.gz")
mat <- readMM(file = matrix.path)
feature.names = read.delim(features.path, header = FALSE, stringsAsFactors = FALSE)
barcode.names = read.delim(barcode.path, header = FALSE, stringsAsFactors = FALSE)
colnames(mat) = barcode.names$V1
rownames(mat) = feature.names$V1

write.csv(mat, "exp6_adt.csv")
