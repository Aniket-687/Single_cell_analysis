# script to perform standard workflow steps to analyze single cell RNA-Seq data
# data: 20k Mixture of NSCLC DTCs from 7 donors, 3' v3.1
# data source: https://www.10xgenomics.com/resources/datasets/10-k-human-pbm-cs-multiome-v-1-0-chromium-controller-1-standard-2-0-0     

setwd("C:/Users/dasan/Downloads/Bioinformatics_Project/Single_Cell_Analysis/scripts")

# load libraries
library(Seurat)
library(tidyverse)
library(hdf5r)
library(ggplot2)

# Load the NSCLC dataset
nsclc.sparse.m <- Read10X_h5(
  filename = "C:/Users/dasan/Downloads/Bioinformatics_Project/Single_Cell_Analysis/Data/20k_NSCLC_DTC_3p_nextgem_donor_1_count_sample_feature_bc_matrix.h5")

str(nsclc.sparse.m)
cts <-  nsclc.sparse.m$`Gene Expression`

# Initialize the Seurat object with the raw (non-normalized data).
nsclc.seurat.obj <- CreateSeuratObject(counts = cts, project = "NSCLC", min.cells = 3, min.features = 200)
str(nsclc.seurat.obj)
nsclc.seurat.obj
# 21478 features across 1575 samples

# 1. QC -------
View(nsclc.seurat.obj@meta.data)
# % MT reads
nsclc.seurat.obj[["percent.mt"]] <- PercentageFeatureSet(nsclc.seurat.obj, pattern = "^MT-")
View(nsclc.seurat.obj@meta.data)

#Output folder for plots
out_dir <- "C:/Users/dasan/Downloads/Bioinformatics_Project/Single_Cell_Analysis/Seurat_Plots/"
dir.create(out_dir, showWarnings = FALSE)  # create folder if not exists

p1 <- VlnPlot(nsclc.seurat.obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
p2 <- FeatureScatter(nsclc.seurat.obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth(method = 'lm')

ggsave(filename = "Violin_Plot.png", plot = p1, path = out_dir, width = 6, height = 5, dpi = 300)
ggsave(filename = "Scatter_Plot.png", plot = p2, path = out_dir, width = 6, height = 5, dpi = 300)


# 2. Filtering -----------------
nsclc.seurat.obj <- subset(nsclc.seurat.obj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & 
                             percent.mt < 5)

# 3. Normalize data ----------
#nsclc.seurat.obj <- NormalizeData(nsclc.seurat.obj, normalization.method = "LogNormalize", scale.factor = 10000)
# OR
nsclc.seurat.obj <- NormalizeData(nsclc.seurat.obj)
str(nsclc.seurat.obj)

# 4. Identify highly variable features --------------
nsclc.seurat.obj <- FindVariableFeatures(nsclc.seurat.obj, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(nsclc.seurat.obj), 10)

# plot variable features with and without labels
p3 <- VariableFeaturePlot(nsclc.seurat.obj)
p4 <- LabelPoints(plot = p3, points = top10, repel = TRUE)
combined_plot <- p3 + p4
ggsave(filename = "Variable_features_with_&_without_labels.png", plot = combined_plot, path = out_dir, width = 12, height = 6, dpi = 300)

# 5. Scaling -------------
all.genes <- rownames(nsclc.seurat.obj)
nsclc.seurat.obj <- ScaleData(nsclc.seurat.obj, features = all.genes)

str(nsclc.seurat.obj)

# 6. Perform Linear dimensionality reduction --------------
nsclc.seurat.obj <- RunPCA(nsclc.seurat.obj, features = VariableFeatures(object = nsclc.seurat.obj))

# visualize PCA results
print(nsclc.seurat.obj[["pca"]], dims = 1:5, nfeatures = 5)

png(filename = paste0(out_dir, "Dimentional_Heatmap.png"), width = 2000, height = 1600, res = 300)
DimHeatmap(nsclc.seurat.obj, dims = 1, cells = 500, balanced = TRUE)
dev.off()

# determine dimensionality of the dat
p6 <- ElbowPlot(nsclc.seurat.obj)
ggsave(filename = "Elbow_Plot.png", plot = p6, path = out_dir, width = 6, height = 5, dpi = 300)

# 7. Clustering ------------
nsclc.seurat.obj <- FindNeighbors(nsclc.seurat.obj, dims = 1:15)

# understanding resolution
nsclc.seurat.obj <- FindClusters(nsclc.seurat.obj, resolution = c(0.1,0.3, 0.5, 0.7, 1))
View(nsclc.seurat.obj@meta.data)

p7 <- DimPlot(nsclc.seurat.obj, group.by = "RNA_snn_res.0.5", label = TRUE)
ggsave(filename = "RNA_snn_res_0.5.png", plot = p7, path = out_dir, width = 6, height = 5, dpi = 300)

# setting identity of clusters
Idents(nsclc.seurat.obj)
Idents(nsclc.seurat.obj) <- "RNA_snn_res.0.1"
Idents(nsclc.seurat.obj)

# non-linear dimensionality reduction --------------
# If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
# 'umap-learn')
nsclc.seurat.obj <- RunUMAP(nsclc.seurat.obj, dims = 1:15)
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
p8 <- DimPlot(nsclc.seurat.obj, reduction = "umap")
ggsave(filename = "UMAP_clustering.png", plot = p8, path = out_dir, width = 6, height = 5, dpi = 300)
