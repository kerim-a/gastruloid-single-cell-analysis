#load required packages
library(ggplot2)
library(cowplot)
library(dplyr)
library(Seurat)
library(Matrix)
library(viridis)
library(limma)
library(clusterProfiler)
library(org.Mm.eg.db)
library(scCustomize)
library(dittoSeq)
library(SeuratWrappers)

#load and filter 0h mESC dataset from this study to perform integration (via fastMNN) with a previously published mESC dataset of E14 mESCs grown in DMEM+serum+LIF
matrix_0h.data <- Read10X(data.dir = 'C:\\Users\\anlas\\Documents\\Phd\\RNAseq_analysis\\2020\\0h\\filtered_feature_bc_matrix', gene.column = 2, unique.features = TRUE)
g_0h <- CreateSeuratObject(counts = matrix_0h.data, project = "scRNA_gast0h_2020", min.cells = 3, min.features = 200)
g_0h[["percent.mt"]] <- PercentageFeatureSet(g_0h, pattern = "^mt-")
g0h_filtered <- subset(g_0h, subset = nFeature_RNA > 2500 & percent.mt < 10 & nCount_RNA < 150000)

#load, assess and filter alternative 0h mESC dataset (Alda-Catalinas et al., 2020)
matrix_E14_data <- Read10X(data.dir = 'Y:\\Kerim_Anlas\\Gastruloids_SingleCellAnalysis\\internet_datasets\\barbraham_data\\E14_mESCs')
E14_SL <- CreateSeuratObject(counts = matrix_E14_data, project = "ESL0h_alt", min.cells = 3, min.features = 200)
E14_SL[["percent.mt"]] <- PercentageFeatureSet(E14_SL, pattern = "^mt-")
VlnPlot(E14_SL, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
E14_SL_filtered <- subset(E14_SL, subset = nFeature_RNA > 3000 & percent.mt < 10 & nCount_RNA > 15000)
VlnPlot(E14_SL_filtered, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#define an ident to easily distinguish the 2 datasets
g0h_filtered$merge.ident <- "BraGFP"
E14_SL_filtered$merge.ident <- "E14_comp"

#merge them for integration, remove any "scale.data" information, normalize together and find variable features/genes
all_objects <- merge(g0h_filtered, y = c(E14_SL_filtered), add.cell.ids = c('BraGFP_0h','E14_0h'), project = "mESCcomp")
all_objects@assays$RNA@scale.data <- as.matrix(0)
all_objects <- NormalizeData(all_objects)
all_objects <- FindVariableFeatures(all_objects, selection.method = "vst", nfeatures = 3000)

#use these variable features for integration, however filter out known cell cycle associated genes (based on seurats built-in lists, converted to mouse)
features <- SelectIntegrationFeatures(object.list = SplitObject(all_objects, split.by = "orig.ident"), nfeatures = 3000)
features_filtered1 <- setdiff(features, g2m_genes)
features_filtered2 <- setdiff(features_filtered1, s_genes)
#run the fastMNN integration
all_objects <- RunFastMNN(object.list = SplitObject(all_objects, split.by = "orig.ident"), features = features_filtered2)

#perform clustering analysis on the integrated dataset, save etc.
all_objects <- RunUMAP(all_objects, reduction = "mnn", dims = 1:25)
all_objects <- FindNeighbors(all_objects, reduction = "mnn", dims = 1:25)
all_objects <- FindClusters(all_objects, resolution = 0.3)

p <- DimPlot(all_objects, reduction = "umap")
AugmentPlot(plot = p)
plot(p)
p <- DimPlot(all_objects, reduction = "umap", group.by = "merge.ident")
AugmentPlot(plot = p)
plot(p)

saveRDS(all_objects, file = "Y:\\Kerim_Anlas\\Gastruloids_SingleCellAnalysis\\seurat_objects\\brabrahamE14_vsTGFP_best.rds")
all_objects <- readRDS(file = "Y:\\Kerim_Anlas\\Gastruloids_SingleCellAnalysis\\seurat_objects\\brabrahamE14_vsTGFP_best.rds")

#explore the dataset (plot selected genes, find variable features, generate a heatmap of most variable genes per cluster etc)
DefaultAssay(all_objects) <- "RNA"

p <- FeaturePlot(all_objects, features = c("Zfp42","Sox2","Klf2","Pou3f1","Pim2","Nanog","T","Eomes","Wnt3","Dnmt3b","Krt8","Krt18"), order = TRUE, cols = c("lightgrey", "dodgerblue4"), pt.size = 0.5)
AugmentPlot(plot = p)
plot(p)

plasma <- viridis(30, direction = 1, option = "D") #use B or D
p <- FeaturePlot(all_objects, features = c("Krt18"), order = TRUE, cols = plasma, pt.size = 0.5) 
AugmentPlot(plot = p)
plot(p)

all_objects.markers <- FindAllMarkers(all_objects, only.pos = TRUE, min.pct = 0.5, logfc.threshold = 0.5)
top10 <- all_objects.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

DefaultAssay(all_objects) <- "mnn.reconstructed"
all_objects <- ScaleData(all_objects)
DoHeatmap(all_objects, features = top5$gene, slot = "scale.data") + scale_fill_viridis()

ggplot(all_objects@meta.data, aes(x=seurat_clusters, fill=merge.ident)) + geom_bar(position = "fill")
