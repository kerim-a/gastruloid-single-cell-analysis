#load required packages
library(ggplot2)
library(cowplot)
library(dplyr)
library(Seurat)
library(Matrix)
library(viridis)
library(limma)

#load 72h dataset
matrix_72h.data <- Read10X(data.dir = 'C:\\Users\\anlas\\Documents\\Phd\\RNAseq_analysis\\72h\\filtered_feature_bc_matrix', gene.column = 2, unique.features = TRUE)
g_72h <- CreateSeuratObject(counts = matrix_72h.data, project = "scRNA_gast72h", min.cells = 3, min.features = 200)
g_72h[["percent.mt"]] <- PercentageFeatureSet(g_72h, pattern = "^mt-")
#perform QC, normalization, and identification of most variable genes
g72h_filtered <- subset(g_72h, subset = nFeature_RNA > 2500 & percent.mt < 20 & percent.mt > 1 & nCount_RNA < 135000)
g72h_filtered <- NormalizeData(g72h_filtered, normalization.method = "LogNormalize", scale.factor = 10000)
g72h_filtered <- FindVariableFeatures(g72h_filtered, selection.method = "vst", nfeatures = 3000)
#add ID
g72h_filtered$merge.ident <- "72h"
#save
saveRDS(g72h_filtered, file = "C:\\Users\\anlas\\Documents\\Phd\\RNAseq_analysis\\outs\\g72h_filtered_forint.rds")

############
#load the 24h dataset (2nd replicate)
matrix_24h_2.data <- Read10X(data.dir = 'C:\\Users\\anlas\\Documents\\Phd\\RNAseq_analysis\\2020\\for_merging\\24h\\filtered_feature_bc_matrix', gene.column = 2, unique.features = TRUE)
g_24h_2 <- CreateSeuratObject(counts = matrix_24h_2.data, project = "scRNA_gast24h_2020", min.cells = 3, min.features = 200)
g_24h_2[["percent.mt"]] <- PercentageFeatureSet(g_24h_2, pattern = "^mt-")
#add information about replicate
g_24h_2$replicate <- "rep2"
#perform QC, normalization, and identification of most variable genes
g_24h_2 <- subset(g_24h_2, subset = nFeature_RNA > 3750 & percent.mt < 15 & percent.mt > 1 & nCount_RNA < 150000)
VlnPlot(g_24h_2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
g_24h_2 <- NormalizeData(g_24h_2, normalization.method = "LogNormalize", scale.factor = 10000)
g_24h_2 <- FindVariableFeatures(g_24h_2, selection.method = "vst", nfeatures = 3000)
#save
saveRDS(g_24h_2, file = "C:\\Users\\anlas\\Documents\\Phd\\RNAseq_analysis\\outs\\g_24h_2.rds")

#load the other 24h dataset (1st replicate)
matrix_24h.data <- Read10X(data.dir = 'C:\\Users\\anlas\\Documents\\Phd\\RNAseq_analysis\\24h\\filtered_feature_bc_matrix', gene.column = 2, unique.features = TRUE)
g_24h <- CreateSeuratObject(counts = matrix_24h.data, project = "scRNA_gast24h", min.cells = 3, min.features = 200)
g_24h[["percent.mt"]] <- PercentageFeatureSet(g_24h, pattern = "^mt-")
#add information about replicate
g_24h$replicate <- "rep1"
#perform QC, normalization, and identification of most variable genes
g_24h <- subset(g_24h, subset = nFeature_RNA > 2750 & percent.mt < 15 & percent.mt > 1 & nCount_RNA < 100000)
VlnPlot(g_24h, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
g_24h <- NormalizeData(g_24h, normalization.method = "LogNormalize", scale.factor = 10000)
g_24h <- FindVariableFeatures(g_24h, selection.method = "vst", nfeatures = 3000)
#save
saveRDS(g_24h, file = "C:\\Users\\anlas\\Documents\\Phd\\RNAseq_analysis\\outs\\g_24h.rds")

#integrate 24h datasets
g_24hc.anchors <- FindIntegrationAnchors(object.list = list(g_24h, g_24h_2), dims = 1:50)
g_24hc.combined <- IntegrateData(anchorset = g_24hc.anchors, dims = 1:50)
#add ID
g_24hc.combined$merge.ident <- "24h"
#perform cell cycle correction on the integrated 24h dataset
#convert human genes into mouse orthologues
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
require("biomaRt")
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
convertHumanGeneList <- function(x){
  genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
  humanx <- unique(genesV2[, 2])
  # Print the first 6 genes found to the screen
  print(head(humanx))
  return(humanx)
}
m.s.genes <- convertHumanGeneList(cc.genes.updated.2019$s.genes)
m.g2m.genes <- convertHumanGeneList(cc.genes.updated.2019$g2m.genes)
#perform scoring
g_24hc.combined <- CellCycleScoring(g_24hc.combined, s.features = m.s.genes, g2m.features = m.g2m.genes, set.ident = TRUE)
#regress out difference between G2M and S phase scores (maintain differences in dividing versus non-dividing cells)
#this was chosen because differentiation in division activity may be meaningful or a result in the context of cell differentiation dynamics
g_24hc.combined$CC.Difference <- g_24hc.combined$S.Score - g_24hc.combined$G2M.Score
g_24hc.combined <- ScaleData(g_24hc.combined, vars.to.regress = "CC.Difference", features = rownames(g_24hc.combined))
#save
saveRDS(g_24hc.combined, file = "C:\\Users\\anlas\\Documents\\Phd\\RNAseq_analysis\\outs\\g_24hc.combined_forint.rds")

#####################
#load the 48h dataset (2nd replicate)
matrix_48h_2.data <- Read10X(data.dir = 'C:\\Users\\anlas\\Documents\\Phd\\RNAseq_analysis\\2020\\for_merging\\48h\\filtered_feature_bc_matrix', gene.column = 2, unique.features = TRUE)
g_48h_2 <- CreateSeuratObject(counts = matrix_48h_2.data, project = "scRNA_gast48h_2020", min.cells = 3, min.features = 200)
g_48h_2[["percent.mt"]] <- PercentageFeatureSet(g_48h_2, pattern = "^mt-")
#add information about replicate
g_48h_2$replicate <- "rep2_48h"
#perform QC, normalization, and identification of most variable genes
g_48h_2 <- subset(g_48h_2, subset = nFeature_RNA > 3500 & percent.mt < 15 & percent.mt > 1 & nCount_RNA < 120000)
VlnPlot(g_48h_2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
g_48h_2 <- NormalizeData(g_48h_2, normalization.method = "LogNormalize", scale.factor = 10000)
g_48h_2 <- FindVariableFeatures(g_48h_2, selection.method = "vst", nfeatures = 3000)
#save
saveRDS(g_48h_2, file = "C:\\Users\\anlas\\Documents\\Phd\\RNAseq_analysis\\outs\\g_48h_2.rds")

#load the 48h dataset (1st replicate)
matrix_48h.data <- Read10X(data.dir = 'C:\\Users\\anlas\\Documents\\Phd\\RNAseq_analysis\\48h\\filtered_feature_bc_matrix', gene.column = 2, unique.features = TRUE)
g_48h <- CreateSeuratObject(counts = matrix_48h.data, project = "scRNA_gast48h", min.cells = 3, min.features = 200)
g_48h[["percent.mt"]] <- PercentageFeatureSet(g_48h, pattern = "^mt-")
#add information about replicate
g_48h$replicate <- "rep1_48h"
#perform QC, normalization, and identification of most variable genes
g_48h <- subset(g_48h, subset = nFeature_RNA > 2750 & percent.mt < 15 & percent.mt > 1 & nCount_RNA < 100000)
VlnPlot(g_48h, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
g_48h <- NormalizeData(g_48h, normalization.method = "LogNormalize", scale.factor = 10000)
g_48h <- FindVariableFeatures(g_48h, selection.method = "vst", nfeatures = 3000)
#save
saveRDS(g_48h, file = "C:\\Users\\anlas\\Documents\\Phd\\RNAseq_analysis\\outs\\g_48h.rds")

#integrate 48h datasets
g_48hc.anchors <- FindIntegrationAnchors(object.list = list(g_48h, g_48h_2), dims = 1:50)
g_48hc.combined <- IntegrateData(anchorset = g_48hc.anchors, dims = 1:50)
#add ID
g_48hc.combined$merge.ident <- "48h"
#save
saveRDS(g_48hc.combined, file = "C:\\Users\\anlas\\Documents\\Phd\\RNAseq_analysis\\outs\\g_48hc.combined_forint.rds")

###########################
#integrate all datasets
all.anchors <- FindIntegrationAnchors(object.list = list(g_24hc.combined,g_48hc.combined,g72h_filtered), dims = 1:50)
all_data.combined <- IntegrateData(anchorset = all.anchors, dims = 1:50)

#scale data, perform PCA and clustering analyses, generate UMAP
DefaultAssay(all_data.combined) <- "integrated"
all.genes <- rownames(all_data.combined)
all_data.combined <- ScaleData(all_data.combined, features = all.genes)
all_data.combined <- RunPCA(all_data.combined, features = VariableFeatures(object = all_data.combined))
all_data.combined <- FindNeighbors(all_data.combined, dims = 1:50)
all_data.combined <- FindClusters(all_data.combined, resolution = 0.5)
all_data.combined <- RunUMAP(all_data.combined, dims = 1:50)

#generate UMAP plots (specify commands as desired)
p <- DimPlot(all_data.combined, reduction = "umap")
AugmentPlot(plot = p)
plot(p)
p <- DimPlot(all_data.combined, reduction = "umap", group.by = "merge.ident", cols = c('24h' = 'lightsteelblue3', '48h' = 'steelblue', '72h' = 'orchid4'), pt.size = 0.25)
AugmentPlot(plot = p)
plot(p)
#recolor clusters as seen in the publication
p <- DimPlot(all_data.combined, reduction = "umap",  cols = c('0' = 'skyblue2', '1' = 'slategrey', '2' = 'steelblue', '3' = 'skyblue4', '4' = 'royalblue4', '5' = 'lightsteelblue3', '6' = 'black', '7' = 'aquamarine4', '8' = 'indianred3', '9' = 'thistle3', '10' = 'lightcyan2'))
AugmentPlot(plot = p)
plot(p)
#save
saveRDS(all_data.combined, file = "C:\\Users\\anlas\\Documents\\Phd\\RNAseq_analysis\\outs\\all_data.combined_2.rds")

#generate a heatmap with top 5 differentially expressed genes for each cluster
all_data.combined.markers <- FindAllMarkers(all_data.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top5 <- all_data.combined.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
DoHeatmap(all_data.combined, features = top5$gene) + NoLegend() + scale_fill_viridis()

#calculate differentially expressed marker genes for each cluster
#change assay to RNA
DefaultAssay(all_data.combined) <- "RNA"
#calculate genes
all_data.combined.markers <- FindAllMarkers(all_data.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
#save the table
write.csv(all_data.combined.markers,"C:\\Users\\anlas\\Documents\\Phd\\RNAseq_analysis\\outs\\markers_clusters.csv", row.names = FALSE)

#generate plots to check expression of genes on UMAP
#change assay
DefaultAssay(all_data.combined) <- "RNA"
#set color code
plasma <- viridis(30, direction = 1, option = "D") 
#generate plot(s)
p <- FeaturePlot(all_data.combined, features = c("Sox1"), order = TRUE, cols = plasma, pt.size = 0.5) 
AugmentPlot(plot = p)
plot(p)

#generate plots for only at a certain time
Idents(object=all_data.combined) <- "merge.ident"
all_data.combined_24h <- subset(x = all_data.combined, idents = "24h")
all_data.combined_48h <- subset(x = all_data.combined, idents = "48h")
all_data.combined_72h <- subset(x = all_data.combined, idents = "72h")
DefaultAssay(all_data.combined_24h) <- "RNA"
DefaultAssay(all_data.combined_48h) <- "RNA"
DefaultAssay(all_data.combined_72h) <- "RNA"
#example for 24h (adapt as required)
p <- FeaturePlot(all_data.combined_24h, features = c("T"), order = TRUE, cols = plasma, pt.size = 0.5) 
AugmentPlot(plot = p)
plot(p)

#################
#perform GO term analysis
#load required packages
library(clusterProfiler)
library(org.Mm.eg.db)
#set assay to RNA
DefaultAssay(all_data.combined) <- "RNA" 
#calculate variable genes
all_data.combined_markers2 <- FindAllMarkers(all_data.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
#prepare data structure for GO term analysis, i.e. make separate lists for variable genes for each cluster
filtered_dataset = all_data.combined_markers2 %>%
  dplyr::select(cluster, gene)

cluster_genes_tot = filtered_dataset %>%
  dplyr::select(gene)
cluster_total = cluster_genes_tot$gene

cluster0 = filtered_dataset %>% filter(cluster == 0)
cluster0 = cluster0$gene

cluster1 = filtered_dataset %>% filter(cluster == 1)
cluster1 = cluster1$gene

cluster2 = filtered_dataset %>% filter(cluster == 2)
cluster2 = cluster2$gene

cluster3 = filtered_dataset %>% filter(cluster == 3)
cluster3 = cluster3$gene

cluster4 = filtered_dataset %>% filter(cluster == 4)
cluster4 = cluster4$gene

cluster5 = filtered_dataset %>% filter(cluster == 5)
cluster5 = cluster5$gene

cluster6 = filtered_dataset %>% filter(cluster == 6)
cluster6 = cluster6$gene

cluster7 = filtered_dataset %>% filter(cluster == 7)
cluster7 = cluster7$gene

cluster8 = filtered_dataset %>% filter(cluster == 8)
cluster8 = cluster8$gene

cluster9 = filtered_dataset %>% filter(cluster == 9)
cluster9 = cluster9$gene

cluster10 = filtered_dataset %>% filter(cluster == 10)
cluster10 = cluster10$gene

#assign ENTREZ-IDs
eg0 =bitr(cluster0, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
head(eg0)
eg0_en= eg0$ENTREZID

eg1 =bitr(cluster1, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
head(eg1)
eg1_en= eg1$ENTREZID

eg2 = bitr(cluster2, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
head(eg2)
eg2_en= eg2$ENTREZID

eg3 = bitr(cluster3, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
head(eg3)
eg3_en= eg3$ENTREZID

eg4 = bitr(cluster4, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
head(eg4)
eg4_en= eg4$ENTREZID

eg5 = bitr(cluster5, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
head(eg5)
eg5_en= eg5$ENTREZID

eg6 = bitr(cluster6, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
head(eg6)
eg6_en= eg6$ENTREZID

eg7 = bitr(cluster7, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
head(eg7)
eg7_en= eg7$ENTREZID

eg8 = bitr(cluster8, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
head(eg8)
eg8_en= eg8$ENTREZID

eg9 = bitr(cluster9, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
head(eg9)
eg9_en= eg9$ENTREZID

eg10 = bitr(cluster10, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
head(eg10)
eg10_en= eg10$ENTREZID

eg_universe = bitr(cluster_total, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
head(eg_universe)
eg_uni_en = eg_universe$ENTREZID

#perform GO term analysis
ego <- enrichGO(gene          = eg0_en,
                universe      = names(eg_uni_en),
                OrgDb         = org.Mm.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)
head(ego)

ego1 <- enrichGO(gene          = eg1_en,
                 universe      = names(eg_uni_en),
                 OrgDb         = org.Mm.eg.db,
                 ont           = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.01,
                 qvalueCutoff  = 0.05,
                 readable      = TRUE)
head(ego1)

ego2 <- enrichGO(gene          = eg2_en,
                 universe      = names(eg_uni_en),
                 OrgDb         = org.Mm.eg.db,
                 ont           = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.01,
                 qvalueCutoff  = 0.05,
                 readable      = TRUE)
head(ego2)

ego3 <- enrichGO(gene          = eg3_en,
                 universe      = names(eg_uni_en),
                 OrgDb         = org.Mm.eg.db,
                 ont           = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.01,
                 qvalueCutoff  = 0.05,
                 readable      = TRUE)
head(ego3)

ego4 <- enrichGO(gene          = eg4_en,
                 universe      = names(eg_uni_en),
                 OrgDb         = org.Mm.eg.db,
                 ont           = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.01,
                 qvalueCutoff  = 0.05,
                 readable      = TRUE)
head(ego4)

ego5 <- enrichGO(gene          = eg5_en,
                 universe      = names(eg_uni_en),
                 OrgDb         = org.Mm.eg.db,
                 ont           = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.01,
                 qvalueCutoff  = 0.05,
                 readable      = TRUE)
head(ego5)

ego6 <- enrichGO(gene          = eg6_en,
                 universe      = names(eg_uni_en),
                 OrgDb         = org.Mm.eg.db,
                 ont           = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.01,
                 qvalueCutoff  = 0.05,
                 readable      = TRUE)
head(ego6)

ego7 <- enrichGO(gene          = eg7_en,
                 universe      = names(eg_uni_en),
                 OrgDb         = org.Mm.eg.db,
                 ont           = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.01,
                 qvalueCutoff  = 0.05,
                 readable      = TRUE)
head(ego7)

ego8 <- enrichGO(gene          = eg8_en,
                 universe      = names(eg_uni_en),
                 OrgDb         = org.Mm.eg.db,
                 ont           = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.01,
                 qvalueCutoff  = 0.05,
                 readable      = TRUE)
head(ego8)

ego9 <- enrichGO(gene          = eg9_en,
                 universe      = names(eg_uni_en),
                 OrgDb         = org.Mm.eg.db,
                 ont           = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.01,
                 qvalueCutoff  = 0.05,
                 readable      = TRUE)
head(ego9)

ego10 <- enrichGO(gene         = eg10_en,
                  universe      = names(eg_uni_en),
                  OrgDb         = org.Mm.eg.db,
                  ont           = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.01,
                  qvalueCutoff  = 0.05,
                  readable      = TRUE)
head(ego10)

#save results
write.csv(ego@result, file = "C:\\Users\\anlas\\Desktop\\results_goterm\\c1.csv")
write.csv(ego1@result, file = "C:\\Users\\anlas\\Desktop\\results_goterm\\c2.csv")
write.csv(ego2@result, file = "C:\\Users\\anlas\\Desktop\\results_goterm\\c3.csv")
write.csv(ego3@result, file = "C:\\Users\\anlas\\Desktop\\results_goterm\\c4.csv")
write.csv(ego4@result, file = "C:\\Users\\anlas\\Desktop\\results_goterm\\c5.csv")
write.csv(ego5@result, file = "C:\\Users\\anlas\\Desktop\\results_goterm\\c6.csv")
write.csv(ego6@result, file = "C:\\Users\\anlas\\Desktop\\results_goterm\\c7.csv")
write.csv(ego7@result, file = "C:\\Users\\anlas\\Desktop\\results_goterm\\c8.csv")
write.csv(ego8@result, file = "C:\\Users\\anlas\\Desktop\\results_goterm\\c9.csv")
write.csv(ego9@result, file = "C:\\Users\\anlas\\Desktop\\results_goterm\\c10.csv")
write.csv(ego10@result, file = "C:\\Users\\anlas\\Desktop\\results_goterm\\c11.csv")

