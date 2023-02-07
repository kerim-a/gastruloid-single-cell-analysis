#load required packages
library(ggplot2)
library(cowplot)
library(dplyr)
library(Seurat)
library(Matrix)
library(viridis)
library(limma)

#load dataset
g72h_filtered <- readRDS(file = "C:\\Users\\anlas\\Documents\\Phd\\RNAseq_analysis\\seurat_objects\\g72h_filtered_forint.rds")

#scale data, run PCA and generate a UMAP
all.genes <- rownames(g72h_filtered)
g72h_filtered <- ScaleData(g72h_filtered, features = all.genes)
g72h_filtered <- RunPCA(g72h_filtered, features = VariableFeatures(object = g72h_filtered))
g72h_filtered <- FindNeighbors(g72h_filtered, dims = 1:50)
g72h_filtered <- FindClusters(g72h_filtered, resolution = 0.7)
g72h_filtered <- RunUMAP(g72h_filtered, dims = 1:50)

#generate a UMAP plot
p <- DimPlot(g72h_filtered, reduction = "umap")
AugmentPlot(plot = p)
plot(p)

#set color code for plots
plasma <- viridis(30, direction = 1, option = "D") 
#generate plot(s) to check expression of select marker genes on UMAP space
p <- FeaturePlot(g72h_filtered, features = c("Sox1"), order = TRUE, cols = plasma, pt.size = 0.5) 
AugmentPlot(plot = p)
plot(p)

#identify most variable genes and save results
g72h_filtered.markers <- FindAllMarkers(g72h_filtered, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(g72h_filtered.markers, file = "C:\\Users\\anlas\\Desktop\\72h_filtered_allmarkers.csv")

#extract top 5 differentially expressed genes per cluster and generate a heatmap
top5 <- g72h_filtered.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
DoHeatmap(g72h_filtered, features = top5$gene) + NoLegend() + scale_fill_viridis()

#save
saveRDS(g72h_filtered, file = "C:\\Users\\anlas\\Documents\\Phd\\RNAseq_analysis\\outs\\g72h_filtered_forint.rds")

# perform GO-term analysis to further characterize clusters
library(clusterProfiler)
library(org.Mm.eg.db)

g72h_filtered.markers2 <- FindAllMarkers(g72h_filtered, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
filtered_dataset = g72h_filtered.markers2 %>%
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

eg_universe = bitr(cluster_total, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
head(eg_universe)
eg_uni_en = eg_universe$ENTREZID

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

#save results
write.csv(ego@result, file = "C:\\Users\\anlas\\Desktop\\results_goterm72h\\c1.csv")
write.csv(ego1@result, file = "C:\\Users\\anlas\\Desktop\\results_goterm72h\\c2.csv")
write.csv(ego2@result, file = "C:\\Users\\anlas\\Desktop\\results_goterm72h\\c3.csv")
write.csv(ego3@result, file = "C:\\Users\\anlas\\Desktop\\results_goterm72h\\c4.csv")
write.csv(ego4@result, file = "C:\\Users\\anlas\\Desktop\\results_goterm72h\\c5.csv")
write.csv(ego5@result, file = "C:\\Users\\anlas\\Desktop\\results_goterm72h\\c6.csv")
write.csv(ego6@result, file = "C:\\Users\\anlas\\Desktop\\results_goterm72h\\c7.csv")
write.csv(ego7@result, file = "C:\\Users\\anlas\\Desktop\\results_goterm72h\\c8.csv")
write.csv(ego8@result, file = "C:\\Users\\anlas\\Desktop\\results_goterm72h\\c9.csv")
write.csv(ego9@result, file = "C:\\Users\\anlas\\Desktop\\results_goterm72h\\c10.csv")


###########################
#load the 0h (mESC) data and perform QC
matrix_0h.data <- Read10X(data.dir = 'C:\\Users\\anlas\\Documents\\Phd\\RNAseq_analysis\\2020\\0h\\filtered_feature_bc_matrix', gene.column = 2, unique.features = TRUE)
g_0h <- CreateSeuratObject(counts = matrix_0h.data, project = "scRNA_gast0h_2020", min.cells = 3, min.features = 200)
g_0h[["percent.mt"]] <- PercentageFeatureSet(g_0h, pattern = "^mt-")
g0h_filtered <- subset(g_0h, subset = nFeature_RNA > 2500 & percent.mt < 10 & nCount_RNA < 150000)

#normalize data, identify most variable gene and perform data scaling
g0h_filtered <- NormalizeData(g0h_filtered, normalization.method = "LogNormalize", scale.factor = 10000)
g0h_filtered <- FindVariableFeatures(g0h_filtered, selection.method = "vst", nfeatures = 3000)
all.genes <- rownames(g0h_filtered)
g0h_filtered <- ScaleData(g0h_filtered, features = all.genes)

#add ID for merging with other data
g0h_filtered$merge.ident <- "0h"

#run PCA, clustering and generate a UMAP plot
g0h_filtered <- RunPCA(g0h_filtered, features = VariableFeatures(object = g0h_filtered))
g0h_filtered <- FindNeighbors(g0h_filtered, dims = 1:50)
g0h_filtered <- FindClusters(g0h_filtered, resolution = 0.5)
g0h_filtered <- RunUMAP(g0h_filtered, dims = 1:50)
p <- DimPlot(g0h_filtered, reduction = "umap")
AugmentPlot(plot = p)
plot(p)

#perform cell cycle correction
#convert human gene lists into mouse orthologues
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

#cell cycle scoring
g0h_filtered <- CellCycleScoring(g0h_filtered, s.features = m.s.genes, g2m.features = m.g2m.genes, set.ident = TRUE)
#regress out cell cycle components
g0h_filtered <- ScaleData(g0h_filtered, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(g0h_filtered))

#run PCA and re-generate a UMAP plot with the re-scaled data
g0h_filtered <- RunPCA(g0h_filtered, features = VariableFeatures(object = g0h_filtered))
g0h_filtered <- FindNeighbors(g0h_filtered, dims = 1:50)
g0h_filtered <- FindClusters(g0h_filtered, resolution = 0.5)
g0h_filtered <- RunUMAP(g0h_filtered, dims = 1:50)
p <- DimPlot(g0h_filtered, reduction = "umap")
AugmentPlot(plot = p)
plot(p)

#save seurat object
saveRDS(g0h_filtered, file = "C:\\Users\\anlas\\Documents\\Phd\\RNAseq_analysis\\outs\\g0h_filtered_cc.rds")

#identify differentially expressed genes and save gene list
g0h_filtered.markers <- FindAllMarkers(g0h_filtered, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(g0h_filtered.markers,"C:\\Users\\anlas\\Documents\\Phd\\Papers_experimental\\For_paper_scRNAseq\\r_4_0\\top30_0h_alt2_2.csv", row.names = TRUE)

#plot a heatmap with the top 5 differentially expressed genes
top5 <- g0h_filtered.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
DoHeatmap(g0h_filtered, features = top5$gene) + NoLegend() + scale_fill_viridis()

#further characterize clusters via GO term analysis
library(clusterProfiler)
library(org.Mm.eg.db)

g0h_filtered.markers2 <- FindAllMarkers(g0h_filtered, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
filtered_dataset = g0h_filtered.markers2 %>%
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


eg_universe = bitr(cluster_total, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
head(eg_universe)
eg_uni_en = eg_universe$ENTREZID

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
#save results
write.csv(ego@result, file = "C:\\Users\\anlas\\Desktop\\results_goterm0h\\c1.csv")
write.csv(ego1@result, file = "C:\\Users\\anlas\\Desktop\\results_goterm0h\\c2.csv")
write.csv(ego2@result, file = "C:\\Users\\anlas\\Desktop\\results_goterm0h\\c3.csv")
write.csv(ego3@result, file = "C:\\Users\\anlas\\Desktop\\results_goterm0h\\c4.csv")
write.csv(ego4@result, file = "C:\\Users\\anlas\\Desktop\\results_goterm0h\\c5.csv")
write.csv(ego5@result, file = "C:\\Users\\anlas\\Desktop\\results_goterm0h\\c6.csv")
write.csv(ego6@result, file = "C:\\Users\\anlas\\Desktop\\results_goterm0h\\c7.csv")
