#load required packages
library(ggplot2)
library(cowplot)
library(dplyr)
library(Seurat)
library(Matrix)
library(viridis)
library(limma)

#load datasets
g_24hc.combined <- readRDS(file = "C:\\Users\\anlas\\Documents\\Phd\\RNAseq_analysis\\seurat_objects\\g_24hc.combined_forint.rds")
g_48hc.combined <- readRDS(file = "C:\\Users\\anlas\\Documents\\Phd\\RNAseq_analysis\\seurat_objects\\g_48hc.combined_forint.rds")
g72h_filtered <- readRDS(file = "C:\\Users\\anlas\\Documents\\Phd\\RNAseq_analysis\\seurat_objects\\g72h_filtered_forint.rds")
#set assays to RNA
DefaultAssay(g_24hc.combined) <- "RNA"
DefaultAssay(g_48hc.combined) <- "RNA"
DefaultAssay(g72h_filtered) <- "RNA"

#extract T+ and T- cells, identify differentially expressed genes and plot.
#example here for 48h dataset

#a ridgeplot can be generated to study distribution of expression of gene-of-interest
RidgePlot(g_48hc.combined,features = c("T"))
#extract T+ and T- cells and generate a new object with them, set assay to RNA for the new object
Tpos <- subset(g_48hc.combined, subset = T > 1.5)
Tneg <- subset(g_48hc.combined, subset = T <= 0.05)
g48h_tcomp <- merge(Tpos, y = Tneg, add.cell.ids = c("Tpos", "Tneg"), project = "48h_tcomp", merge.data = TRUE)
DefaultAssay(g48h_tcomp) <- "RNA"

#set Idents for the new object
Tpos <- subset(g48h_tcomp, subset = T > 1.5)
Tneg <- subset(g48h_tcomp, subset = T <= 0.05)
cells.use <- WhichCells(object = Tpos)
g48h_tcomp <- SetIdent(object = g48h_tcomp, cells = cells.use, value = 'Tpos')
cells.use2 <- WhichCells(object = Tneg)
g48h_tcomp <- SetIdent(object = g48h_tcomp, cells = cells.use2, value = 'Tneg')
head(x = Idents(object = g48h_tcomp))
DefaultAssay(g48h_tcomp) <- "RNA"

#identify marker genes and save as csv table
g48h_tcomp.markers <- FindAllMarkers(g48h_tcomp, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(g48h_tcomp.markers, file = "C:\\Users\\anlas\\Desktop\\new_dotplots\\48h_tpos-tneg_allmarkers.csv")
#extract top 20 differentially expressed genes
top20 <- g48h_tcomp.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
#generate dotplot
DotPlot(g48h_tcomp, features = top20$gene, cols = c("grey95","dodgerblue4"), scale = FALSE) + RotatedAxis()

#perform GO term analysis
#example here for 48h dataset
#load required packages
library(clusterProfiler)
library(org.Mm.eg.db)

#identify differentially expressed genes 
g48h_tcomp.markers2 <- FindAllMarkers(g48h_tcomp, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
filtered_dataset = g48h_tcomp.markers2 %>%
  dplyr::select(cluster, gene)

#prepare data structure, assign ENTREZ-IDs and run GO term analysis
cluster_genes_tot = filtered_dataset %>%
  dplyr::select(gene)
cluster_total = cluster_genes_tot$gene

cluster0 = filtered_dataset %>% filter(cluster == "Tpos")
cluster0 = cluster0$gene

cluster1 = filtered_dataset %>% filter(cluster == "Tneg")
cluster1 = cluster1$gene

eg0 =bitr(cluster0, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
head(eg0)
eg0_en= eg0$ENTREZID

eg1 =bitr(cluster1, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
head(eg1)
eg1_en= eg1$ENTREZID

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

#save results
write.csv(ego@result, file = "C:\\Users\\anlas\\Desktop\\new_dotplots\\48h_tpos_new.csv")
write.csv(ego1@result, file = "C:\\Users\\anlas\\Desktop\\new_dotplots\\48h_tneg_new.csv")

########################
#do the same as above but for Sox+ and Sox- cells
#example again for the 48h dataset

#extractSox+ and Sox- cells and generate a new object
RidgePlot(g_48hc.combined,features = c("Sox2"))
Sox2pos <- subset(g_48hc.combined, subset = Sox2 > 1.0)
Sox2neg <- subset(g_48hc.combined, subset = Sox2 <= 0.05)
g48h_sox2comp <- merge(Sox2pos, y = Sox2neg, add.cell.ids = c("Sox2pos", "Sox2neg"), project = "48h_sox2comp", merge.data = TRUE)
DefaultAssay(g48h_sox2comp) <- "RNA"

#set Idents
Sox2pos <- subset(g48h_sox2comp, subset = Sox2 > 1.0)
Sox2neg <- subset(g48h_sox2comp, subset = Sox2 <= 0.05)
cells.use <- WhichCells(object = Sox2pos)
g48h_sox2comp <- SetIdent(object = g48h_sox2comp, cells = cells.use, value = 'Sox2pos')
cells.use2 <- WhichCells(object = Sox2neg)
g48h_sox2comp <- SetIdent(object = g48h_sox2comp, cells = cells.use2, value = 'Sox2neg')
head(x = Idents(object = g48h_sox2comp))
DefaultAssay(g48h_sox2comp) <- "RNA"

#find variable genes and save the list as csv
g48h_sox2comp.markers <- FindAllMarkers(g48h_sox2comp, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top20 <- g48h_sox2comp.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
write.csv(g48h_sox2comp.markers, file = "C:\\Users\\anlas\\Desktop\\new_dotplots\\48h_sox2posneg_allmarkers.csv")

#generate dotplot
DotPlot(g48h_sox2comp, features = top20$gene, cols = c("grey95","dodgerblue4"), scale = FALSE) + RotatedAxis()

#perform GO term analysis
#example here for 48h dataset
#load required packages
g48h_sox2comp.markers2 <- FindAllMarkers(g48h_sox2comp, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
filtered_dataset = g48h_sox2comp.markers2 %>%
  dplyr::select(cluster, gene)

cluster_genes_tot = filtered_dataset %>%
  dplyr::select(gene)
cluster_total = cluster_genes_tot$gene

cluster0 = filtered_dataset %>% filter(cluster == "Sox2pos")
cluster0 = cluster0$gene

cluster1 = filtered_dataset %>% filter(cluster == "Sox2neg")
cluster1 = cluster1$gene

eg0 =bitr(cluster0, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
head(eg0)
eg0_en= eg0$ENTREZID

eg1 =bitr(cluster1, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
head(eg1)
eg1_en= eg1$ENTREZID

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

#save results
write.csv(ego@result, file = "C:\\Users\\anlas\\Desktop\\new_dotplots\\48h_sox2pos_new.csv")
write.csv(ego1@result, file = "C:\\Users\\anlas\\Desktop\\new_dotplots\\48h_sox2neg_new.csv")