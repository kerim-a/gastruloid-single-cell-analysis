library(ggplot2)
library(monocle3)
library(Seurat)
library(SeuratWrappers)
library(tidyverse)
library(viridis)
library(reshape2)
library(tidyr)

#load seurat object and convert to cell_data_set
all_data.combined <- readRDS(file = "C:\\Users\\anlas\\Documents\\Phd\\RNAseq_analysis\\seurat_objects\\all_data.combined_2_integrated.rds")
DefaultAssay(all_data.combined) <- "RNA"

#annotate cds dataset, transfer gene names, cluster partitions etc
fData(cds)$gene_short_name <- rownames(fData(cds))

recreate.partition <- c(rep(1,length(cds@colData@rownames)))
names(recreate.partition) <- cds@colData@rownames
recreate.partition <- as.factor(recreate.partition)
cds@clusters$UMAP$partitions <- recreate.partition

list_cluster <- all_data.combined@active.ident
cds@clusters$UMAP$clusters <- list_cluster

cds@int_colData@listData$reducedDims$UMAP <- all_data.combined@reductions$umap@cell.embeddings

#build trajectory
cds <- learn_graph(cds, use_partition = TRUE, close_loop = TRUE)

plot_cells(cds,
                color_cells_by = 'cluster',
                label_groups_by_cluster = FALSE,
                trajectory_graph_color = "black",
                trajectory_graph_segment_size = 1.5,
                label_branch_points = FALSE,
                label_principal_points = FALSE,
                label_roots = FALSE,
                label_leaves = FALSE,
                group_label_size = 5)

#order cells, the root cells are the most pluripotent cluster in this case.
cds <- order_cells(cds, reduction_method = 'UMAP', root_cells = colnames(cds[,clusters(cds) == 6])) 

plot_cells(cds,
                color_cells_by = 'pseudotime',
                label_groups_by_cluster = FALSE,
                trajectory_graph_color = "black",
                trajectory_graph_segment_size = 1.5,
                label_branch_points = FALSE,
                label_roots = FALSE,
                label_leaves = FALSE,
                group_label_size = 5)

#transfer pseudotime information into seurat
all_data.combined$pseudotime <- pseudotime(cds)
FeaturePlot(all_data.combined, features = 'pseudotime', cols = c("lightgrey", "dodgerblue4"), pt.size = 0.5)

#plot line profiles from genes of interest
#make more plots, save as 7,3 size
for_gg <- FetchData(object = all_data.combined, vars = c("pseudotime","T","Wnt3a"))
for_gg_new <- gather(for_gg, gene, expression, T:Wnt3a, factor_key=TRUE)

ggplot(for_gg_new, aes(pseudotime, expression))+
  geom_smooth(aes(color = gene, fill = gene)) + 
  scale_color_viridis(discrete = TRUE, option = "D") +
  scale_fill_viridis(discrete = TRUE) +
  theme_light()
