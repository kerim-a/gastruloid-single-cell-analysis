#load required packages
library(ggplot2)
library(cowplot)
library(dplyr)
library(Seurat)
library(Matrix)
library(viridis)
library(limma)

#load datasets
g0h_filtered_cc <- readRDS(file = "C:\\Users\\anlas\\Documents\\Phd\\RNAseq_analysis\\seurat_objects\\g0h_filtered_cc.rds")
g_24h_2 <- readRDS(file = "C:\\Users\\anlas\\Documents\\Phd\\RNAseq_analysis\\outs\\g_24h_2.rds")
g_24h <- readRDS(file = "C:\\Users\\anlas\\Documents\\Phd\\RNAseq_analysis\\outs\\g_24h.rds")
g_48h_2 <- readRDS(file = "C:\\Users\\anlas\\Documents\\Phd\\RNAseq_analysis\\outs\\g_48h_2.rds")
g_48h <- readRDS(file = "C:\\Users\\anlas\\Documents\\Phd\\RNAseq_analysis\\outs\\g_48h.rds")
g72h_filtered <- readRDS(file = "C:\\Users\\anlas\\Documents\\Phd\\RNAseq_analysis\\seurat_objects\\g72h_filtered_forint.rds")

#merge 24h and 48h replicates
g24h_merged <- merge(g_24h, y = g_24h_2, add.cell.ids = c("24h_1", "24h_2"), project = "24h_merge", merge.data = TRUE)
g48h_merged <- merge(g_48h, y = g_48h_2, add.cell.ids = c("48h_1", "48h_2"), project = "48h_merge", merge.data = TRUE)
#add ID
g24h_merged$merge.ident <- "24h"
g48h_merged$merge.ident <- "48h"

#merge all datasets
alldata_merged <- merge(g0h_filtered_cc, y = c(g24h_merged, g48h_merged, g72h_filtered), add.cell.ids = c("0h", "24h", "48h", "72h"), project = "alldata_merge", merge.data = TRUE)
#set assay to RNA
DefaultAssay(alldata_merged) <- "RNA"
#save
saveRDS(alldata_merged, file = "C:\\Users\\anlas\\Documents\\Phd\\RNAseq_analysis\\seurat_objects\\alldata_merged_0-72h_oldpaper_new.rds")
alldata_merged <- readRDS(file = "C:\\Users\\anlas\\Documents\\Phd\\RNAseq_analysis\\seurat_objects\\alldata_merged_0-72h_oldpaper_new.rds")

#make gene lists
list_mesoderm <- c("T","Eomes","Gsc","Wnt3","Wnt3a","Mixl1","Snai1","Cer1","Lefty1","Meox1","Tbx6","Mesp1","Hes7","Msgn1","Foxc2","Osr1","Six2","Lhx1","Hand1","Hand2","Myl7","Pdgfra","Kdr")
list_endoderm <- c("Sox17","Foxa2","Gata4")
list_ectoderm <- c("Sox1","Sox2","Sox3","Sox11","Ncam1","Epha5","Ptn","Prtg","Otx2","Zic1","Neurog2","Pax2","Pax6","Nes","Pou3f1","Utf1")
list_pluripotent <- c("Zfp42","Sox2","Nanog","Esrrb","Klf4")
list_mechanics <- c("Cdh1","Cdh2","Cdh11","Krt8","Krt18","Fn1","Vim")
list_primedpluri <- c("Pim2","Dnmt3b","Pou3f1","Fgf5","Fst")

#generate dotplots (optionally flip coordinates)
DotPlot(alldata_merged, features = list_mesoderm, cols = c("lightgrey","dodgerblue4"), group.by = "merge.ident", scale.by = "size", scale = TRUE) #+ coord_flip()
DotPlot(alldata_merged, features = list_endoderm, cols = c("lightgrey","dodgerblue4"), group.by = "merge.ident", scale.by = "size", scale = TRUE) #+ coord_flip()
DotPlot(alldata_merged, features = list_pluripotent, cols = c("lightgrey","dodgerblue4"), group.by = "merge.ident", scale.by = "size", scale = TRUE) #+ coord_flip()
DotPlot(alldata_merged, features = list_ectoderm, cols = c("lightgrey","dodgerblue4"), group.by = "merge.ident", scale.by = "size", scale = TRUE) #+ coord_flip()
DotPlot(alldata_merged, features = list_mechanics, cols = c("lightgrey","dodgerblue4"), group.by = "merge.ident", scale.by = "size", scale = TRUE) #+ coord_flip()
DotPlot(alldata_merged, features = list_primedpluri, cols = c("lightgrey","dodgerblue4"), group.by = "merge.ident", scale.by = "size", scale = TRUE) #+ coord_flip()
