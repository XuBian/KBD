suppressPackageStartupMessages({
library(doParallel)
library(Seurat)
library(SeuratWrappers)
library(harmony)
library(Rcpp)
library(ggplot2)
library(gridExtra)
library(dplyr)
})
rm(list = ls())
dim <- 30
resolution <- 0.9

# load data
cat("loading data")
ifnb_list <- readRDS("ifnb_list.rds")
combined <- merge(x=ifnb_list[[1]],y=ifnb_list[-1],project='Cart_KBD',add.cell.ids = names(ifnb_list),merge.data = T)
DefaultAssay(combined) <- "RNA"
# set parameters
features <- 2500

# run harmony
cat("running harmony")
combined <- NormalizeData(combined)
combined <- FindVariableFeatures(combined, selection.method = "vst", nfeatures = features) 
combined <- ScaleData(combined)
combined <- RunPCA(combined,reduction.name='pca',npcs = dim) 
combined <- RunHarmony(combined,reduction = 'pca',group.by.vars ='orig.ident', reduction.save = 'harmony')
combined <- RunUMAP(combined,reduction='harmony',dims=1:dim) 
combined <- FindNeighbors(combined, reduction = "harmony", dims = 1:dim)
combined <- FindClusters(combined, resolution = resolution)
p1 <- DimPlot(combined, group.by = "orig.ident",raster = F)
p2 <- DimPlot(combined, group.by = "seurat_clusters",label=T, label.size =10,raster = F )
combined_plot <- grid.arrange(p1, p2, ncol = 2)
ggsave(filename = paste0("dim_", dim, "_res_", resolution, ".harmony.umap.pdf"),
       plot = combined_plot, width = 14, height = 7)
saveRDS(combined, file = paste0("dim_",dim, "_res_",resolution,".harmony.rds"))

# Findmarkers
cat("running FindMarkers")

library(openxlsx)
count <- 1
markers_wilcox.list = list()
cat("running All FindMarkers")
for (i in levels(combined@meta.data$seurat_clusters)) {
  res <- FindMarkers(combined,slot="data",
                     ident.1=i,
                     ident.2=NULL,
                     logfc.threshold=0.5,
                     test.use='wilcox',
                     min.pct=0.25, 
                     min.diff.pct=-Inf, 
                     group.by='seurat_clusters')
  markers_wilcox.list[[count]] <- res[order(res$avg_log2FC, decreasing=TRUE), ]
  count <- count + 1
}## end for
names(markers_wilcox.list) <- levels(combined@meta.data$seurat_clusters)
openxlsx::write.xlsx(x = markers_wilcox.list, file =  paste0("dim_",dim, "_res_",resolution,".markers.xlsx"),row.names=TRUE)

