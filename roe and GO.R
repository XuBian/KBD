library(Seurat)
library(scRNAtoolVis)
library(ggplot2)
library(openxlsx)
library(tidydr)
library(RColorBrewer)
library(Startrac)
library(tictoc)
library(ggpubr)
library(ComplexHeatmap)
library(circlize)
library(tidyverse)
library(sscVis)
library(pheatmap)

seurat <- readRDS("annotated_data.rds")
cell.levels <- c("InfC","preInfC","FC","preFC","HTC","preHTC","HomC","RepC","RegC","EC","ProC")
marker_gene<-c("GPR183","CD74","CXCL8","IFI27","IFI16","COL1A1","COL1A2","MMP2",
               "WWP2","COL27A1","PLCG2","COL10A1","IBSP","SPP1",
               "CRTAC1","ABI3BP","PRG4","HSPA1A","JUN","DDIT3","OGN","CILP","CILP2",
               "CHI3L2","CHI3L1","SAA1","CYTL1","FRZB","CHRDL2","HMGA1","BMP2","C11orf96")

#### 1. marker genes dot plot
pdf("markers.dotplot.pdf",width = 10,height = 5)
DotPlot(
  seurat, group.by = "celltype",
  features = marker_gene
) +
  scale_color_gradient2(low = "navy", mid = "white", high = "firebrick3") +
  theme(axis.text.x = element_text(size = 12, angle = 90, hjust = 1, face = "italic")) +
  theme(axis.text.y = element_text(size = 10)) +
  theme(legend.title = element_text(face = "bold", size = 12)) +
  xlab("") +
  ylab("")
dev.off()

#### 2. plot roe for celltype enrichment using STARTAC
meta.data <- seurat@meta.data
R_oe <- calTissueDist(meta.data,
                      byPatient = F,
                      colname.cluster = "celltype",
                      colname.patient = "sampleID",
                      colname.tissue = "group",
                      method = "chisq",
                      min.rowSum = 0)
data <- as.matrix(R_oe)
data <- t(data)
colors <- colorRampPalette(brewer.pal(n = 9, name = "YlOrBr"))(100)

symbol_matrix <- ifelse(data > 1, "+++",
                        ifelse(data > 0.8 & data <= 1, "++",
                               ifelse(data >= 0.2 & data <= 0.8, "+",
                                      ifelse(data > 0 & data < 0.2, "+/-",
                                             ifelse(data == 0, "-", "")))))
data <- as.matrix(data)
p <- pheatmap(data,
              color = colors,
              cluster_rows = FALSE,
              cluster_cols = FALSE,
              scale = "none",
              show_rownames = TRUE,
              show_colnames = TRUE,
              display_numbers = symbol_matrix,  
              number_color = "black",          
              fontsize_number = 10,            
              clustering_distance_rows = "euclidean",
              clustering_distance_cols = "euclidean",
              clustering_method = "complete",
              annotation_legend = TRUE,
              name = "Ro/e")
pdf("ROE.pdf",width = 6,height = 3)
p
dev.off()

#### 3. GO analysis
library(clusterProfiler)
library(org.Hs.eg.db)

GO_analysis=function(gene,universe=NULL,go_padj=0.05)
{
  id <- bitr(gene,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = "org.Hs.eg.db")
  geneList=id[,2]
  names(geneList)=as.character(geneList)
  if(!is.null(universe))
  {
    id2 <- bitr(universe,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = "org.Hs.eg.db")
    universe=id2[,2]
    names(universe)=as.character(universe)
  }
  GO <- enrichGO(OrgDb="org.Hs.eg.db", 
                 gene = names(geneList), 
                 pvalueCutoff = go_padj,
                 keyType = "ENTREZID",
                 pAdjustMethod = "BH", 
                 qvalueCutoff = 0.2, 
                 ont = "ALL",
                 minGSSize = 10, 
                 maxGSSize = 500, 
                 pool=T,
                 universe=universe,
                 readable=TRUE)
  return(GO)
}	
universe <- rownames(seurat)
celltype <- levels(seurat@meta.data$celltype)

GO.list <- list()
for (i in celltype){
  DE<-openxlsx::read.xlsx("markers.xlsx", sheet=i,rowNames=T)
  DE<-DE[DE$avg_log2FC>1&DE$p_val_adj<0.05,]
  DE<-DE[order(DE$avg_log2FC,decreasing = T),]
  gene=rownames(DE)
  go<-GO_analysis(gene,universe)
  go<- as.data.frame(go@result)
  GO.list[[i]] <- go
}

openxlsx::write.xlsx(GO.list,file='GO.xlsx')

#### 4. GO dotplot
## step 1. get top GO in each group enriched cell type
group_to_celltypes <- list(
  Control_EC = "EC", Control_RegC = "RegC",
  OA_preInfC = "preInfC", OA_InfC = "InfC",
  KBD_ProC = "ProC", KBD_HomC = "HomC"
)

go_sele_per_group <- list()

for (group_celltype in names(group_to_celltypes)) {
  group <- strsplit(group_celltype, "_")[[1]][1]  
  cell_type <- group_to_celltypes[[group_celltype]]  
  
  file_path <- paste0(group,".GO.xlsx")
  data <- tryCatch({
    openxlsx::read.xlsx(file_path, sheet = cell_type)
  }, error = function(e) {
    message(paste("Cannot read", group, cell_type, "sheet"))
    return(NULL)
  })
  
  if (!is.null(data)) {
    go_ids <- data$ID[1:10] 
    go_sele_per_group[[group_celltype]] <- unique(go_ids)
  }
}
## step 2. get matched data in all three groups
matched_data <- data.frame()
celltype_to_group <- list(
  EC = "Control_EC", RegC = "Control_RegC",
  preInfC = "OA_preInfC", InfC = "OA_InfC",
  ProC = "KBD_ProC", HomC = "KBD_HomC"
)
cell_types <- names(celltype_to_group)  
groups <- c("Control", "OA", "KBD")   

for (cell_type in cell_types) {
  matched_group_celltype <- celltype_to_group[[cell_type]]  
  
  for (group in groups) {
    file_path <- paste0(group, ".GO.xlsx")
    
    data <- tryCatch({
      openxlsx::read.xlsx(file_path, sheet = cell_type)
    }, error = function(e) {
      message(paste("Cannot read", group, cell_type, "sheet"))
      return(NULL)
    })
    
    if (!is.null(data)) {
      matched <- subset(data, ID %in% go_sele_per_group[[matched_group_celltype]]) 
      matched$Group <- group
      matched$CellType <- cell_type
      matched_data <- rbind(matched_data, matched)
    }
  }
}

head(matched_data)
matched_data <- as.data.frame(matched_data)

## 3. plot
res <- matched_data
res$logP <- -log10(res$p.adjust)
my.colors <- colorRampPalette(c("white", "#e3b0a3", "#d92523"))(100)

normalize <- function(x) {
  return((x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE)))
}

# normalize by group
res <- res %>%
  group_by(Description) %>%
  mutate(logP_norm_Descr = normalize(logP)) %>%
  ungroup()

res$Description <- stringr::str_to_sentence(res$Description)

res$Group <- factor(res$Group, levels = c( "KBD", "OA","Control"))
p <- ggplot(res, aes(x = Group, y = Description, size = Count, color = logP_norm_Descr)) +
  geom_point() +
  scale_color_gradientn(colors = my.colors) +  
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 15),  
    axis.text.y = element_text(size = 15),  
    legend.title = element_text(size = 12, face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    axis.line = element_line(colour = "black", size = 0.5), 
    axis.ticks = element_line(colour = "black", size = 0.5) 
  ) +
  labs(x = "Group", y = "Description", size = "Gene Count", color = "-logP_norm")
print(p)

pdf("GO_plot.pdf",width = 12,height = 12)
p
dev.off()

