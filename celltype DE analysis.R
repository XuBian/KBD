library(Seurat)
library(cacoa)
library(tidyverse)
library(magrittr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggrepel)
library(ggplot2)
library(R.utils)
library(reshape2)

seurat <- readRDS("annotated_data.rds")

#### 1. differential expressed gene analysis using cacoa, KBD vs. OA as an example
sample.groups=unique(seurat@meta.data$orig.ident)
sample.groups=as.character(sample.groups)
names(sample.groups)=sample.groups
sample.groups[c(1:8)]='OA'
sample.groups[c(9:14)]='KBD'
cell.groups=seurat@meta.data$celltype
names(cell.groups)=colnames(seurat)
sample.per.cell=seurat@meta.data$orig.ident
names(sample.per.cell)=colnames(seurat)
ref.level='OA'
target.level='KBD'
cm=seurat@assays$RNA@counts
embedding=seurat@reductions$umap@cell.embeddings
rownames(embedding)=colnames(cm)
cao<-Cacoa$new(as.matrix(cm),sample.groups=sample.groups, cell.groups=cell.groups, sample.per.cell=sample.per.cell, 
               ref.level=ref.level, target.level=target.level, embedding=embedding)

cao$plot.params <- list(size=0.1, alpha=0.1, font.size=c(2, 3))
cao$plot.theme <- cao$plot.theme + theme(legend.background=element_blank())

possible.tests <- c('DESeq2.Wald', 'DESeq2.LRT', 
                    'Wilcoxon.edgeR','t-test.edgeR',
                    'Wilcoxon.DESeq2', 't-test.DESeq2',
                    'limma-voom')
possible.resampling <- c('loo')
possible.tests.loo = paste(possible.tests, possible.resampling, sep = '.')   

f.changes = FALSE
for(test in possible.tests){
  # message(test)
  for(resampling in possible.resampling){
    name <- paste(test, resampling, sep='.')
    if(!(name %in% names(cao$test.results))){
      print(name)
      cao$estimateDEPerCellType(test = test, resampling.method = resampling,
                                name = name)
      f.changes = TRUE
    }
  }
}  

saveRDS(cao,"KBDvsOA.cao.rds")


#### 2. plot differential expressed gene 
cell.levels <- c("preFC", "HTC", "preHTC", "HomC", "RepC", "RegC", "EC", "ProC")
cell.color <- c("#F8C7B4", "#DF5406", "#2EB9CC", "#C657B3", "#ADD187", "#159931", "#CB2B1F", "#E99407")
cell.color <- setNames(cell.color, cell.levels)

# Function to process differential expression results
process_DE_results <- function(cao_data) {
  DE.list <- list()
  sub_levels <- rev(c("preFC", "HTC", "preHTC", "HomC", "RepC", "RegC", "EC", "ProC"))
  
  for (i in sub_levels) {
    DE <- as.data.frame(cao_data$test.results$DESeq2.Wald.loo[[i]]$res)
    DE <- DE[order(DE$log2FoldChange, decreasing = TRUE), ]
    DE <- DE[which(abs(DE$log2FoldChange) > 1 & DE$padj < 0.05), ]
    DE.list[[i]] <- DE
  }
  
  return(DE.list)
}

# Load and process each comparison
KBD_OA.DE <- process_DE_results(readRDS("KBDvsOA.cao.rds"))
KBD_CON.DE <- process_DE_results(readRDS("KBDvsControl.cao.rds"))
OA_CON.DE <- process_DE_results(readRDS("OAvsControl.cao.rds"))

# Function to summarize gene counts for upregulated and downregulated genes
summarize_DE <- function(DE_list) {
  up <- sapply(DE_list, function(DE) length(DE$Gene[DE$log2FoldChange > 0]))
  down <- sapply(DE_list, function(DE) length(DE$Gene[DE$log2FoldChange < 0]))
  return(list(up = up, down = down))
}

# Summarize DE counts for each comparison
KBD_OA_counts <- summarize_DE(KBD_OA.DE)
KBD_CON_counts <- summarize_DE(KBD_CON.DE)
OA_CON_counts <- summarize_DE(OA_CON.DE)

# Combine data into long format for plotting
create_bar_data <- function(up, down, comparison_name) {
  bar_data <- data.frame(
    CellType = rep(names(up), 2),
    Direction = rep(c(paste(comparison_name, "Up", sep = "_"), paste(comparison_name, "Down", sep = "_")), each = length(up)),
    GeneCount = c(up, -down)
  )
  bar_data_long <- melt(bar_data, id.vars = c("CellType", "Direction"), variable.name = "Category", value.name = "GeneCount")
  return(bar_data_long)
}

# Prepare data for plotting
bar_data_KBD_CON <- create_bar_data(KBD_CON_counts$up, KBD_CON_counts$down, "KBD_CON")
bar_data_KBD_OA <- create_bar_data(KBD_OA_counts$up, KBD_OA_counts$down, "KBD_OA")
bar_data_OA_CON <- create_bar_data(OA_CON_counts$up, OA_CON_counts$down, "OA_CON")

# Function to create the plot
create_plot <- function(bar_data_long, title) {
  bar_data_long$CellType <- factor(bar_data_long$CellType, levels = rev(cell.levels))
  
  ggplot(bar_data_long, aes(x = CellType, y = GeneCount, fill = Direction)) +
    geom_bar(stat = "identity", position = "stack") +
    geom_text(aes(label = paste0 ("# ",abs(GeneCount)), 
                  y = ifelse(GeneCount > 0, GeneCount + 15, GeneCount - 15)),
              size = 3.5, fontface = "bold") +
    
    labs(x = "Cell Types", y = "Number of Genes (Upregulated vs Downregulated)", title = title) +
    theme_minimal(base_size = 14) + 
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(color = "black"),
      axis.ticks = element_line(color = "black"),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
      axis.text.y = element_text(size = 14),
      axis.title.x = element_text(size = 16, face = "bold"),
      axis.title.y = element_text(size = 16, face = "bold")
    ) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black")
}


# Plot and save
pdf("KBD_CON_DE_count_bar.pdf", width = 6, height = 4)# adjust
print(
  create_plot(bar_data_KBD_CON, "Number of Upregulated and Downregulated Genes in KBD vs Control") +
    scale_fill_manual(values = c( "#7AC3DF","#CC79A7"))
)
dev.off()

pdf("KBD_OA_DE_count_bar.pdf", width = 6, height = 4)# adjust
print(
  create_plot(bar_data_KBD_OA, "Number of Upregulated and Downregulated Genes in KBD vs OA") +
    scale_fill_manual(values = c("#EB7E60","#CC79A7"))
)
dev.off()

pdf("OA_CON_DE_count_bar.pdf", width = 6, height = 4)# adjust
print(
  create_plot(bar_data_OA_CON, "Number of Upregulated and Downregulated Genes in OA vs Control") +
    scale_fill_manual(values = c("#66cc66","#EB7E60"))
)
dev.off()
