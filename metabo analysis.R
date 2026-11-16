library(Seurat)
library(scRNAtoolVis)
library(ggplot2)
library(openxlsx)
library(tidydr)
library(ggrepel)
library(RColorBrewer)
library(clusterProfiler)
library(org.Hs.eg.db)
library(pheatmap)
library(fgsea)
library(ComplexHeatmap)
library(matrixStats)
library(dplyr)

seurat <- readRDS("annotated_data.rds")
cell.levels <- c("InfC","preInfC","FC","preFC","HTC","preHTC","HomC","RepC","RegC","EC","ProC")
cell.color <- c("#717200","#763931","#2561A0","#F8C7B4","#DF5406","#2EB9CC","#C657B3","#ADD187",
                "#159931","#CB2B1F","#E99407")
cell.color <- setNames(cell.color,cell.levels)
group.levels <- c("Control","OA","KBD")
group.color <- c("#469393","#F1A93B","#C24976")
group.color <- setNames(group.color,group.levels)

#### 1. correlation analysis of metabolism-related genes
mgene=gmtPathways("KEGG_metabolism.gmt")
mgene<- unique(as.vector(unname(unlist(mgene))))
meta.genes.present <- intersect(mgene, rownames(seurat))  
seurat <- subset(seurat, features = meta.genes.present)  
seurat$group_celltype <- paste0(seurat$group,"_",seurat$celltype)

calculate_celltype_correlation <- function(seurat_obj) {
  expr_matrix <- GetAssayData(seurat_obj, slot = "data") 
  celltypes <- seurat_obj$group_celltype
  
  avg_expr <- sapply(unique(celltypes), function(ct) {
    cells <- colnames(seurat_obj)[celltypes == ct]
    rowMeans(expr_matrix[, cells, drop = FALSE])
  })
  colnames(avg_expr) <- unique(celltypes)
  cor_matrix <- cor(avg_expr, method = "spearman")
  return(cor_matrix)
}

cor_data <- calculate_celltype_correlation(seurat)

# normalize correlation data
normalize_0_1 <- function(mat) {
  rng <- range(mat, na.rm = TRUE)   
  (mat - rng[1]) / (rng[2] - rng[1])
}

cor_data_norm <- normalize_0_1(cor_data)
range(cor_data_norm)

# add annotation
group_labels <- factor(sub("_.*", "", order), levels = group.levels, ordered = TRUE)  
celltype_labels <- factor(sub(".*_", "", order), levels = rev(cell.levels), ordered = TRUE)  
row_ha <- rowAnnotation(
  group = anno_simple(group_labels, col = group.color),
  celltype = anno_simple(celltype_labels, col = cell.color)
)
col_ha <- HeatmapAnnotation(
  group = group_labels, 
  celltype = celltype_labels,
  col = list(group = group.color, celltype = cell.color),
  annotation_name_side = "left"  
)
my.colors <- colorRampPalette(c("#2c3385", "#9499c0", "white", "#e3b0a3", "#d92523"))(100)

# plot
ht <- Heatmap(
  cor_data_norm,
  name = "Correlation",  
  col = my.colors,  
  cluster_rows = FALSE, 
  cluster_columns = FALSE,
  show_row_names = TRUE, 
  show_column_names = TRUE,
  row_names_side = "left", 
  column_names_side = "top",
  row_dend_side = "left",
  column_dend_side = "top",
  top_annotation = col_ha,  
  left_annotation = row_ha   
)
pdf("group_celltype.spearman.correlation.pdf", width = 12, height = 10)#adjust
print(ht)
dev.off()

#### 2. metabolic activity analysis, take ProC metabolic activity as an example
## 1. setp 1, prepare pathway gene weights
pathway_gene <- list()
gene_pathways <- list()

gmt_file <- "KEGG_metabolism.gmt"
gmt_lines <- readLines(gmt_file)

for (line in gmt_lines) {
  parts <- unlist(strsplit(line, "\t"))
  pathway_name <- parts[1]
  genes <- parts[3:length(parts)] %>% 
    .[. != ""] %>%    
    unique()          
  
  pathway_gene[[pathway_name]] <- genes
  
  for (gene in genes) {
    gene_pathways[[gene]] <- unique(c(gene_pathways[[gene]], pathway_name))
  }
}

gene_weight <- lapply(gene_pathways, function(x) 1/length(x))

## step 2: load data and construct results matrix
cat("load data and construct results matrix")
idx1 <- which(seurat@meta.data$celltype=="ProC")
seurat <- seurat[,idx1]
all_cell_types <- seurat@meta.data$group 
cell_types <- unique(all_cell_types)
pathways <- pathway_gene 
pathway_names <- names(pathways)
gene_pathway_number <- unlist(gene_weight) 
# load data
norm_tpm <- seurat@assays$RNA@data
# construct result matrix
cat("construct result matrix")
mean_expression_shuffle <- matrix(NA, nrow=length(pathway_names), ncol=length(cell_types), 
                                  dimnames = list(pathway_names, cell_types))
mean_expression_noshuffle <- matrix(NA, nrow=length(pathway_names), ncol=length(cell_types),
                                    dimnames = list(pathway_names, cell_types))
pvalues_mat <- matrix(NA, nrow=length(pathway_names), ncol=length(cell_types),
                      dimnames = list(pathway_names, cell_types))

# calculate pathway activity
cat("calculate pathway activity")
set.seed(123)
for(p in pathway_names){
  genes <- pathways[[p]]
  genes_comm <- intersect(genes, rownames(norm_tpm))
  if(length(genes_comm) < 5) next 
  
  # get expression matrix
  pathway_expr <- norm_tpm[genes_comm, ]
  pathway_expr <- pathway_expr[rowSums(pathway_expr) > 0, ]
  
  # calculate mean gene expression in all groups (celltypes) 
  mean_exp <- apply(pathway_expr, 1, function(x) by(x, all_cell_types, mean))
  keep <- colnames(mean_exp)[colAlls(mean_exp > 0.001)]
  
  if(length(keep) < 3) next
  pathway_expr <- pathway_expr[keep, ]
  
  # replace 0 using min(non-zero) values to avoid extreme values
  pathway_expr <- t(apply(pathway_expr, 1, function(x) {
    x[x <= 0] <- min(x[x > 0])
    x
  }))
  
  # get gene weights in each pathway
  weights <- gene_pathway_number[keep]
  
  # calculate relative expression 
  mean_exp <- apply(pathway_expr, 1, function(x) by(x, all_cell_types, mean))
  ratio_exp <- t(mean_exp) / colMeans(mean_exp)
  
  # remove abnormal values
  col_quantile <- apply(ratio_exp, 2, quantile, na.rm=TRUE)
  col_upper <- col_quantile["75%", ] * 3
  col_lower <- col_quantile["25%", ] / 3
  outliers <- apply(ratio_exp, 1, function(x) any((x > col_upper) | (x < col_lower)))
  
  if(sum(!outliers) < 3) next
  keep <- names(outliers)[!outliers]
  pathway_expr <- pathway_expr[keep, ]
  weights <- gene_pathway_number[keep]
  
  # calculate weighted average expression
  mean_exp <- apply(pathway_expr, 1, function(x) by(x, all_cell_types, mean))
  ratio_exp <- t(mean_exp) / colMeans(mean_exp)
  pathway_activity <- apply(ratio_exp, 2, function(x) weighted.mean(x, weights/sum(weights)))
  
  # Random permutation test
  cat("Random permutation test")
  shuffle_times <- 10000
  
  mean_expression_shuffle[p, ] <- pathway_activity[cell_types]
  mean_expression_noshuffle[p, ] <- pathway_activity[cell_types]
  
  weight_values <- weights/sum(weights)
  shuffle_results <- matrix(NA, nrow=shuffle_times, ncol=length(cell_types))
  colnames(shuffle_results) <- cell_types
  
  for(i in 1:shuffle_times){
    # shuffle celltype labels
    shuffled_types <- sample(all_cell_types)
    shuffle_exp <- apply(pathway_expr, 1, function(x) by(x, shuffled_types, mean))
    shuffle_ratio <- t(shuffle_exp) / colMeans(shuffle_exp)
    shuffle_results[i, ] <- apply(shuffle_ratio, 2, function(x) weighted.mean(x, weight_values))
  }
  
  # calculate p values
  cat("calculate p values")
  for(c in cell_types){
    if(is.na(mean_expression_shuffle[p, c])) next
    if(mean_expression_shuffle[p, c] > 1){
      pval <- sum(shuffle_results[, c] > mean_expression_shuffle[p, c]) / shuffle_times
    } else {
      pval <- sum(shuffle_results[, c] < mean_expression_shuffle[p, c]) / shuffle_times
    }
    if(pval > 0.05) mean_expression_shuffle[p, c] <- NA 
    pvalues_mat[p, c] <- pval
  }
}

# remove all NA rows
valid_rows <- rowSums(!is.na(mean_expression_shuffle)) > 0
mean_expression_shuffle <- mean_expression_shuffle[valid_rows, ]

# save result
cat("save result")
write.table(mean_expression_noshuffle, "ProC_Pathway_Activity_Raw.txt", sep="\t")
write.table(mean_expression_shuffle, "ProC_Pathway_Activity_Filtered.txt", sep="\t")
write.table(pvalues_mat, "ProC_Pathway_Activity_pvalues.txt", sep="\t")








