library(Seurat)
library(ggplot2)
library(openxlsx)
library(stringr)
library(dplyr)
library(ggridges)
library(ComplexHeatmap)


seurat <- readRDS("annotated_data.rds")
cell.levels <- c("InfC","preInfC","FC","preFC","HTC","preHTC","HomC","RepC","RegC","EC","ProC")
cell.color <- c("#717200","#763931","#2561A0","#F8C7B4","#DF5406","#2EB9CC","#C657B3","#ADD187",
                "#159931","#CB2B1F","#E99407")
cell.color <- setNames(cell.color,cell.levels)
group.levels <- c("Control","OA","KBD")
group.color <- c("#469393","#F1A93B","#C24976")
group.color <- setNames(group.color,group.levels)

gene_list <- readRDS("gene_list.rds") # pathway genes for apoptosis, ferroptosis, necroptosis, autophagy, pyroptosis, ER stress, oxidative stress, DNA damage, response to H2O2, iron accumulation, NLRP3 inflammasome, inflammatory response

##### AddModuleScore for target gene sets
seurat <- AddModuleScore(
  object = seurat,
  features = gene_list,
  assay = "RNA",
  slot = "data",
  name = "PathwayScore",
  ctrl = 100
)

##### ridge plot by group_celltype
score_columns <- colnames(seurat@meta.data)[c(10:21)]
for (score in score_columns) {
  data <- FetchData(seurat, vars = c("celltype", "group", score))
  data$celltype <- factor(data$celltype, levels = cell.levels)
  data$group <- factor(data$group, levels = group.levels)
  
  stat_lines <- data %>%
    group_by(celltype, group) %>%
    summarise(stat_val = mean(get(score), na.rm = TRUE), .groups = "drop") 
  
  stat_lines <- stat_lines %>%
    mutate(offset = case_when(
      group == "Control" ~  0,#adjust
      group == "OA"      ~  0,#adjust
      group == "KBD"     ~  0 #adjust
    ))
  
  p <- ggplot(data, aes(x = get(score), y = celltype, fill = group)) +
    geom_density_ridges(alpha = 0.8, scale = 0.9, rel_min_height = 0.01,
                        color = "black", size = 0.3) +
    geom_point(data = stat_lines,
               aes(x = stat_val, y = celltype, fill = group), 
               color = "black",  
               shape = 21, size = 3, stroke = 1.2,  
               position = position_nudge(x = stat_lines$offset)) +
    scale_fill_manual(values = group.color) +
    scale_color_manual(values = group.color) +
    labs(x = score, y = NULL) +
    theme_minimal(base_size = 14) +
    theme(
      legend.position = "right",
      axis.text = element_text(size = 12, color = "black"),  
      axis.ticks = element_line(color = "black"),  
      axis.title.x = element_text(size = 14, color = "black"),  
      axis.title.y = element_text(size = 14, color = "black"),  
      axis.line = element_line(color = "black", size = 0.5),  
      panel.grid = element_blank()
    ) +
    scale_x_continuous(limits = c(-0.125, 0.55))  
  
  pdf(paste0("group_celltype.", score, ".ridgeplot.pdf"), width = 7, height = 8)
  print(p)
  dev.off()
}

#### calculate earth mover's distance (EMD) between two conditions
library(tidyverse)
library(emdist)
library(transport)# R4.3

score_col <- "Apoptosis.score1"# take apotpsis as an example
data <- FetchData(seurat, vars = c("celltype", "group", score_col))

compute_emd_w1d <- function(x, y, p = 1) {
  x <- x[is.finite(x)]
  y <- y[is.finite(y)]
  transport::wasserstein1d(x, y, p = p)  # x and y as vector
}

emd_results <- data %>%
  group_split(celltype) %>%
  map_dfr(function(df_ct){
    gs <- sort(unique(df_ct$group))
    combn(gs, 2, simplify = FALSE) %>%
      map_dfr(function(gp){
        x <- df_ct %>% filter(group == gp[1]) %>% pull(!!sym(score_col))
        y <- df_ct %>% filter(group == gp[2]) %>% pull(!!sym(score_col))
        tibble(
          celltype = unique(df_ct$celltype),
          group1 = gp[1], group2 = gp[2],
          EMD_w1 = compute_emd_w1d(x, y, p = 1)
        )
      })
  })

as.data.frame(emd_results)
res <- as.data.frame(emd_results)
res <- res[order(res$EMD_w1,decreasing = T),]

#### correlation between cell death and inducers, take apoptosis as an example
score_columns <- colnames(seurat@meta.data)[c(14:21)] # apoptosis and inducers
correlation_matrix <- matrix(NA, nrow = length(score_columns) - 1, ncol = nrow(group_celltype_combinations))
rownames(correlation_matrix) <- score_columns[-1]  
colnames(correlation_matrix) <- apply(group_celltype_combinations, 1, function(x) paste(x, collapse = "_"))

for (i in 1:nrow(group_celltype_combinations)) {
  group <- group_celltype_combinations[i, "group"]
  celltype <- group_celltype_combinations[i, "celltype"]
  
  idx <- which(all$group == group & all$celltype == celltype)
  sub_seurat <- all[, idx]
  
  data <- FetchData(sub_seurat, vars = score_columns)
  
  for (j in 2:length(score_columns)) {  
    cor_test <- cor.test(data$Apoptosis.score1, data[[score_columns[j]]], method = "spearman")
    correlation_matrix[j - 1, i] <- cor_test$estimate
  }
}

# add annotation
order_levels <-c("Control_ProC","OA_ProC","KBD_ProC","Control_EC","OA_EC","KBD_EC","Control_RegC","OA_RegC","KBD_RegC",
                 "Control_RepC","OA_RepC","KBD_RepC","Control_HomC","OA_HomC","KBD_HomC","Control_preHTC","OA_preHTC","KBD_preHTC",
                 "Control_HTC","OA_HTC","KBD_HTC","Control_preFC","OA_preFC","KBD_preFC","OA_FC","KBD_FC",
                 "OA_preInfC","KBD_preInfC","OA_InfC","KBD_InfC")
col_info <- strsplit(order_levels, "_")
celltype <- sapply(col_info, `[`, 2)
group <- sapply(col_info, `[`, 1)

ha_column <- HeatmapAnnotation(
  Group = group,
  CellType = celltype,
  col = list(
    Group = group.color,
    CellType = cell.color
  ),
  annotation_height = unit(c(5, 5), "mm")
)


min_val <- min(correlation_matrix, na.rm = TRUE)
max_val <- max(correlation_matrix, na.rm = TRUE)
# scale data: -1 to 1
scaled_matrix <- 2 * ((correlation_matrix - min_val) / (max_val - min_val)) - 1
col_fun = colorRamp2(c(-1,0,1), c("#2c3385","white","firebrick3"))

ht <- Heatmap(
  as.matrix(scaled_matrix),
  name = "Apoptosis_correlation",
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  column_names_rot = 45,
  top_annotation = ha_column,
  column_split = factor(celltype, levels = rev(cell.levels)),
  gap = unit(1, "mm"),
  border_gp = gpar(col = "black", lty = 2, lwd = 1.5),
  row_names_gp = gpar(fontsize = 12),
  column_names_gp = gpar(fontsize = 10),
  border = TRUE,
  
  cell_fun = function(j, i, x, y, width, height, fill) {
    value <- correlation_matrix[i, j]
    if (!is.na(value) && abs(value) > 0.5) { 
      grid.text(
        sprintf("%.2f", value), 
        x, y,
        gp = gpar(fontsize = 8, col = ifelse(value > 0, "black", "white")) 
      )
    }
  }
)

pdf("Apoptosis_inducer.pdf", width = 12, height = 3.5)
ht
dev.off()
