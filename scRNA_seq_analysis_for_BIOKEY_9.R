##This R file is used for TNBC patient sample BIOKEY_9 analysis. In the codes, BIOKEY_9 is referred as P9.

.libPaths()
getwd()
#Set working directory
setwd("/home/yue1118/patient_scRNA_seq_analysis/")

library(Seurat)
library(dplyr)
library(reshape2)
library(stringr)
library(ggplot2)
library(nnls)
library(msigdbr)
library(clusterProfiler)
library(DoubletFinder)

w <- read.delim("/home/yue1118/TNBC_paper_data/CCLE_TNBC_model/NMF_using_524_genes_selected_using_pan_cancer_cell_lines_including_MYC_on_30_updated_TNBC_cell_lines_normalized_for_samples_12222022/7/w.tsv")
w_matrix_genes <- w$X
empty_df <- data.frame(row.names = w$X)
empty_df$ID <- w$X
w_numeric <- as.matrix(sapply(w[, -1], as.numeric))


metadata_cohort1 <- read.csv('./1872-BIOKEY_metaData_cohort1_web.csv')
rownames(metadata_cohort1) <- metadata_cohort1$Cell
View(metadata_cohort1)

table(metadata_cohort1$expansion)

table(metadata_cohort1[metadata_cohort1$BC_type =='TNBC',]$patient_id)

Belgium_cohort_1 <- readRDS('./1863-counts_cells_cohort1.rds')
Belgium_cohort_1 = CreateSeuratObject(Belgium_cohort_1, project="patient_TNBC_scRNA", min.cells = 100, min.features = 100)
Belgium_cohort_1 = PercentageFeatureSet(object = Belgium_cohort_1, pattern = "^MT-", col.name = "percent.mt")
Belgium_cohort_1 <- AddMetaData(Belgium_cohort_1, metadata = metadata_cohort1)
Belgium_cohort_1 = subset(Belgium_cohort_1, subset = nCount_RNA > 1000 & nFeature_RNA > 500 & percent.mt < 15 & nFeature_RNA < 6000)

#Subset patient no.9
Idents(Belgium_cohort_1) <- Belgium_cohort_1$patient_id
Belgium_cohort_1_P9 <- subset(Belgium_cohort_1, idents = c('BIOKEY_9'))
Belgium_cohort_1_P9 <- SCTransform(Belgium_cohort_1_P9, vars.to.regress = "percent.mt", verbose = FALSE, return.only.var.genes=FALSE, vst.flavor="v2")
Belgium_cohort_1_P9 <- RunPCA(Belgium_cohort_1_P9, npcs = 50, verbose = FALSE)
Belgium_cohort_1_P9 <- FindNeighbors(Belgium_cohort_1_P9, reduction = "pca", dims = 1:50, verbose = FALSE)
Belgium_cohort_1_P9 <- FindClusters(Belgium_cohort_1_P9, resolution = 0.7, verbose = FALSE)
Belgium_cohort_1_P9 <- RunUMAP(Belgium_cohort_1_P9, reduction = "pca", dims = 1:50, verbose = FALSE, seed.use = 7)
#saveRDS(Belgium_cohort_1_P9, file = "Belgium_cohort_1_P9_intra_tumoral_before_and_after_anti_PD1_10232024.rds")

#For Fig.4A
Belgium_cohort_1_P9 <- readRDS('/home/yue1118/patient_scRNA_seq_analysis/Belgium_cohort_1_P9_intra_tumoral_before_and_after_anti_PD1_10232024.rds')
DimPlot(Belgium_cohort_1_P9, group.by = "seurat_clusters", label = TRUE)
ggsave("/home/yue1118/TNBC_scRNA_analysis_workflow/TNBC_paper_data_and_figure/Belgium_P09_UMAP_seurat_clusters_11242024.svg", height = 5, width =5)
DimPlot(Belgium_cohort_1_P9, group.by = "timepoint", label = TRUE)
ggsave("/home/yue1118/TNBC_scRNA_analysis_workflow/TNBC_paper_data_and_figure/Belgium_P09_UMAP_timepoint_11242024.svg", height = 5, width =5)
DimPlot(Belgium_cohort_1_P9, group.by = "cellType", label = TRUE)
ggsave("/home/yue1118/TNBC_scRNA_analysis_workflow/TNBC_paper_data_and_figure/Belgium_P09_UMAP_cell_type_11242024.svg", height = 5, width =6)
Idents(Belgium_cohort_1_P9) <- Belgium_cohort_1_P9$seurat_clusters
DimPlot(Belgium_cohort_1_P9)

Idents(Belgium_cohort_1_P9) <- Belgium_cohort_1_P9$cellType
Belgium_cohort_1_P9_cancer_cell <- subset(Belgium_cohort_1_P9, idents = "Cancer_cell")

#For Fig.4B
Idents(Belgium_cohort_1_P9_cancer_cell) <- "timepoint"
df_state_all <- data.frame(matrix(ncol = 7, nrow = 0))
for (i in c("Pre", "On")) {
  df_state <- data.frame(matrix(ncol = 7, nrow = 0))
  P9_individual <- subset(Belgium_cohort_1_P9_cancer_cell, idents = i)
  sct_data <- GetAssayData(P9_individual, assay = "SCT", slot = "counts")
  common_genes <- intersect(rownames(sct_data), w_matrix_genes)
  sct_data_aligned <- sct_data[common_genes, ]
  sct_data_aligned <- as.data.frame(sct_data_aligned)
  sct_data_aligned$ID <- rownames(sct_data_aligned)
  sct_data_aligned_2 <- merge(empty_df, sct_data_aligned, by="ID", all.x = TRUE)
  rownames(sct_data_aligned_2) <- sct_data_aligned_2$ID
  sct_data_aligned_2 <- sct_data_aligned_2[,-1]
  sct_data_aligned_2[] <- lapply(sct_data_aligned_2, function(x) ifelse(is.na(x), 0, x))
  sct_data_aligned_2 <- sct_data_aligned_2[w_matrix_genes, ]
  sct_data_matrix <- as.matrix(sct_data_aligned_2)
  H_matrix <- matrix(0, nrow = ncol(w_numeric), ncol = ncol(sct_data_matrix))
  for (j in 1:ncol(sct_data_matrix)) {
    nnls_result <- nnls(w_numeric, sct_data_matrix[, j])
    H_matrix[, j] <- nnls_result$x
  }
  H_matrix <- t(H_matrix)
  H_matrix <- as.data.frame(H_matrix)
  rownames(H_matrix) <- colnames(P9_individual)
  df_state <- rbind(df_state, H_matrix)
  zero_rows <- df_state[apply(df_state, 1, function(x) all(x == 0)), ] #0 rows with all 0 values
  df_non_zero_rows <- df_state[apply(df_state, 1, function(x) any(x != 0)), ]
  max_col_index <- apply(df_non_zero_rows, 1, which.max)
  max_col_names <- colnames(df_non_zero_rows)[max_col_index]
  df_max_state <- data.frame(rownames(df_non_zero_rows), max_col_names)
  df1 <- as.data.frame(table(df_max_state$max_col_names))
  df_state_percentage <- data.frame(category = df1$Var1, count = df1$Freq)
  df_state_percentage$percentage <- (df_state_percentage$count / sum(df_state_percentage$count)) *100
  ggplot(df_state_percentage, aes(x = "", y = count, fill = category)) +
    geom_bar(width = 1, stat = "identity", color = "white") +
    coord_polar("y", start = 0) +
    theme_void() +
    scale_fill_manual(values = c("V1" = "lightblue", "V2" = "blue",  "V3" = "chocolate1", "V4" = "darkred","V5"="brown2", "V6" = "darksalmon",  "V7"="aquamarine")) #+ 
  #labs(title = paste("single-cell state map"), fill = "Category")
  ggsave(paste0("/home/yue1118/TNBC_scRNA_analysis_workflow/TNBC_paper_data_and_figure/Belgium_TNBC_P09_", i, "_cancer_cells_pie_chart_11142024.svg"),dpi = 300, height = 5, width = 5)
  df_state_all <- rbind(df_state_all, df_state)
}
write.csv(df_state_all, '/home/yue1118/TNBC_scRNA_analysis_workflow/TNBC_paper_data_and_figure/Belgium_P09_single_cell_states_by_NMF_04212024.csv')

#For Fig 4C.
df_state <- read.csv('/home/yue1118/TNBC_scRNA_analysis_workflow/TNBC_paper_data_and_figure/Belgium_P09_single_cell_states_by_NMF_04212024.csv')
View(df_state)
df_state$cell <- df_state$X
df_state_long <- melt(df_state, variable.name = "Variable", value.name = "Value")
View(df_state_long)
df_state_long$condition <- sapply(str_split(df_state_long$cell, "_[ATCG]"), `[`, 1)
df_state_long$condition <- factor(df_state_long$condition, levels = c("BIOKEY_9_Pre", "BIOKEY_9_On"))
ggplot(df_state_long, aes(x = Variable, y = Value,fill = condition)) + 
  geom_boxplot() +
  labs(title = "Boxplot of Values by Group",
       x = "Group",
       y = "Value") +
  theme_minimal()
ggsave(paste0("/home/yue1118/TNBC_scRNA_analysis_workflow/TNBC_paper_data_and_figure/Belgium_P09_boxplot_04252024.svg"), height = 5, width =10)

View(df_state_long)
df_state_long_pre <- df_state_long[df_state_long$condition == "BIOKEY_9_Pre",]
df_state_long_on <- df_state_long[df_state_long$condition == "BIOKEY_9_On",]

wilcox.test(df_state_long_pre[df_state_long_pre$Variable == "V1",]$Value, df_state_long_on[df_state_long_on$Variable == "V1",]$Value) #p-value = 5.416e-08
wilcox.test(df_state_long_pre[df_state_long_pre$Variable == "V2",]$Value, df_state_long_on[df_state_long_on$Variable == "V2",]$Value) #p-value < 2.2e-16
wilcox.test(df_state_long_pre[df_state_long_pre$Variable == "V3",]$Value, df_state_long_on[df_state_long_on$Variable == "V3",]$Value) #p-value < 2.2e-16
wilcox.test(df_state_long_pre[df_state_long_pre$Variable == "V4",]$Value, df_state_long_on[df_state_long_on$Variable == "V4",]$Value) #p-value = 1.435e-09
wilcox.test(df_state_long_pre[df_state_long_pre$Variable == "V5",]$Value, df_state_long_on[df_state_long_on$Variable == "V5",]$Value) #p-value = 0.04233
wilcox.test(df_state_long_pre[df_state_long_pre$Variable == "V6",]$Value, df_state_long_on[df_state_long_on$Variable == "V6",]$Value) #p-value = 0.6711
wilcox.test(df_state_long_pre[df_state_long_pre$Variable == "V7",]$Value, df_state_long_on[df_state_long_on$Variable == "V7",]$Value) #p-value = 3.756e-11


#For Fig 4D
w <- read.delim("/home/yue1118/TNBC_paper_data/CCLE_TNBC_model/NMF_using_524_genes_selected_using_pan_cancer_cell_lines_including_MYC_on_30_updated_TNBC_cell_lines_normalized_for_samples_12222022/7/w.tsv")
w_matrix_genes <- w$X
w_2 <- w
rownames(w_2) <- w_2$X
w_2 <- w_2[, c("F0", "F1", "F2", "F3", "F4", "F5", "F6")]
df_w_zscore <- as.data.frame(t(apply(w_2, 1, function(x) scale(x, center = TRUE))))
colnames(df_w_zscore) <- c("F0", "F1", "F2", "F3", "F4", "F5", "F6")
df_w_zscore <- df_w_zscore[,c("F0", "F1", "F2", "F3", "F4", "F5")]
View(df_w_zscore)

df_w_zscore$max_attribute <- apply(df_w_zscore, 1, function(x) names(df_w_zscore)[which.max(x)])
table(df_w_zscore$max_attribute)
View(df_w_zscore)
df_w_zscore$gene <- rownames(df_w_zscore)
custom_gene_set <- df_w_zscore[, c("max_attribute","gene")]
colnames(custom_gene_set) <- c("term", "gene")
rownames(custom_gene_set) <- NULL

df_deg <- read.csv('/home/yue1118/patient_scRNA_seq_analysis/Belgium_P9_cancer_cell_DEG_before_and_after_treatment_04092024.csv')
all_deg_gene_vector <- setNames(df_deg$avg_log2FC, df_deg$X)
all_deg_gene_vector <- sort(all_deg_gene_vector, decreasing = TRUE)

set.seed(7)
gsea_result <- GSEA(
  geneList = all_deg_gene_vector,         # Ranked gene list
  TERM2GENE = custom_gene_set,  # Custom gene set
  pvalueCutoff = 0.05,          # Adjust p-value cutoff for significance
  verbose = FALSE,        # Use "ENTREZID" if using Entrez IDs
)
print(gsea_result)
dotplot(gsea_result)

enrichplot::gseaplot2(gsea_result, geneSetID = "F2")
ggsave("/home/yue1118/TNBC_scRNA_analysis_workflow/TNBC_paper_data_and_figure/Belgium_P09_F2_GSEA_enrichment_plot.svg", height = 5, width =8)
enrichplot::gseaplot2(gsea_result, geneSetID = "F1")
ggsave("/home/yue1118/TNBC_scRNA_analysis_workflow/TNBC_paper_data_and_figure/Belgium_P09_F1_GSEA_enrichment_plot.svg", height = 5, width =8)


#For Fig 4E
estrogen_early_gene_sets <- read.gmt('/home/yue1118/TNBC_scRNA_analysis_workflow/TNBC_paper_data_and_figure/HALLMARK_ESTROGEN_RESPONSE_EARLY.gmt') #MsigDB Hallmark estrogen response early pathway genes
estrogen_early_gene_list <- estrogen_early_gene_sets$gene

estrogen_late_gene_sets <- read.gmt('/home/yue1118/TNBC_scRNA_analysis_workflow/TNBC_paper_data_and_figure/HALLMARK_ESTROGEN_RESPONSE_LATE.gmt') #MsigDB Hallmark estrogen response late pathway genes
estrogen_late_gene_list <- estrogen_late_gene_sets$gene

list_estrogen <- c(estrogen_early_gene_list, estrogen_late_gene_list)
list_estrogen <- unique(list_estrogen)
list_estrogen <- list(list_estrogen)

NFKB_gene_sets <- read.gmt('/home/yue1118/TNBC_scRNA_analysis_workflow/TNBC_paper_data_and_figure/HALLMARK_TNFA_SIGNALING_VIA_NFKB.gmt') #MsigDB Hallmark estrogen response late pathway genes
NFKB_gene_list <- NFKB_gene_sets$gene
list_NFKB <- list(NFKB_gene_list)

EMT_gene_sets <- read.gmt('/home/yue1118/TNBC_scRNA_analysis_workflow/TNBC_paper_data_and_figure/HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION.gmt')
EMT_gene_list <- EMT_gene_sets$gene
list_EMT <- list(EMT_gene_list)

Belgium_cohort_1_P9_cancer_cell <- AddModuleScore(Belgium_cohort_1_P9_cancer_cell, features = list_estrogen, pool = NULL,nbin = 24,ctrl = 5,k = FALSE,assay = NULL,name = "estrogen_score",seed = 1,search = FALSE,slot = "data")
Belgium_cohort_1_P9_cancer_cell <- AddModuleScore(Belgium_cohort_1_P9_cancer_cell, features = list_NFKB, pool = NULL,nbin = 24,ctrl = 5,k = FALSE,assay = NULL,name = "NFKB_score",seed = 1,search = FALSE,slot = "data")
Belgium_cohort_1_P9_cancer_cell <- AddModuleScore(Belgium_cohort_1_P9_cancer_cell, features = list_EMT, pool = NULL,nbin = 24,ctrl = 5,k = FALSE,assay = NULL,name = "EMT_score",seed = 1,search = FALSE,slot = "data")
df_plot <- Belgium_cohort_1_P9_cancer_cell@meta.data[,c("timepoint", "estrogen_score1", "NFKB_score1", "EMT_score1")]
df_plot$cells <- rownames(df_plot)
df_plot2 <- df_plot[, c(1,2,3,4)]
data_long <- melt(df_plot2, id.vars = "timepoint", variable.name = "pathway", value.name = "scores")
data_long$timepoint <- factor(data_long$timepoint, levels = c("Pre", "On"))
View(data_long)
ggplot(data_long, aes(x = pathway, y = scores, fill = timepoint)) +
  geom_violin() +
  labs(title = "Expression Levels of Three Pathways", x = "Pathway", y = "pathway score level") +
  theme_minimal() + 
  theme(
    panel.grid.major.x = element_blank(),               # Remove major grid lines on x-axis if needed
    panel.grid.major.y = element_blank()                # Remove major grid lines on y-axis if needed
  )
ggsave("/home/yue1118/TNBC_scRNA_analysis_workflow/TNBC_paper_data_and_figure/Belgium_TNBC_P09_MsigDB_pathway_scores_violin_plots_10282024.svg",dpi = 300, height = 5, width = 10)

#For Fig.4F
w <- read.delim("/home/yue1118/TNBC_paper_data/CCLE_TNBC_model/NMF_using_524_genes_selected_using_pan_cancer_cell_lines_including_MYC_on_30_updated_TNBC_cell_lines_normalized_for_samples_12222022/7/w.tsv")
w_matrix_genes <- w$X
w_2 <- w
rownames(w_2) <- w_2$X
w_2 <- w_2[, c("F0", "F1", "F2", "F3", "F4", "F5", "F6")]
df_w_zscore <- as.data.frame(t(apply(w_2, 1, function(x) scale(x, center = TRUE))))
colnames(df_w_zscore) <- c("F0", "F1", "F2", "F3", "F4", "F5", "F6")
df_w_zscore <- df_w_zscore[,c("F0", "F1", "F2", "F3", "F4", "F5")]

df_w_zscore$max_attribute <- apply(df_w_zscore, 1, function(x) names(df_w_zscore)[which.max(x)])
table(df_w_zscore$max_attribute)

create_pie_decreased_DEG_from_table <- function(input_df) {
  # Step 1: Create a frequency table
  freq_table <- table(df_w_zscore[rownames(df_w_zscore) %in% input_df[input_df$avg_log2FC < 0,]$X,]$max_attribute)
  
  # Step 2: Convert the frequency table to a dataframe
  df <- as.data.frame(freq_table)
  df_state_percentage <- data.frame(category = df$Var1, count = df$Freq)
  df_state_percentage$percentage <- (df_state_percentage$count / sum(df_state_percentage$count)) *100
  
  # Step 3: Create a pie chart from the dataframe
  ggplot(df_state_percentage, aes(x = "", y = count, fill = category)) +
    geom_bar(width = 1, stat = "identity", color = "white") +
    coord_polar("y", start = 0) +
    theme_void() +
    scale_fill_manual(values = c("F0" = "lightblue", "F1" = "blue",  "F2" = "chocolate1", "F3" = "darkred","F4"="brown2", "F5" = "darksalmon",  "F6"="aquamarine"))
}


create_pie_increased_DEG_from_table <- function(input_df) {
  # Step 1: Create a frequency table
  freq_table <- table(df_w_zscore[rownames(df_w_zscore) %in% input_df[input_df$avg_log2FC > 0,]$X,]$max_attribute)
  
  # Step 2: Convert the frequency table to a dataframe
  df <- as.data.frame(freq_table)
  df_state_percentage <- data.frame(category = df$Var1, count = df$Freq)
  df_state_percentage$percentage <- (df_state_percentage$count / sum(df_state_percentage$count)) *100
  
  # Step 3: Create a pie chart from the dataframe
  ggplot(df_state_percentage, aes(x = "", y = count, fill = category)) +
    geom_bar(width = 1, stat = "identity", color = "white") +
    coord_polar("y", start = 0) +
    theme_void() +
    scale_fill_manual(values = c("F0" = "lightblue", "F1" = "blue",  "F2" = "chocolate1", "F3" = "darkred","F4"="brown2", "F5" = "darksalmon",  "F6"="aquamarine"))
}


volcano_plot_function <- function(input_df) {
  # Convert p-value to -log10 p-value for the volcano plot
  input_df$neg_log10_pvalue <- -log10(input_df$p_val_adj)
  
  # Add a column to classify significance (adjust thresholds as needed)
  input_df$significant <- ifelse(input_df$avg_log2FC > 0 & input_df$p_val_adj < 0.05, "Upregulated",
                                 ifelse(input_df$avg_log2FC < 0 & input_df$p_val_adj < 0.05, "Downregulated", "Not Significant"))
  
  # Create the volcano plot
  ggplot(input_df, aes(x = avg_log2FC, y = neg_log10_pvalue)) +
    geom_point(aes(color = significant)) + # Color points by significance
    scale_color_manual(values = c("blue", "gray", "red")) + # Custom colors
    labs(title = "Volcano Plot", x = "Log2 Fold Change", y = "-log10(adj p-value)") +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") + # Add dashed lines for log2FC thresholds
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") + 
    theme(
      panel.grid.major.x = element_blank(),               # Remove major grid lines on x-axis if needed
      panel.grid.major.y = element_blank()                # Remove major grid lines on y-axis if needed
    )# Horizontal line for p-value threshold
}

P09_DEG <- read.csv('../patient_scRNA_seq_analysis/Belgium_P9_cancer_cell_DEG_before_and_after_treatment_04092024.csv')
df_P09_2 <- P09_DEG[(P09_DEG$p_val_adj < 0.05) & (P09_DEG$X %in% w_matrix_genes),]

create_pie_increased_DEG_from_table(df_P09_2)
ggsave("/home/yue1118/TNBC_scRNA_analysis_workflow/TNBC_paper_data_and_figure/Belgium_P09_pie_chart_increased_DEG_11142024.svg", height = 5, width =5)
create_pie_decreased_DEG_from_table(df_P09_2)
ggsave("/home/yue1118/TNBC_scRNA_analysis_workflow/TNBC_paper_data_and_figure/Belgium_P09_pie_chart_decreased_DEG_11142024.svg", height = 5, width =5)
volcano_plot_function(P09_DEG)
ggsave("/home/yue1118/TNBC_scRNA_analysis_workflow/TNBC_paper_data_and_figure/Belgium_P09_DEG_volcano_plot_10282024.svg", height = 5, width =6)



#For Fig.4G
for (i in c(0,3,6)) {
  for (h in c("Pre", "On")){
    Idents(Belgium_cohort_1_P9_cancer_cell) <- Belgium_cohort_1_P9_cancer_cell$timepoint
    P9_individual <- subset(Belgium_cohort_1_P9_cancer_cell, idents = h)
    Idents(P9_individual) <- P9_individual$seurat_clusters
    P9_individual <- subset(P9_individual, idents = i)
    print(c(i, h, dim(P9_individual)))
    df_state <- data.frame(matrix(ncol = 7, nrow = 0))
    sct_data <- GetAssayData(P9_individual, assay = "SCT", slot = "counts")
    common_genes <- intersect(rownames(sct_data), w_matrix_genes)
    sct_data_aligned <- sct_data[common_genes, ]
    sct_data_aligned <- as.data.frame(sct_data_aligned)
    sct_data_aligned$ID <- rownames(sct_data_aligned)
    sct_data_aligned_2 <- merge(empty_df, sct_data_aligned, by="ID", all.x = TRUE)
    rownames(sct_data_aligned_2) <- sct_data_aligned_2$ID
    sct_data_aligned_2 <- sct_data_aligned_2[,-1]
    sct_data_aligned_2[] <- lapply(sct_data_aligned_2, function(x) ifelse(is.na(x), 0, x))
    sct_data_aligned_2 <- sct_data_aligned_2[w_matrix_genes, ]
    sct_data_matrix <- as.matrix(sct_data_aligned_2)
    H_matrix <- matrix(0, nrow = ncol(w_numeric), ncol = ncol(sct_data_matrix))
    for (j in 1:ncol(sct_data_matrix)) {
      nnls_result <- nnls(w_numeric, sct_data_matrix[, j])
      H_matrix[, j] <- nnls_result$x
    }
    H_matrix <- t(H_matrix)
    H_matrix <- as.data.frame(H_matrix)
    rownames(H_matrix) <- colnames(P9_individual)
    df_state <- rbind(df_state, H_matrix)
    zero_rows <- df_state[apply(df_state, 1, function(x) all(x == 0)), ] #0 rows with all 0 values 
    df_non_zero_rows <- df_state[apply(df_state, 1, function(x) any(x != 0)), ]
    max_col_index <- apply(df_non_zero_rows, 1, which.max)
    max_col_names <- colnames(df_non_zero_rows)[max_col_index]
    df_max_state <- data.frame(rownames(df_non_zero_rows), max_col_names)
    df1 <- as.data.frame(table(df_max_state$max_col_names))
    df_state_percentage <- data.frame(category = df1$Var1, count = df1$Freq)
    df_state_percentage$percentage <- (df_state_percentage$count / sum(df_state_percentage$count)) *100
    ggplot(df_state_percentage, aes(x = "", y = count, fill = category)) +
      geom_bar(width = 1, stat = "identity", color = "white") +
      coord_polar("y", start = 0) +
      theme_void() +
      scale_fill_manual(values = c("V1" = "lightblue", "V2" = "blue",  "V3" = "chocolate1", "V4" = "darkred","V5"="brown2", "V6" = "darksalmon",  "V7"="aquamarine"))# + 
    #labs(title = paste("single-cell state map"), fill = "Category")
    ggsave(paste0("/home/yue1118/TNBC_scRNA_analysis_workflow/TNBC_paper_data_and_figure/Belgium_TNBC_P09_intra_tumoral_condition_",h,"_cluster_", i, "_cancer_cells_pie_chart_11142024.svg"),dpi = 300, height = 5, width = 5)
  }
}

#For the data for P9 in Fig.4H (data for other patients was generated in the same way)
df_state <- read.csv('/home/yue1118/TNBC_scRNA_analysis_workflow/TNBC_paper_data_and_figure/Belgium_P09_single_cell_states_by_NMF_04212024.csv')
df_state$cell <- df_state$X
df_state_pre <- subset(df_state, grepl("Pre", cell))
df_state_on <- subset(df_state, grepl("On", cell))

df_pre_non_zero_rows <- df_state_pre[apply(df_state_pre, 1, function(x) any(x != 0)), ]
df_on_non_zero_rows <- df_state_on[apply(df_state_on, 1, function(x) any(x != 0)), ]

max_col_index <- apply(df_pre_non_zero_rows, 1, which.max)
max_col_names <- colnames(df_pre_non_zero_rows)[max_col_index]
df_pre_max_state <- data.frame(rownames(df_pre_non_zero_rows), max_col_names)
df1 <- as.data.frame(table(df_pre_max_state$max_col_names))
df_state_percentage <- data.frame(category = df1$Var1, count = df1$Freq)
df_state_percentage$percentage <- (df_state_percentage$count / sum(df_state_percentage$count)) *100
df_state_percentage$condition <- "Pre"
View(df_state_percentage)

max_col_index <- apply(df_on_non_zero_rows, 1, which.max)
max_col_names <- colnames(df_on_non_zero_rows)[max_col_index]
df_on_max_state <- data.frame(rownames(df_on_non_zero_rows), max_col_names)
df1 <- as.data.frame(table(df_on_max_state$max_col_names))
df_state_percentage_2 <- data.frame(category = df1$Var1, count = df1$Freq)
df_state_percentage_2$percentage <- (df_state_percentage_2$count / sum(df_state_percentage_2$count)) *100
df_state_percentage_2$condition <- "On"


df_state_percentage
df_state_percentage_2

data <- matrix(c(21,  45,  39, 111,  54,  14,   # counts for "before" condition
                 4,   9, 124,  38,   2,   1),  # counts for "after" condition
               nrow = 6, ncol = 2,
               byrow = FALSE)
rownames(data) <- c("F0", "F1", "F2", "F3", "F5", "F6")
colnames(data) <- c("Before", "After")

data

result <- chisq.test(data)
print(result)

result$stdres


#For the data for P9 in Fig.4I (data for other patients was generated in the same way)
Idents(Belgium_cohort_1_P9) <- Belgium_cohort_1_P9$cellType
Belgium_cohort_1_P9_cancer_cell <- subset(Belgium_cohort_1_P9, idents = "Cancer_cell")

df_state_percentage_all <- data.frame(matrix(ncol = 4, nrow = 0))
for (i in c(0,3,6)) {
  for (h in c("Pre", "On")){
    Idents(Belgium_cohort_1_P9_cancer_cell) <- Belgium_cohort_1_P9_cancer_cell$timepoint
    P9_individual <- subset(Belgium_cohort_1_P9_cancer_cell, idents = h)
    Idents(P9_individual) <- P9_individual$seurat_clusters
    P9_individual <- subset(P9_individual, idents = i)
    print(c(i,h, dim(P9_individual)))
    df_state <- data.frame(matrix(ncol = 7, nrow = 0))
    sct_data <- GetAssayData(P9_individual, assay = "SCT", slot = "counts")
    common_genes <- intersect(rownames(sct_data), w_matrix_genes)
    sct_data_aligned <- sct_data[common_genes, ]
    sct_data_aligned <- as.data.frame(sct_data_aligned)
    sct_data_aligned$ID <- rownames(sct_data_aligned)
    sct_data_aligned_2 <- merge(empty_df, sct_data_aligned, by="ID", all.x = TRUE)
    rownames(sct_data_aligned_2) <- sct_data_aligned_2$ID
    sct_data_aligned_2 <- sct_data_aligned_2[,-1]
    sct_data_aligned_2[] <- lapply(sct_data_aligned_2, function(x) ifelse(is.na(x), 0, x))
    sct_data_aligned_2 <- sct_data_aligned_2[w_matrix_genes, ]
    sct_data_matrix <- as.matrix(sct_data_aligned_2)
    H_matrix <- matrix(0, nrow = ncol(w_numeric), ncol = ncol(sct_data_matrix))
    for (j in 1:ncol(sct_data_matrix)) {
      nnls_result <- nnls(w_numeric, sct_data_matrix[, j])
      H_matrix[, j] <- nnls_result$x
    }
    H_matrix <- t(H_matrix)
    H_matrix <- as.data.frame(H_matrix)
    rownames(H_matrix) <- colnames(P9_individual)
    df_state <- rbind(df_state, H_matrix)
    zero_rows <- df_state[apply(df_state, 1, function(x) all(x == 0)), ] #0 rows with all 0 values 
    df_non_zero_rows <- df_state[apply(df_state, 1, function(x) any(x != 0)), ]
    max_col_index <- apply(df_non_zero_rows, 1, which.max)
    max_col_names <- colnames(df_non_zero_rows)[max_col_index]
    df_max_state <- data.frame(rownames(df_non_zero_rows), max_col_names)
    df1 <- as.data.frame(table(df_max_state$max_col_names))
    df_state_percentage <- data.frame(category = df1$Var1, count = df1$Freq)
    df_state_percentage$percentage <- (df_state_percentage$count / sum(df_state_percentage$count)) *100
    df_state_percentage$condition <- h
    df_state_percentage$cluster_number <- i
    df_state_percentage_all <- rbind(df_state_percentage_all, df_state_percentage)
  }
}

write.csv(df_state_percentage_all, '/home/yue1118/TNBC_scRNA_analysis_workflow/TNBC_paper_data_and_figure/P09_cancer_cell_sub_clusters_state_assignment_11072024.csv')

View(df_state_percentage_all[(df_state_percentage_all$category=='V3') & (df_state_percentage_all$condition=='Pre'),])

View(df_state_percentage_all[(df_state_percentage_all$category=='V3') & (df_state_percentage_all$condition=='On'),])

df_state_percentage_all[(df_state_percentage_all$category=='V3') & (df_state_percentage_all$condition=='Pre'),]
df_state_percentage_all[(df_state_percentage_all$category=='V3') & (df_state_percentage_all$condition=='On'),]

data <- matrix(c(16,9,14,   # counts for "before" condition
                 1,15,108),  # counts for "after" condition
               nrow = 3, ncol = 2,
               byrow = FALSE)
rownames(data) <- c("cluster_0", "cluster_3", "cluster_6")
colnames(data) <- c("Before", "After")

data

result <- chisq.test(data)
print(result)

result$stdres




