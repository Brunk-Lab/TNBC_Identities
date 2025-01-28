.libPaths()
#Set working directory
setwd("/home/yue1118/TNBC_scRNA_analysis_workflow/")
#Libraries
library(Seurat)
library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)
library(reshape2)
library(stringr)
library(clusterProfiler)

#Building the Seurat object for HCC1143

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


#HCC1143 Control 24hr 
#Cells in this condition were treated with DMSO for 24hr

HCC1143_control_24hr_count_data <- read.table("/home/yue1118/patient_scRNA_seq_analysis/GSM4131388_HCC1143_24h_Ctrl_UMI_count.txt", header = TRUE, row.names = 1, sep = "\t")
HCC1143_control_24hr <- CreateSeuratObject(counts = HCC1143_control_24hr_count_data, min.cells = 100, min.features = 100)
HCC1143_control_24hr = RenameCells(HCC1143_control_24hr, add.cell.id = "HCC1143_control_24hr")
HCC1143_control_24hr = PercentageFeatureSet(object = HCC1143_control_24hr, pattern = "^MT-", col.name = "percent.mt")
HCC1143_control_24hr = subset(HCC1143_control_24hr, subset = nCount_RNA > 1000 & nFeature_RNA > 500 & percent.mt < 10)
HCC1143_control_24hr$dataset <- "control_24hr"

#Pseudobulk for control 24hr condition:
#HCC1143_control_24hr <- SCTransform(HCC1143_control_24hr, vars.to.regress = "percent.mt", verbose = FALSE, return.only.var.genes=FALSE, vst.flavor="v2")
#HCC1143_control_24hr_pseudo_bulk <- as.data.frame(AggregateExpression(HCC1143_control_24hr, assays = "SCT",return.seurat = F))
#colnames(HCC1143_control_24hr_pseudo_bulk) <- c('normalized_counts_HCC1143_control_24hr')
#write.csv(HCC1143_control_24hr_pseudo_bulk, './GSE139129_HCC1143_control_24hr_pseudo_bulk_03312024.csv')


#HCC1143 Treatment 24hr
#Cells in this condition were treated with Paclitaxel for 24hr

HCC1143_treated_24hr_count_data <- read.table("/home/yue1118/patient_scRNA_seq_analysis/GSM4131389_HCC1143_24h_Treated_UMI_count.txt", header = TRUE, row.names = 1, sep = "\t")
HCC1143_treated_24hr <- CreateSeuratObject(counts = HCC1143_treated_24hr_count_data, min.cells = 100, min.features = 100)
HCC1143_treated_24hr = RenameCells(HCC1143_treated_24hr, add.cell.id = "HCC1143_treated_24hr")
HCC1143_treated_24hr = PercentageFeatureSet(object = HCC1143_treated_24hr, pattern = "^MT-", col.name = "percent.mt")
HCC1143_treated_24hr = subset(HCC1143_treated_24hr, subset = nCount_RNA > 1000 & nFeature_RNA > 500 & percent.mt < 10)
HCC1143_treated_24hr$dataset <- "treated_24hr"


#HCC1143 Control 72hr 
#Cells in this condition were treated with DMSO for 72hr

HCC1143_control_72hr_count_data <- read.table("/home/yue1118/patient_scRNA_seq_analysis/GSM4131390_HCC1143_72h_Ctrl_UMI_count.txt", header = TRUE, row.names = 1, sep = "\t")
HCC1143_control_72hr <- CreateSeuratObject(counts = HCC1143_control_72hr_count_data, min.cells = 100, min.features = 100)
HCC1143_control_72hr = RenameCells(HCC1143_control_72hr, add.cell.id = "HCC1143_control_72hr")
HCC1143_control_72hr = PercentageFeatureSet(object = HCC1143_control_72hr, pattern = "^MT-", col.name = "percent.mt")
HCC1143_control_72hr = subset(HCC1143_control_72hr, subset = nCount_RNA > 1000 & nFeature_RNA > 500 & percent.mt < 10)
HCC1143_control_72hr$dataset <- "control_72hr"


#HCC1143 Treatment 72hr
#Cells in this condition were treated with Paclitaxel for 72hr

HCC1143_treated_72hr_count_data <- read.table("/home/yue1118/patient_scRNA_seq_analysis/GSM4131391_HCC1143_72h_Treated_UMI_count.txt", header = TRUE, row.names = 1, sep = "\t")
HCC1143_treated_72hr <- CreateSeuratObject(counts = HCC1143_treated_72hr_count_data, min.cells = 100, min.features = 100)
HCC1143_treated_72hr = RenameCells(HCC1143_treated_72hr, add.cell.id = "HCC1143_treated_72hr")
HCC1143_treated_72hr = PercentageFeatureSet(object = HCC1143_treated_72hr, pattern = "^MT-", col.name = "percent.mt")
HCC1143_treated_72hr = subset(HCC1143_treated_72hr, subset = nCount_RNA > 1000 & nFeature_RNA > 500 & percent.mt < 10)
HCC1143_treated_72hr$dataset <- "treated_72hr"


HCC1143.list = as.list(c(HCC1143_control_24hr, HCC1143_treated_24hr, HCC1143_control_72hr, HCC1143_treated_72hr))
names(HCC1143.list) = c("HCC1143_control_24hr", "HCC1143_treated_24hr", "HCC1143_control_72hr", "HCC1143_treated_72hr")


HCC1143.merged <- HCC1143.list[[1]]
for(i in 2:length(HCC1143.list)){
  HCC1143.merged <- merge(HCC1143.merged,HCC1143.list[[i]])
}

HCC1143.merged <- SCTransform(HCC1143.merged, vars.to.regress = "percent.mt", verbose = FALSE, return.only.var.genes=FALSE, vst.flavor="v2")

HCC1143.merged <- RunPCA(HCC1143.merged, npcs = 50, verbose = FALSE)
HCC1143.merged <- FindNeighbors(HCC1143.merged, reduction = "pca", dims = 1:50, verbose = FALSE)
HCC1143.merged <- FindClusters(HCC1143.merged, resolution = 0.7)
HCC1143.merged <- RunUMAP(HCC1143.merged, dims = 1:50, reduction = "pca", verbose = FALSE, seed.use = 7)
Idents(HCC1143.merged) <- HCC1143.merged$dataset

DefaultAssay(HCC1143.merged) <- "SCT"
#HCC1143_2 <- subset(HCC1143.merged, idents = c("control_24hr","treated_24hr"))

HCC1143.merged <- AddModuleScore(HCC1143.merged, features = list_estrogen, pool = NULL,nbin = 24,ctrl = 5,k = FALSE,assay = NULL,name = "estrogen_score",seed = 1,search = FALSE,slot = "data")
HCC1143.merged <- AddModuleScore(HCC1143.merged, features = list_NFKB, pool = NULL,nbin = 24,ctrl = 5,k = FALSE,assay = NULL,name = "NFKB_score",seed = 1,search = FALSE,slot = "data")
HCC1143.merged <- AddModuleScore(HCC1143.merged, features = list_EMT, pool = NULL,nbin = 24,ctrl = 5,k = FALSE,assay = NULL,name = "EMT_score",seed = 1,search = FALSE,slot = "data")
HCC1143.merged$condition <- sapply(str_split(colnames(HCC1143.merged), "_[ATCG]"), `[`, 1)
#saveRDS(HCC1143.merged, "/home/yue1118/patient_scRNA_seq_analysis/HCC1143_chemo_control_24hr_72hr_merged_10292024.rds")


#For Figure 3B & 3C
w <- read.delim("/home/yue1118/TNBC_paper_data/CCLE_TNBC_model/NMF_using_524_genes_selected_using_pan_cancer_cell_lines_including_MYC_on_30_updated_TNBC_cell_lines_normalized_for_samples_12222022/7/w.tsv")
w_matrix_genes <- w$X
empty_df <- data.frame(row.names = w$X)
empty_df$ID <- w$X
w_numeric <- as.matrix(sapply(w[, -1], as.numeric))

Idents(HCC1143.merged) <- "dataset"
for (i in c("control_24hr","treated_24hr","control_72hr", "treated_72hr")) {
  df_state <- data.frame(matrix(ncol = 7, nrow = 0))
  HCC1143_individual <- subset(HCC1143.merged, idents = i)
  sct_data <- GetAssayData(HCC1143_individual, assay = "SCT", slot = "counts")
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
  rownames(H_matrix) <- colnames(HCC1143_individual)
  df_state <- rbind(df_state, H_matrix)
  max_col_index <- apply(df_state, 1, which.max)
  max_col_names <- colnames(df_state)[max_col_index]
  df_max_state <- data.frame(rownames(df_state), max_col_names)
  df1 <- as.data.frame(table(df_max_state$max_col_names))
  df_state_percentage <- data.frame(category = df1$Var1, count = df1$Freq)
  df_state_percentage$percentage <- (df_state_percentage$count / sum(df_state_percentage$count)) *100
  ggplot(df_state_percentage, aes(x = "", y = count, fill = category)) +
    geom_bar(width = 1, stat = "identity", color = "white") +
    coord_polar("y", start = 0) +
    theme_void() +
    scale_fill_manual(values = c("V1" = "lightblue", "V2" = "blue",  "V3" = "chocolate1", "V4" = "darkred","V5"="brown2", "V6" = "darksalmon",  "V7"="aquamarine"))# + 
  #labs(title = paste("single-cell state map"), fill = "Category")
  ggsave(paste0("/home/yue1118/TNBC_scRNA_analysis_workflow/TNBC_paper_data_and_figure/GSE139129_HCC1143_", i, "_cancer_cells_pie_chart_11142024.svg"),dpi = 300, height = 5, width = 5)
}


#For Figure 3D
Idents(HCC1143.merged) <- "dataset"
df_state <- data.frame(matrix(ncol = 7, nrow = 0))
for (i in c("control_24hr","treated_24hr","control_72hr", "treated_72hr")) {
  HCC1143_individual <- subset(HCC1143.merged, idents = i)
  sct_data <- GetAssayData(HCC1143_individual, assay = "SCT", slot = "counts")
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
  rownames(H_matrix) <- colnames(HCC1143_individual)
  df_state <- rbind(df_state, H_matrix)
}
write.csv(df_state, '/home/yue1118/TNBC_scRNA_analysis_workflow/TNBC_paper_data_and_figure/HCC1143_single_cell_states_by_NMF_04182024.csv')

df_state <- read.csv('/home/yue1118/TNBC_scRNA_analysis_workflow/TNBC_paper_data_and_figure/HCC1143_single_cell_states_by_NMF_04182024.csv')
df_state$cell <- df_state$X
df_state_long <- melt(df_state, variable.name = "Variable", value.name = "Value")
df_state_long$condition <- sapply(str_split(df_state_long$cell, "hr_"), `[`, 1)
df_state_long$condition <- paste0(df_state_long$condition, "hr")
df_state_long_24hr <- df_state_long[grepl("_24hr", df_state_long$condition),]
df_state_long_72hr <- df_state_long[grepl("_72hr", df_state_long$condition),]
ggplot(df_state_long_24hr, aes(x = Variable, y = Value,fill = condition)) + 
  geom_boxplot() +
  labs(title = "Boxplot of Values by Group",
       x = "Group",
       y = "Value") +
  theme_minimal()
ggsave(paste0("./TNBC_paper_data_and_figure/GSE139129_HCC1143_24hr_boxplot_04242024.svg"), height = 3, width =10)

df_state_long_control_24hr <- df_state_long_24hr[df_state_long_24hr$condition == "HCC1143_control_24hr",]
df_state_long_treated_24hr <- df_state_long_24hr[df_state_long_24hr$condition == "HCC1143_treated_24hr",]
View(df_state_long_control_24hr)

wilcox.test(df_state_long_control_24hr[df_state_long_control_24hr$Variable == "V1",]$Value, df_state_long_treated_24hr[df_state_long_treated_24hr$Variable == "V1",]$Value) #p-value = 3.644e-05
wilcox.test(df_state_long_control_24hr[df_state_long_control_24hr$Variable == "V2",]$Value, df_state_long_treated_24hr[df_state_long_treated_24hr$Variable == "V2",]$Value) #p-value = 0.2898
wilcox.test(df_state_long_control_24hr[df_state_long_control_24hr$Variable == "V3",]$Value, df_state_long_treated_24hr[df_state_long_treated_24hr$Variable == "V3",]$Value) #p-value < 2.2e-16
wilcox.test(df_state_long_control_24hr[df_state_long_control_24hr$Variable == "V4",]$Value, df_state_long_treated_24hr[df_state_long_treated_24hr$Variable == "V4",]$Value) #p-value = 1.402e-07
wilcox.test(df_state_long_control_24hr[df_state_long_control_24hr$Variable == "V5",]$Value, df_state_long_treated_24hr[df_state_long_treated_24hr$Variable == "V5",]$Value) #p-value = 0.7472
wilcox.test(df_state_long_control_24hr[df_state_long_control_24hr$Variable == "V6",]$Value, df_state_long_treated_24hr[df_state_long_treated_24hr$Variable == "V6",]$Value) #p-value = 0.002788
wilcox.test(df_state_long_control_24hr[df_state_long_control_24hr$Variable == "V7",]$Value, df_state_long_treated_24hr[df_state_long_treated_24hr$Variable == "V7",]$Value) #p-value = 0.01786



#For Figure 3E
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
df_w_zscore$gene <- rownames(df_w_zscore)
custom_gene_set <- df_w_zscore[, c("max_attribute","gene")]
colnames(custom_gene_set) <- c("term", "gene")
rownames(custom_gene_set) <- NULL


df_deg <- read.csv('/home/yue1118/patient_scRNA_seq_analysis/HCC1143_cancer_cell_DEG_before_and_after_treatment_04092024.csv')
all_deg_gene_vector <- setNames(df_deg$avg_log2FC, df_deg$X)
all_deg_gene_vector <- sort(all_deg_gene_vector, decreasing = TRUE)

set.seed(7) # It is important to set seed to make sure you have reproductivity!!!
gsea_result_HCC1143 <- GSEA(
  geneList = all_deg_gene_vector,         # Ranked gene list
  TERM2GENE = custom_gene_set,  # Custom gene set
  pvalueCutoff = 0.05,          # Adjust p-value cutoff for significance
  verbose = FALSE,        # Use "ENTREZID" if using Entrez IDs
)

print(gsea_result_HCC1143)

gseaplot2(gsea_result_HCC1143, 
          geneSetID = "F0", 
          title = "", 
          # color = "red", 
          base_size = 14, 
          rel_heights = c(1.5, 0.5, 1))
ggsave("/home/yue1118/TNBC_scRNA_analysis_workflow/TNBC_paper_data_and_figure/HCC1143_chemo_24hr_F0_GSEA_enrichment_plot.svg", height = 5, width =8,dpi = 300)
gseaplot2(gsea_result_HCC1143, 
          geneSetID = "F1", 
          title = "", 
          # color = "red", 
          base_size = 14, 
          rel_heights = c(1.5, 0.5, 1))
ggsave("/home/yue1118/TNBC_scRNA_analysis_workflow/TNBC_paper_data_and_figure/HCC1143_chemo_24hr_F1_GSEA_enrichment_plot.svg", height = 5, width =8,dpi = 300)
gseaplot2(gsea_result_HCC1143, 
          geneSetID = "F2", 
          title = "", 
          # color = "red", 
          base_size = 14, 
          rel_heights = c(1.5, 0.5, 1))
ggsave("/home/yue1118/TNBC_scRNA_analysis_workflow/TNBC_paper_data_and_figure/HCC1143_chemo_24hr_F2_GSEA_enrichment_plot.svg", height = 5, width =8,dpi = 300)



#For Figure 3F
HCC1143.merged <- readRDS('/home/yue1118/patient_scRNA_seq_analysis/HCC1143_chemo_control_24hr_72hr_merged_10292024.rds')
df_plot <- HCC1143.merged@meta.data[,c("condition", "estrogen_score1", "NFKB_score1", "EMT_score1")]
df_plot$cells <- rownames(df_plot)
df_plot2 <- df_plot[, c(1,2,3,4)]
data_long <- melt(df_plot2, id.vars = "condition", variable.name = "pathway", value.name = "scores")
data_long$condition <- factor(data_long$condition, levels = c("HCC1143_control_24hr", "HCC1143_treated_24hr", "HCC1143_control_72hr", "HCC1143_treated_72hr"))
View(data_long)
data_long <- data_long[data_long$condition %in% c("HCC1143_control_24hr", "HCC1143_treated_24hr"),]
ggplot(data_long, aes(x = pathway, y = scores, fill = condition)) +
  geom_violin() +
  labs(title = "Expression Levels of Three Pathways", x = "Pathway", y = "pathway score level") +
  theme_minimal()
ggsave("/home/yue1118/TNBC_scRNA_analysis_workflow/TNBC_paper_data_and_figure/HCC143_chemo_24hr_MsigDB_pathways_violin_plots.svg", height = 5, width =10,dpi = 300)


Idents(HCC1143.merged) <- HCC1143.merged$dataset
HCC1143_pre <- subset(HCC1143.merged, idents = "control_24hr")
HCC1143_on <- subset(HCC1143.merged, idents = "treated_24hr")
wilcox.test(HCC1143_pre$estrogen_score1, HCC1143_on$estrogen_score1, paired = FALSE)
wilcox.test(HCC1143_pre$NFKB_score1, HCC1143_on$NFKB_score1, paired = FALSE)
wilcox.test(HCC1143_pre$EMT_score1, HCC1143_on$EMT_score1, paired = FALSE)

#For Figure 3G
w <- read.delim("/home/yue1118/TNBC_paper_data/CCLE_TNBC_model/NMF_using_524_genes_selected_using_pan_cancer_cell_lines_including_MYC_on_30_updated_TNBC_cell_lines_normalized_for_samples_12222022/7/w.tsv")
w_matrix_genes <- w$X
w_2 <- w
rownames(w_2) <- w_2$X
w_2 <- w_2[, c("F0", "F1", "F2", "F3", "F4", "F5", "F6")]
df_w_zscore <- as.data.frame(t(apply(w_2, 1, function(x) scale(x, center = TRUE))))
colnames(df_w_zscore) <- c("F0", "F1", "F2", "F3", "F4", "F5", "F6")
df_w_zscore <- df_w_zscore[,c("F0", "F1", "F2", "F3", "F4", "F5")] #Here we get rid of F6 state in that a lot of genes in F2 are overlapping with F6.

df_w_zscore$max_attribute <- apply(df_w_zscore, 1, function(x) names(df_w_zscore)[which.max(x)])
table(df_w_zscore$max_attribute)

#To create a function to make pie charts for upregulated genes
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

#To create a function to make pie charts for downregulated genes
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

#To create a function to make volcano plots for differentially expressed genes after treatment
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

#To compute differentially expressed genes after paclitaxel treatment
Idents(HCC1143.merged) <- HCC1143.merged$dataset
HCC1143.merged <- PrepSCTFindMarkers(HCC1143.merged)
DEG_HCC1143_cancer_cell <- FindMarkers(HCC1143.merged, ident.1 = "treated_24hr", ident.2 = "control_24hr", min.pct = 0.25)
write.csv(DEG_HCC1143_cancer_cell, '../patient_scRNA_seq_analysis/HCC1143_cancer_cell_DEG_before_and_after_treatment_04092024.csv')

HCC1143_DEG <- read.csv('/home/yue1118/patient_scRNA_seq_analysis/HCC1143_cancer_cell_DEG_before_and_after_treatment_04092024.csv') #treated_24hr vs control_24hr
df_HCC1143 <- HCC1143_DEG[HCC1143_DEG$p_val_adj < 0.05,]
df_HCC1143_2 <- HCC1143_DEG[(HCC1143_DEG$p_val_adj < 0.05) & (HCC1143_DEG$X %in% w_matrix_genes),]

create_pie_increased_DEG_from_table(df_HCC1143_2)
ggsave("/home/yue1118/TNBC_scRNA_analysis_workflow/TNBC_paper_data_and_figure/HCC1143_paclitaxel_24hr_pie_chart_increased_DEG_11142024.svg", height = 5, width =5)
create_pie_decreased_DEG_from_table(df_HCC1143_2)
ggsave("/home/yue1118/TNBC_scRNA_analysis_workflow/TNBC_paper_data_and_figure/HCC1143_paclitaxel_24hr_pie_chart_decreased_DEG_11142024.svg", height = 5, width =5)
volcano_plot_function(HCC1143_DEG)
ggsave("/home/yue1118/TNBC_scRNA_analysis_workflow/TNBC_paper_data_and_figure/HCC1143_paclitaxel_24hr_DEG_volcano_plot_10292024.svg", height = 5, width =6)


#To calculate identity change residual score for Fig 3I
df_state <- read.csv('/home/yue1118/TNBC_scRNA_analysis_workflow/TNBC_paper_data_and_figure/HCC1143_single_cell_states_by_NMF_04182024.csv')
df_state$cell <- df_state$X
df_state_control <- subset(df_state, grepl("HCC1143_control_24hr", cell))
df_state_treated <- subset(df_state, grepl("HCC1143_treated_24hr", cell))

df_control_non_zero_rows <- df_state_control[apply(df_state_control, 1, function(x) any(x != 0)), ]
df_treated_non_zero_rows <- df_state_treated[apply(df_state_treated, 1, function(x) any(x != 0)), ]

max_col_index <- apply(df_control_non_zero_rows, 1, which.max)
max_col_names <- colnames(df_control_non_zero_rows)[max_col_index]
df_control_max_state <- data.frame(rownames(df_control_non_zero_rows), max_col_names)
df1 <- as.data.frame(table(df_control_max_state$max_col_names))
df_state_percentage <- data.frame(category = df1$Var1, count = df1$Freq)
df_state_percentage$percentage <- (df_state_percentage$count / sum(df_state_percentage$count)) *100
df_state_percentage$condition <- "Control"
View(df_state_percentage)

max_col_index <- apply(df_treated_non_zero_rows, 1, which.max)
max_col_names <- colnames(df_treated_non_zero_rows)[max_col_index]
df_treated_max_state <- data.frame(rownames(df_treated_non_zero_rows), max_col_names)
df1 <- as.data.frame(table(df_treated_max_state$max_col_names))
df_state_percentage_2 <- data.frame(category = df1$Var1, count = df1$Freq)
df_state_percentage_2$percentage <- (df_state_percentage_2$count / sum(df_state_percentage_2$count)) *100
df_state_percentage_2$condition <- "Treated"


df_state_percentage
df_state_percentage_2

data <- matrix(c(16, 29, 10, 903,   # counts for "before" condition
                 25, 44, 30, 583),  # counts for "after" condition
               nrow = 4, ncol = 2,
               byrow = FALSE)
rownames(data) <- c("F0", "F1", "F2", "F3")
colnames(data) <- c("Before", "After")


