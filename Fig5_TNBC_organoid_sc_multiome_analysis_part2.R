#Set working directory
.libPaths()
setwd("/home/yue1118/TNBC_scRNA_analysis_workflow/")

library(Signac)
library(Seurat)
library(GenomicRanges)
library(future)
library(EnsDb.Hsapiens.v86)
library(dplyr)
library(tidyr)
library(ggplot2)
library(chromVAR)
library(JASPAR2020)
library(TFBSTools)
library(motifmatchr)
library(BSgenome.Hsapiens.UCSC.hg38)
library(MACSr)
library(reshape2)
library(pheatmap)
library(clusterProfiler)
library(stringr)
library(nnls)
library(biomaRt)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)  # Change to correct genome (e.g., hg19 or mm10)
library(org.Hs.eg.db) 
library(svglite)
library(msigdbr)
library(ChIPseeker)


NFKB_gene_sets <- read.gmt('/home/yue1118/TNBC_scRNA_analysis_workflow/TNBC_paper_data_and_figure/HALLMARK_TNFA_SIGNALING_VIA_NFKB.gmt') #MsigDB Hallmark NFKB genes
NFKB_gene_list <- NFKB_gene_sets$gene
list_NFKB <- list(NFKB_gene_list)

#To load the organoid seurat object
MS177_organoid_3_conditions <- readRDS('/home/yue1118/TNBC_scRNA_analysis_workflow/Organoid_TNBC_JQ1_MS177_24hr_3_conditions_with_gene_activity_with_cluster_identities_with_doublets_classification_01082025.rds')
Idents(MS177_organoid_3_conditions) <- "DF.classification"
MS177_organoid_3_conditions <- subset(MS177_organoid_3_conditions, idents = "Singlet")

#For Fig 5I.

#To get the differential TF motifs after MS177 treatment
DefaultAssay(MS177_organoid_3_conditions) <- 'chromvar'
Idents(MS177_organoid_3_conditions) <- MS177_organoid_3_conditions$dataset
MS177_motif_markers <- FindMarkers(MS177_organoid_3_conditions, ident.1 = "MS177_24hr", ident.2 =  "control_24hr", min.pct = 0.25, logfc.threshold = 0.10, recorrect_umi=F, mean.fxn = rowMeans,fc.name = "avg_diff")

#To get all the upregulated genes after MS177
DefaultAssay(MS177_organoid_3_conditions) <- 'SCT'
DEG_organoid_MS177 <- FindMarkers(MS177_organoid_3_conditions, ident.1 = "MS177_24hr", ident.2 = "control_24hr", min.pct = 0.25)
DEG_organoid_MS177_all_up_regulated_genes <- DEG_organoid_MS177[(DEG_organoid_MS177$p_val_adj < 0.05) & (DEG_organoid_MS177$avg_log2FC > 0),]
write.csv(DEG_organoid_MS177_all_up_regulated_genes, '/home/yue1118/TNBC_scRNA_analysis_workflow/TNBC_paper_data_and_figure/all_upregulated_genes_after_MS177_in_the_organoid_doublet_removed_pval_adj_lower_than_0.05_and_avg_log2FC_more_than_0_01082025.csv')
#The genes included in the dataset above were input into EnrichR database to look for potentially involved TFs
#The ChEA 2022 database was used to extract all the potentially involved TFs
#The result file is named "ChEA_2022_table_using_all_the_upregulated_genes_after_MS177_doublet_removed_01082025.txt"


#To get all the upregulated genes in NFKB hallmark pathway after MS177
DEG_organoid_MS177_all_up_regulated_NFKB_genes <- DEG_organoid_MS177[(DEG_organoid_MS177$p_val_adj < 0.05) & (DEG_organoid_MS177$avg_log2FC > 0) & (rownames(DEG_organoid_MS177) %in% NFKB_gene_sets$gene),]
write.csv(DEG_organoid_MS177_all_up_regulated_NFKB_genes, '/home/yue1118/TNBC_scRNA_analysis_workflow/TNBC_paper_data_and_figure/all_upregulated_NFKB_genes_after_MS177_in_the_organoid_doublet_removed_pval_adj_lower_than_0.05_and_avg_log2FC_more_than_0_01082025.csv')
#The genes included in the dataset above were input into EnrichR database to look for potentially involved TFs
#The ChEA 2022 database was used to extract all the potentially involved TFs
#The result file is named "ChEA_2022_table_using_all_the_upregulated_NFKB_genes_after_MS177_doublet_removed_01082025.txt"


pfms <- getMatrixSet(JASPAR2020,opts = list(species = 9606, all_versions = FALSE))
i = "MA0105.4"
pfms@listData[[i]]@name

list_names <- c()
list_motifs <- rownames(MS177_organoid_3_conditions@assays$chromvar@data)
for (i in list_motifs) {
  list_names <- c(list_names, pfms@listData[[i]]@name)
}
length(list_names)
length(unique(list_names))
df_tf <- data.frame(Name = list_names, Motif = list_motifs)

df1 <- MS177_motif_markers[(MS177_motif_markers$p_val_adj < 0.05) & (MS177_motif_markers$avg_diff > 0.25),]
df1$Motif <- rownames(df1)

df2 <- merge(df1, df_tf, by= "Motif")
df2 <- df2[,c("Name", "avg_diff", "Motif", "p_val_adj")]
colnames(df2) <- c("TF", "avg_diff", "Motif", "motif_p_val_adj")
View(df2)

list_tf_with_motf_changes <- df2$TF
list_NFKB_motifs <- intersect(NFKB_gene_sets$gene, list_tf_with_motf_changes)

#For all the up-regulated genes after MS177 in general  #To get the df_chea dataset as shown below, we get all the up-regulated
df_chea <- read.table('../TNBC_scRNA_analysis_workflow/TNBC_paper_data_and_figure/ChEA_2022_table_using_all_the_upregulated_genes_after_MS177_doublet_removed_01082025.txt',header = TRUE, sep = "\t", fileEncoding = "UTF-8")
df_chea$TF <- sapply(str_split(df_chea$Term, " "), `[`, 1)
df_chea <- df_chea[df_chea$Adjusted.P.value < 0.05,]
df_chea <- df_chea[grep("Human", df_chea$Term),]
df_chea_2 <- df_chea[df_chea$TF %in% df2$TF,]
df_chea_2 <- df_chea_2[!duplicated(df_chea_2$TF), ]
df_chea_2 <- df_chea_2[, c("Term", "Adjusted.P.value", "TF")]
colnames(df_chea_2) <- c("global_DEG_TF_Term", "global_DEG_TF_enrichment_adjusted.P.value", "TF")
df3 <- merge(df2, df_chea_2, by= "TF", all.x=T)
View(df3)
View(df_chea)

df_chea_nfkb <- read.table('../TNBC_scRNA_analysis_workflow/TNBC_paper_data_and_figure/ChEA_2022_table_using_all_the_upregulated_NFKB_genes_after_MS177_doublet_removed_01082025.txt',header = TRUE, sep = "\t", fileEncoding = "UTF-8")
df_chea_nfkb$TF <- sapply(str_split(df_chea_nfkb$Term, " "), `[`, 1)
df_chea_nfkb <- df_chea_nfkb[df_chea_nfkb$Adjusted.P.value < 0.05,]
df_chea_nfkb <- df_chea_nfkb[grep("Human", df_chea_nfkb$Term),]
df_chea_nfkb_2 <- df_chea_nfkb[df_chea_nfkb$TF %in% df2$TF,]
df_chea_nfkb_2 <- df_chea_nfkb_2[!duplicated(df_chea_nfkb_2$TF), ]
df_chea_nfkb_2 <- df_chea_nfkb_2[, c("Term", "Adjusted.P.value", "TF")]
colnames(df_chea_nfkb_2) <- c("NFKB_DEG_TF_Term", "NFKB_DEG_TF_enrichment_adjusted.P.value", "TF")
df4 <- merge(df3, df_chea_nfkb_2, by= "TF", all.x=T)
View(df4)

df4 <- df4[!(is.na(df4$global_DEG_TF_Term) & is.na(df4$NFKB_DEG_TF_Term)), ]
df4$NFKB_or_not <- ifelse(df4$TF %in% list_NFKB_motifs, "Yes", "No")
df4$NFKB_or_not <- factor(df4$NFKB_or_not, levels = c("Yes", "No"))
df4 <- df4[, c("TF", "Motif","avg_diff", "global_DEG_TF_enrichment_adjusted.P.value", "NFKB_DEG_TF_enrichment_adjusted.P.value", "NFKB_or_not")]
df4 <- df4[order(df4$NFKB_or_not), ]
df4$global_DEG_TF_enrichment_adjusted.P.value <- -log(df4$global_DEG_TF_enrichment_adjusted.P.value)
df4$NFKB_DEG_TF_enrichment_adjusted.P.value <- -log(df4$NFKB_DEG_TF_enrichment_adjusted.P.value)
rownames(df4) <- df4$TF


View(df4)
color_palette <- colorRampPalette(colors = c("white", "red"))(100)
color_palette_2 <- colorRampPalette(colors = c("white", "orange"))(100)
svglite('/home/yue1118/TNBC_scRNA_analysis_workflow/TNBC_paper_data_and_figure/MS177_NFKB_enriched_TFs_doublet_removed_01082025.svg', width = 6, height = 12)
pheatmap(df4[, c("global_DEG_TF_enrichment_adjusted.P.value", "NFKB_DEG_TF_enrichment_adjusted.P.value")], cluster_rows = F, cluster_cols = F, color = color_palette)
dev.off()

svglite('/home/yue1118/TNBC_scRNA_analysis_workflow/TNBC_paper_data_and_figure/MS177_upregulated_TF_motifs_doublet_removed_01082025.svg', width = 6, height = 12)
pheatmap(as.data.frame(df4[, c("avg_diff")]), cluster_rows = F, cluster_cols = F, color = color_palette_2)
dev.off()

NFKB_gene_sets <- read.gmt('/home/yue1118/TNBC_scRNA_analysis_workflow/TNBC_paper_data_and_figure/HALLMARK_TNFA_SIGNALING_VIA_NFKB.gmt') #MsigDB Hallmark NFKB genes
NFKB_gene_list <- NFKB_gene_sets$gene
list_NFKB <- list(NFKB_gene_list)




#For Fig 5J.
df1 <- read.csv('/home/yue1118/TNBC_scRNA_analysis_workflow/TNBC_paper_data_and_figure/TF_binding_motifs_target_genes_increased_after_MS177_avg_difference_more_than_0.25.csv')
View(df1)
rela_gene_sets <-  unlist(df1[8, 4:ncol(df1)]) #ChIP-Seq FIBROSARCOMA Human  ChEA_2022
list_overlap_rela <- list(intersect(rela_gene_sets, NFKB_gene_list))
DefaultAssay(MS177_organoid_3_conditions) <- "SCT"
MS177_organoid_3_conditions <- AddModuleScore(MS177_organoid_3_conditions, features = list_overlap_rela, pool = NULL,nbin = 10,ctrl = 5,k = FALSE,assay = NULL,name = "rela_target_score",seed = 1,search = FALSE,slot = "data")

MS177_organoid_3_conditions_subset_1 <- subset(MS177_organoid_3_conditions, idents = c("control_24hr"))
MS177_organoid_3_conditions_subset_2 <- subset(MS177_organoid_3_conditions, idents = c("JQ1_24hr"))
MS177_organoid_3_conditions_subset_3 <- subset(MS177_organoid_3_conditions, idents = c("MS177_24hr"))

#For the control condition: 
svglite("/home/yue1118/TNBC_scRNA_analysis_workflow/TNBC_paper_data_and_figure/control_RELA_target_correlation_doublet_removed_01082025.svg", width = 8, height = 6)
plot(MS177_organoid_3_conditions_subset_1@assays$chromvar@data["MA0107.1",], MS177_organoid_3_conditions_subset_1$rela_target_score1)
chromvar_data <- MS177_organoid_3_conditions_subset_1@assays$chromvar@data
independent_variable <- chromvar_data["MA0107.1",]
model <- lm(MS177_organoid_3_conditions_subset_1$rela_target_score1 ~ as.vector(MS177_organoid_3_conditions_subset_1@assays$chromvar@data["MA0107.1",]))
summary_model <- summary(model)
summary_model$r.squared
summary_model$coefficients
abline(lm(MS177_organoid_3_conditions_subset_1$rela_target_score1 ~ as.vector(MS177_organoid_3_conditions_subset_1@assays$chromvar@data["MA0107.1",])), col = "hotpink2", lwd = 4)
dev.off()

#For the JQ1 condition:
svglite("/home/yue1118/TNBC_scRNA_analysis_workflow/TNBC_paper_data_and_figure/JQ1_RELA_target_correlation_doublet_removed_01082025.svg", width = 8, height = 6)
plot(MS177_organoid_3_conditions_subset_2@assays$chromvar@data["MA0107.1",], MS177_organoid_3_conditions_subset_2$rela_target_score1)
chromvar_data <- MS177_organoid_3_conditions_subset_2@assays$chromvar@data
independent_variable <- chromvar_data["MA0107.1",]
model <- lm(MS177_organoid_3_conditions_subset_2$rela_target_score1 ~ as.vector(MS177_organoid_3_conditions_subset_2@assays$chromvar@data["MA0107.1",]))
summary_model <- summary(model)
summary_model$r.squared
summary_model$coefficients
abline(lm(MS177_organoid_3_conditions_subset_2$rela_target_score1 ~ as.vector(MS177_organoid_3_conditions_subset_2@assays$chromvar@data["MA0107.1",])), col = "hotpink2", lwd = 4)
dev.off()

#For the MS177 condition:
svglite("/home/yue1118/TNBC_scRNA_analysis_workflow/TNBC_paper_data_and_figure/MS177_RELA_target_correlation_doublet_removed_01082025.svg", width = 8, height = 6)
plot(MS177_organoid_3_conditions_subset_3@assays$chromvar@data["MA0107.1",], MS177_organoid_3_conditions_subset_3$rela_target_score1)
chromvar_data <- MS177_organoid_3_conditions_subset_3@assays$chromvar@data
independent_variable <- chromvar_data["MA0107.1",]
model <- lm(MS177_organoid_3_conditions_subset_3$rela_target_score1 ~ as.vector(MS177_organoid_3_conditions_subset_3@assays$chromvar@data["MA0107.1",]))
summary_model <- summary(model)
summary_model$r.squared
summary_model$coefficients
abline(lm(MS177_organoid_3_conditions_subset_3$rela_target_score1 ~ as.vector(MS177_organoid_3_conditions_subset_3@assays$chromvar@data["MA0107.1",])), col = "hotpink2", lwd = 4)
dev.off()


#For the violin plot in Fig 5J
custom_colors <- c("epithelial_like" = "#078992", "stem-cell_like" = "#9E6300")
DefaultAssay(MS177_organoid_3_conditions) <- "chromvar"
Idents(MS177_organoid_3_conditions) <- MS177_organoid_3_conditions$long_condition
MS177_organoid_3_conditions$group <- factor(MS177_organoid_3_conditions$long_condition, levels = c("control_24hr_epithelial_like",  "control_24hr_stem-cell_like", "JQ1_24hr_epithelial_like", "JQ1_24hr_stem-cell_like", "MS177_24hr_epithelial_like",  "MS177_24hr_stem-cell_like"))
Idents(MS177_organoid_3_conditions) <- MS177_organoid_3_conditions$dataset
VlnPlot(MS177_organoid_3_conditions, features = c("MA0107.1"), split.by="group", pt.size = 0, cols = custom_colors)
ggsave("/home/yue1118/TNBC_scRNA_analysis_workflow/TNBC_paper_data_and_figure/RELA_motif_accessibility_2_clusters_3_conditions_doublet_removed_01082025.svg", width = 10, height = 5, dpi = 300)


Idents(MS177_organoid_3_conditions) <- "long_condition"
MS177_organoid_3_conditions_subset_1 <- subset(MS177_organoid_3_conditions, idents = c("control_24hr_epithelial_like"))
MS177_organoid_3_conditions_subset_2 <- subset(MS177_organoid_3_conditions, idents = c("control_24hr_stem-cell_like"))
MS177_organoid_3_conditions_subset_3 <- subset(MS177_organoid_3_conditions, idents = c("JQ1_24hr_epithelial_like"))
MS177_organoid_3_conditions_subset_4 <- subset(MS177_organoid_3_conditions, idents = c("JQ1_24hr_stem-cell_like"))
MS177_organoid_3_conditions_subset_5 <- subset(MS177_organoid_3_conditions, idents = c("MS177_24hr_epithelial_like"))
MS177_organoid_3_conditions_subset_6 <- subset(MS177_organoid_3_conditions, idents = c("MS177_24hr_stem-cell_like"))

vector1 <- as.vector(MS177_organoid_3_conditions_subset_1@assays$chromvar@data["MA0107.1",])
vector2 <- as.vector(MS177_organoid_3_conditions_subset_3@assays$chromvar@data["MA0107.1",])
wilcox.test(vector1, vector2, paired = F)

vector1 <- as.vector(MS177_organoid_3_conditions_subset_1@assays$chromvar@data["MA0107.1",])
vector2 <- as.vector(MS177_organoid_3_conditions_subset_5@assays$chromvar@data["MA0107.1",])
wilcox.test(vector1, vector2, paired = F)

vector1 <- as.vector(MS177_organoid_3_conditions_subset_2@assays$chromvar@data["MA0107.1",])
vector2 <- as.vector(MS177_organoid_3_conditions_subset_4@assays$chromvar@data["MA0107.1",])
wilcox.test(vector1, vector2, paired = F)

vector1 <- as.vector(MS177_organoid_3_conditions_subset_2@assays$chromvar@data["MA0107.1",])
vector2 <- as.vector(MS177_organoid_3_conditions_subset_6@assays$chromvar@data["MA0107.1",])
wilcox.test(vector1, vector2, paired = F)

#For Fig.5K

#To systematically identify TFs that are involved in NFKB expansion after MS177 treatment, we focus on the TFs whose chromatin binding motifs are significantly changed after either MS177 or JQ1 treatment (this include
#TF motifs are increased or decreased after MS177 or JQ1 treatment)
#We only focus on motifs associated with single TFs. TF heterodimers such as BATF::JUN are not considered in this analysis
#Then for each TF, we go to the EnrichR database to look for their target genes based on different TF datasets.
#The collected dataset is stored at '/home/yue1118/TNBC_scRNA_analysis_workflow/TNBC_paper_data_and_figure/TF_binding_motifs_target_genes_increased_after_MS177_avg_difference_more_than_0.25_12282024.csv'
#and '/home/yue1118/TNBC_scRNA_analysis_workflow/TNBC_paper_data_and_figure/TF_motifs_decreased_more_than_0.25_avg_difference_after_MS177_or_JQ1_12282024.csv'
#and '/home/yue1118/TNBC_scRNA_analysis_workflow/TNBC_paper_data_and_figure/TF_motifs_increased_after_JQ1_with_avg_difference_more_than_0.25_12282024.csv'


#To get information for all TF motifs in this organoid object:
pfms <- getMatrixSet(JASPAR2020,opts = list(species = 9606, all_versions = FALSE))
list_names <- c()
list_motifs <- rownames(MS177_organoid_3_conditions@assays$chromvar@data)
for (i in list_motifs) {
  list_names <- c(list_names, pfms@listData[[i]]@name)
}
length(list_names)
length(unique(list_names))
df_tf <- data.frame(Name = list_names, Motif = list_motifs)

#To extract TFs whose chromatin binding motifs are upregulated after MS177 treatment 
DefaultAssay(MS177_organoid_3_conditions) <- 'chromvar'
Idents(MS177_organoid_3_conditions) <- MS177_organoid_3_conditions$dataset
MS177_motif_markers <- FindMarkers(MS177_organoid_3_conditions, ident.1 = "MS177_24hr", ident.2 =  "control_24hr", min.pct = 0.25, logfc.threshold = 0.10, recorrect_umi=F, mean.fxn = rowMeans,fc.name = "avg_diff")
df1 <- MS177_motif_markers[(MS177_motif_markers$p_val_adj < 0.05) & (MS177_motif_markers$avg_diff >0.25),]
df1$Motif <- rownames(df1)
df2 <- merge(df1, df_tf, by= "Motif")
df2 <- df2[,c("Name", "avg_diff", "Motif", "p_val_adj")]
colnames(df2) <- c("TF", "avg_diff", "Motif", "motif_p_val_adj")
df3 <- df2
View(df3)

#To extract TFs whose chromatin binding motifs are upregulated after JQ1 treatment 
DefaultAssay(MS177_organoid_3_conditions) <- 'chromvar'
Idents(MS177_organoid_3_conditions) <- MS177_organoid_3_conditions$dataset
JQ1_motif_markers <- FindMarkers(MS177_organoid_3_conditions, ident.1 = "JQ1_24hr", ident.2 =  "control_24hr", min.pct = 0.25, logfc.threshold = 0.10, recorrect_umi=F, mean.fxn = rowMeans,fc.name = "avg_diff")
df1 <- JQ1_motif_markers[(JQ1_motif_markers$p_val_adj < 0.05) & (JQ1_motif_markers$avg_diff > 0.25),]
df1$Motif <- rownames(df1)
df2 <- merge(df1, df_tf, by= "Motif")
df2 <- df2[,c("Name", "avg_diff", "Motif", "p_val_adj")]
colnames(df2) <- c("TF", "avg_diff", "Motif", "motif_p_val_adj")


#To extract TFs whose chromatin binding motifs are down-regulated after JQ1 or MS177 treatment
DefaultAssay(MS177_organoid_3_conditions) <- 'chromvar'
Idents(MS177_organoid_3_conditions) <- MS177_organoid_3_conditions$dataset
JQ1_motif_markers <- FindMarkers(MS177_organoid_3_conditions, ident.1 = "JQ1_24hr", ident.2 =  "control_24hr", min.pct = 0.25, logfc.threshold = 0.10, recorrect_umi=F, mean.fxn = rowMeans,fc.name = "avg_diff")
MS177_motif_markers <- FindMarkers(MS177_organoid_3_conditions, ident.1 = "MS177_24hr", ident.2 =  "control_24hr", min.pct = 0.25, logfc.threshold = 0.10, recorrect_umi=F, mean.fxn = rowMeans,fc.name = "avg_diff")
df1 <- JQ1_motif_markers[(JQ1_motif_markers$p_val_adj < 0.05) & (JQ1_motif_markers$avg_diff < -0.25),]
df2 <- MS177_motif_markers[(MS177_motif_markers$p_val_adj < 0.05) & (MS177_motif_markers$avg_diff < -0.25),]
df3 <- rbind(df1, df2)
View(df3)

df3$Motif <- rownames(df3)
df3 <- df3 %>% distinct(Motif, .keep_all = TRUE)

df4 <- merge(df3, df_tf, by= "Motif")
df4 <- df4[,c("Name", "avg_diff", "Motif", "p_val_adj")]
colnames(df4) <- c("TF", "avg_diff", "Motif", "motif_p_val_adj")
View(df4[,c("TF", "Motif")])



#MsigDB TNF signaling via NF-kappaB
list_module_NFkB <- c("ABCA1",	"ACKR3",	"AREG",	"ATF3","ATP2B1",	"B4GALT1",	"B4GALT5",	"BCL2A1",	"BCL3",	"BCL6",	"BHLHE40",	"BIRC2",	"BIRC3",	"BMP2",	"BTG1",	"BTG2",	"BTG3",	"CCL2",	"CCL20",	"CCL4","CCL5",
                      "CCN1",	"CCND1",	"CCNL1",	"CCRL2",	"CD44",	"CD69",	"CD80",	"CD83",	"CDKN1A",	"CEBPB",	"CEBPD",	"CFLAR",	"CLCF1",	"CSF1",	"CSF2",	"CXCL1",	"CXCL10",	"CXCL11","CXCL2",	"CXCL3",	"CXCL6",	"DENND5A",	"DNAJB4",	"DRAM1",	"DUSP1",	"DUSP2",	"DUSP4",	"DUSP5",	"EDN1",
                      "EFNA1",	"EGR1",	"EGR2",	"EGR3",	"EHD1",	"EIF1",	"ETS2",	"F2RL1",	"F3",	"FJX1",	"FOS",	"FOSB",	"FOSL1",	"FOSL2",	"FUT4",	"G0S2",	"GADD45A",	"GADD45B",	"GCH1",	"GEM",	"GFPT2","GPR183","HBEGF",	"HES1",	"ICAM1",	"ICOSLG",	"ID2",	"IER2",	"IER3",	"IER5",	"IFIH1",	"IFIT2",	"IFNGR2", "IL12B",	"IL15RA",	"IL18",	"IL1A", "IL1B",	 "IL23A",	"IL6",	"IL6ST",	"IL7R",	
                      "INHBA",	"IRF1",	"IRS2",	"JAG1",	"JUN",	"JUNB",	"KDM6B",	"KLF10",	"KLF2",	"KLF4",	"KLF6", "KLF9", "KYNU",	"LAMB3",	"LDLR", "LIF",	"LITAF", "MAFF",	"MAP2K3",	"MAP3K8",	"MARCKS",	"MCL1", "MSC", "MXD1",	"MYC",	"NAMPT",	"NFAT5", "NFE2L2",	"NFIL3",	"NFKB1",	"NFKB2",	"NFKBIA",	"NFKBIE",	"NINJ1",	"NR4A1", "NR4A2",	"NR4A3",	"OLR1", "PANX1",	"PDE4B", "PDLIM5",	
                      "PER1",	"PFKFB3",	"PHLDA1",	"PHLDA2",	"PLAU",	"PLAUR",	"PLEK",	"PLK2",	"PLPP3",	"PMEPA1",	"PNRC1",	"PPP1R15A", "PTGER4",	"PTGS2",	"PTPRE",	"PTX3",	"RCAN1",	"REL", "RELA",	"RELB",	"RHOB",	"RIGI",	"RIPK2",	"RNF19B",	"SAT1",	"SDC4",	"SERPINB2",	"SERPINB8", "SERPINE1",	"SGK1",	"SIK1",	"SLC16A6",	"SLC2A3",	"SLC2A6",	"SMAD3",	"SNN",	"SOCS3", "SOD2", "SPHK1",	"SPSB1",	"SQSTM1",	"STAT5A",	"TANK",	"TAP1",	"TGIF1",	"TIPARP",	"TLR2",	
                      "TNC",	"TNF",	"TNFAIP2",	"TNFAIP3",	"TNFAIP6",	"TNFAIP8",	"TNFRSF9",	"TNFSF9",	"TNIP1",	"TNIP2", 
                      "TRAF1",	"TRIB1", "TRIP10",	"TSC22D1",	"TUBB2A",	"VEGFA",	"YRDC",	"ZBTB10",	"ZC3H12A",	"ZFP36") 


#To compute the correlation for each section:

df1 <- read.csv('/home/yue1118/TNBC_scRNA_analysis_workflow/TNBC_paper_data_and_figure/TF_binding_motifs_target_genes_increased_after_MS177_avg_difference_more_than_0.25_12282024.csv')
View(df1)
list_r_squared <- c()
list_motif_name <- c()
list_TF_name <- c()
list_overlap_length <- c()
list_condition <- c()
list_pvalue <- c()
DefaultAssay(MS177_organoid_3_conditions) <- "SCT"
Idents(MS177_organoid_3_conditions) <- "dataset"

for(i in c(1:dim(df1)[1])){
  row_vector <- as.vector(unlist(df1[i, 4:ncol(df1)], use.names = FALSE))
  row_vector <- row_vector[grepl("[A-Za-z]", row_vector)]
  list_overlap <- list(intersect(row_vector, list_module_NFkB))
  for (j in c("control_24hr", "JQ1_24hr","MS177_24hr")) {
    try({
      print(i)    
      print(df1[i,2])
      MS177_organoid_3_conditions <- AddModuleScore(MS177_organoid_3_conditions,features = list_overlap, pool = NULL,nbin = 10,ctrl = 5,k = FALSE,assay = NULL,name = "target_gene_score",seed = 1,search = FALSE,slot = "data")
      MS177_organoid_3_conditions_subset <- subset(MS177_organoid_3_conditions, idents = c(j))
      model <- lm(MS177_organoid_3_conditions_subset$target_gene_score1 ~ as.vector(MS177_organoid_3_conditions_subset@assays$chromvar@data[df1[i,1],]))
      summary_model <- summary(model)
      list_r_squared <- c(list_r_squared, summary_model$r.squared)
      list_pvalue <- c(list_pvalue, summary_model$coefficients[2,4])
      list_motif_name <- c(list_motif_name, df1[i,1])
      list_TF_name <- c(list_TF_name, df1[i,2])
      list_overlap_length <- c(list_overlap_length, length(intersect(row_vector, list_module_NFkB)))
      list_condition <- c(list_condition, j)
    }, silent = T)}
}

df1_all <- data.frame(r_squared = list_r_squared, motif_name = list_motif_name, TF_name = list_TF_name, overlap_length = list_overlap_length, condition = list_condition, p_value = list_pvalue)
View(df1_all)
df1_control <- df1_all[df1_all$condition == "control_24hr",]
df1_JQ1 <- df1_all[df1_all$condition == "JQ1_24hr",]
df1_MS177 <- df1_all[df1_all$condition == "MS177_24hr",]

empty_matrix <- matrix(ncol = 2, nrow = 95)
cor_df <- data.frame(empty_matrix)

cor_df$X1 <- df1_JQ1$r_squared - df1_control$r_squared
cor_df$X2 <- df1_MS177$r_squared - df1_control$r_squared
colnames(cor_df) <- c("delta_JQ1", "delta_MS177")
cor_df$TF <- df1_MS177$TF_name
cor_df$correlation <- df1_MS177$r_squared


df2 <- read.csv('/home/yue1118/TNBC_scRNA_analysis_workflow/TNBC_paper_data_and_figure/TF_motifs_increased_after_JQ1_with_avg_difference_more_than_0.25_12282024.csv')
df3 <- df2[!df2$X.1 %in% df1$Gene,]
View(df3)
list_r_squared_2 <- c()
list_motif_name_2 <- c()
list_TF_name_2 <- c()
list_pvalue_2 <- c()
list_overlap_length_2 <- c()
list_condition_2 <- c()
for( i in c(1:dim(df3)[1])){
  row_vector <- as.vector(unlist(df3[i, 4:ncol(df3)], use.names = FALSE))
  row_vector <- row_vector[grepl("[A-Za-z]", row_vector)]
  list_overlap <- list(intersect(row_vector, list_module_NFkB))
  for (j in c("control_24hr", "JQ1_24hr","MS177_24hr")) {
    try({                                                                                                                                                             #loop_NFKB_overlap_score, "
      MS177_organoid_3_conditions <- AddModuleScore(MS177_organoid_3_conditions,features = list_overlap, pool = NULL, nbin = 10,ctrl = 5,k = FALSE,assay = NULL,name = "target_gene_score",seed = 1,search = FALSE,slot = "data")
      MS177_organoid_3_conditions_subset <- subset(MS177_organoid_3_conditions, idents = c(j))
      print(i)
      print(df3[i,2])
      model <- lm(MS177_organoid_3_conditions_subset$target_gene_score1 ~ as.vector(MS177_organoid_3_conditions_subset@assays$chromvar@data[df3[i,1],]))
      summary_model <- summary(model)
      list_r_squared_2 <- c(list_r_squared_2, summary_model$r.squared)
      list_motif_name_2 <- c(list_motif_name_2, df3[i,1])
      list_TF_name_2 <- c(list_TF_name_2, df3[i,2])
      list_pvalue_2 <- c(list_pvalue_2, summary_model$coefficients[2,4])
      list_overlap_length_2 <- c(list_overlap_length_2, length(intersect(row_vector, list_module_NFkB)))
      list_condition_2 <- c(list_condition_2, j)
    }, silent = T)}
}

df2_all <- data.frame(r_squared = list_r_squared_2, motif_name = list_motif_name_2, TF_name = list_TF_name_2, overlap_length = list_overlap_length_2, condition = list_condition_2, p_value = list_pvalue_2)
View(df2_all)
df2_control <- df2_all[df2_all$condition == "control_24hr",]
df2_JQ1 <- df2_all[df2_all$condition == "JQ1_24hr",]
df2_MS177 <- df2_all[df2_all$condition == "MS177_24hr",]


empty_matrix <- matrix(ncol = 2, nrow = 59)
cor_df_2 <- data.frame(empty_matrix)

cor_df_2$X1 <- df2_JQ1$r_squared - df2_control$r_squared
cor_df_2$X2 <- df2_MS177$r_squared - df2_control$r_squared
View(cor_df_2)
colnames(cor_df_2) <- c("delta_JQ1", "delta_MS177")
cor_df_2$TF <- df2_MS177$TF_name
cor_df_2$correlation <- df2_MS177$r_squared
View(cor_df_2)


df_decrease <- read.csv('../TNBC_scRNA_analysis_workflow/TNBC_paper_data_and_figure/TF_motifs_decreased_more_than_0.25_avg_difference_after_MS177_or_JQ1_12282024.csv')
df_decrease <- df_decrease[!(df_decrease$X.1 %in% unique(cor_df_2$TF)), ]
df_decrease <- df_decrease[!(df_decrease$X.1 %in% unique(cor_df$TF)), ]

list_r_squared_3 <- c()
list_motif_name_3 <- c()
list_TF_name_3 <- c()
list_pvalue_3 <- c()
list_overlap_length_3 <- c()
list_condition_3 <- c()
for( i in c(1:dim(df_decrease)[1])){
  row_vector <- as.vector(unlist(df_decrease[i, 4:ncol(df_decrease)], use.names = FALSE))
  row_vector <- row_vector[grepl("[A-Za-z]", row_vector)]
  list_overlap <- list(intersect(row_vector, list_module_NFkB))
  for (j in c("control_24hr", "JQ1_24hr","MS177_24hr")){
    try({                                                                                                                                                             #loop_NFKB_overlap_score, 
      MS177_organoid_3_conditions <- AddModuleScore(MS177_organoid_3_conditions,features = list_overlap, pool = NULL,nbin = 10,ctrl = 5,k = FALSE,assay = NULL,name = "target_gene_score",seed = 1,search = FALSE,slot = "data")
      MS177_organoid_3_conditions_subset <- subset(MS177_organoid_3_conditions, idents = c(j))
      model <- lm(MS177_organoid_3_conditions_subset$target_gene_score1 ~ as.vector(MS177_organoid_3_conditions_subset@assays$chromvar@data[df_decrease[i,1],]))
      summary_model <- summary(model)
      list_r_squared_3 <- c(list_r_squared_3, summary_model$r.squared)
      list_motif_name_3 <- c(list_motif_name_3, df_decrease[i,1])
      list_pvalue_3 <- c(list_pvalue_3, summary_model$coefficients[2,4])
      list_TF_name_3 <- c(list_TF_name_3, df_decrease[i,2])
      list_overlap_length_3 <- c(list_overlap_length_3, length(intersect(row_vector, list_module_NFkB)))
      list_condition_3 <- c(list_condition_3, j)
    }, silent = T)
  }
}

df3_all <- data.frame(r_squared = list_r_squared_3, motif_name = list_motif_name_3, TF_name = list_TF_name_3, overlap_length = list_overlap_length_3, condition = list_condition_3, p_value = list_pvalue_3)
View(df3_all)
df3_control <- df3_all[df3_all$condition == "control_24hr",]
df3_JQ1 <- df3_all[df3_all$condition == "JQ1_24hr",]
df3_MS177 <- df3_all[df3_all$condition == "MS177_24hr",]

empty_matrix <- matrix(ncol = 2, nrow = 39)
cor_df_3 <- data.frame(empty_matrix)

cor_df_3$X1 <- df3_JQ1$r_squared - df3_control$r_squared
cor_df_3$X2 <- df3_MS177$r_squared - df3_control$r_squared
View(cor_df_3)
colnames(cor_df_3) <- c("delta_JQ1", "delta_MS177")
cor_df_3$TF <- df3_MS177$TF_name
cor_df_3$correlation <- df3_MS177$r_squared
cor_df_all <- rbind(cor_df, cor_df_2, cor_df_3)
cor_df_all <- cor_df_all[order(cor_df_all$correlation), ]
normalized_values <- (cor_df_all$correlation - min(cor_df_all$correlation)) / (max(cor_df_all$correlation) - min(cor_df_all$correlation))
colors <- gray(1 - normalized_values) 

write.csv(cor_df_all, "/home/yue1118/TNBC_scRNA_analysis_workflow/TNBC_paper_data_and_figure/TF_only_NFKB_target_genes_correlation_changes_across_different_conditions_JQ1_MS177_doublet_removed_01082025.csv")

df_all <- rbind(df1_all, df2_all, df3_all)
write.csv(df_all, "../TNBC_scRNA_analysis_workflow/TNBC_paper_data_and_figure/TF_only_NFKB_target_genes_correlations_across_different_conditions_JQ1_MS177_doublet_removed_01082025.csv")
df_all <- read.csv("../TNBC_scRNA_analysis_workflow/TNBC_paper_data_and_figure/TF_only_NFKB_target_genes_correlations_across_different_conditions_JQ1_MS177_doublet_removed_01082025.csv")
View(df_all)

svglite("/home/yue1118/TNBC_scRNA_analysis_workflow/TNBC_paper_data_and_figure/Global_TF_correlation_changes_NFKB_target_genes_doublet_removed_01082025.svg", width = 4, height = 4)
par(mar=c(2,2,2,2))
plot(cor_df_all$delta_JQ1, cor_df_all$delta_MS177, col = colors, pch=19, cex=1)
abline(v = 0, col = "green", lwd = 2)
abline(h = 0, col = "green", lwd = 2)
text(cor_df_all$delta_JQ1-0.005, cor_df_all$delta_MS177-0.005, labels = cor_df_all$TF, pos = 4, cex = 1, col = "red")
dev.off()

#Supplementary figures for the organoid

#To plot all potentially involved TFs in MS177 condition according to Fig 5I.
df1 <- read.csv('/home/yue1118/TNBC_scRNA_analysis_workflow/TNBC_paper_data_and_figure/TF_binding_motifs_target_genes_increased_after_MS177_avg_difference_more_than_0.25_12282024.csv')
View(df1)
DefaultAssay(MS177_organoid_3_conditions) <- "SCT"
#MsigDB TNF signaling via NF-kappaB
list_module_NFkB <- c("ABCA1",	"ACKR3",	"AREG",	"ATF3","ATP2B1",	"B4GALT1",	"B4GALT5",	"BCL2A1",	"BCL3",	"BCL6",	"BHLHE40",	"BIRC2",	"BIRC3",	"BMP2",	"BTG1",	"BTG2",	"BTG3",	"CCL2",	"CCL20",	"CCL4","CCL5",
                      "CCN1",	"CCND1",	"CCNL1",	"CCRL2",	"CD44",	"CD69",	"CD80",	"CD83",	"CDKN1A",	"CEBPB",	"CEBPD",	"CFLAR",	"CLCF1",	"CSF1",	"CSF2",	"CXCL1",	"CXCL10",	"CXCL11","CXCL2",	"CXCL3",	"CXCL6",	"DENND5A",	"DNAJB4",	"DRAM1",	"DUSP1",	"DUSP2",	"DUSP4",	"DUSP5",	"EDN1",
                      "EFNA1",	"EGR1",	"EGR2",	"EGR3",	"EHD1",	"EIF1",	"ETS2",	"F2RL1",	"F3",	"FJX1",	"FOS",	"FOSB",	"FOSL1",	"FOSL2",	"FUT4",	"G0S2",	"GADD45A",	"GADD45B",	"GCH1",	"GEM",	"GFPT2","GPR183","HBEGF",	"HES1",	"ICAM1",	"ICOSLG",	"ID2",	"IER2",	"IER3",	"IER5",	"IFIH1",	"IFIT2",	"IFNGR2", "IL12B",	"IL15RA",	"IL18",	"IL1A", "IL1B",	 "IL23A",	"IL6",	"IL6ST",	"IL7R",	
                      "INHBA",	"IRF1",	"IRS2",	"JAG1",	"JUN",	"JUNB",	"KDM6B",	"KLF10",	"KLF2",	"KLF4",	"KLF6", "KLF9", "KYNU",	"LAMB3",	"LDLR", "LIF",	"LITAF", "MAFF",	"MAP2K3",	"MAP3K8",	"MARCKS",	"MCL1", "MSC", "MXD1",	"MYC",	"NAMPT",	"NFAT5", "NFE2L2",	"NFIL3",	"NFKB1",	"NFKB2",	"NFKBIA",	"NFKBIE",	"NINJ1",	"NR4A1", "NR4A2",	"NR4A3",	"OLR1", "PANX1",	"PDE4B", "PDLIM5",	
                      "PER1",	"PFKFB3",	"PHLDA1",	"PHLDA2",	"PLAU",	"PLAUR",	"PLEK",	"PLK2",	"PLPP3",	"PMEPA1",	"PNRC1",	"PPP1R15A", "PTGER4",	"PTGS2",	"PTPRE",	"PTX3",	"RCAN1",	"REL", "RELA",	"RELB",	"RHOB",	"RIGI",	"RIPK2",	"RNF19B",	"SAT1",	"SDC4",	"SERPINB2",	"SERPINB8", "SERPINE1",	"SGK1",	"SIK1",	"SLC16A6",	"SLC2A3",	"SLC2A6",	"SMAD3",	"SNN",	"SOCS3", "SOD2", "SPHK1",	"SPSB1",	"SQSTM1",	"STAT5A",	"TANK",	"TAP1",	"TGIF1",	"TIPARP",	"TLR2",	
                      "TNC",	"TNF",	"TNFAIP2",	"TNFAIP3",	"TNFAIP6",	"TNFAIP8",	"TNFRSF9",	"TNFSF9",	"TNIP1",	"TNIP2", 
                      "TRAF1",	"TRIB1", "TRIP10",	"TSC22D1",	"TUBB2A",	"VEGFA",	"YRDC",	"ZBTB10",	"ZC3H12A",	"ZFP36") 
Idents(MS177_organoid_3_conditions) <- MS177_organoid_3_conditions$dataset

for(i in c(3,19,24,34,46,60,71,72,6,13,16,37,47,49,54,66,69)){ #For the potentially involved TFs in MS177 condition ATF3, EGR1, FOSL1, JUN, KLF4, NFKB1, RELA, RELB, BACH1, CREB1, CTCF, JUND, KLF5, MAF, MITF, PHOX2B,PPARD.
  h=df1[i,2]
  k=df1[i,1]
  row_vector <- as.vector(unlist(df1[i, 4:ncol(df1)], use.names = FALSE))
  row_vector <- row_vector[grepl("[A-Za-z]", row_vector)]
  #list_overlap <- list(intersect(row_vector, list_module_NFkB))
  list_overlap <- list(intersect(row_vector, list_module_NFkB))
  for (j in c("control_24hr", "JQ1_24hr","MS177_24hr")) {
    try({
      print(h)                                                                                                                                                  #"loop_NFKB_overlap_score", ""
      print(k)
      print(j)
      MS177_organoid_3_conditions <- AddModuleScore(MS177_organoid_3_conditions,features = list_overlap, pool = NULL,nbin = 10,ctrl = 5,k = FALSE,assay = NULL,name = "target_gene_score",seed = 1,search = FALSE,slot = "data")
      MS177_organoid_3_conditions_subset <- subset(MS177_organoid_3_conditions, idents = c(j))
      file_name_to_save <- paste0("/home/yue1118/TNBC_scRNA_analysis_workflow/TNBC_paper_data_and_figure/",j, "_",h, "_target_gene_correlation_doublet_removed_01082025.svg")
      svglite(file_name_to_save, width = 8, height = 6)
      par(mar = c(2, 2, 2, 2))
      plot(MS177_organoid_3_conditions_subset@assays$chromvar@data[k,], MS177_organoid_3_conditions_subset$target_gene_score1)
      chromvar_data <- MS177_organoid_3_conditions_subset@assays$chromvar@data
      independent_variable <- chromvar_data[k,]
      model <- lm(MS177_organoid_3_conditions_subset$target_gene_score1 ~ as.vector(MS177_organoid_3_conditions_subset@assays$chromvar@data[k,]))
      summary_model <- summary(model)
      summary_model$r.squared
      summary_model$coefficients
      abline(lm(MS177_organoid_3_conditions_subset$target_gene_score1 ~ as.vector(MS177_organoid_3_conditions_subset@assays$chromvar@data[k,])), col = "hotpink2", lwd = 20)
      dev.off()
      print(paste0("r_squared", summary_model$r.squared))
      print(paste0("pvalue", summary_model$coefficients[2,4]))
    })
  }
}


#To plot the JQ1 version of Fig 5I.

#To get the differential TF motifs after JQ1 treatment
DefaultAssay(MS177_organoid_3_conditions) <- 'chromvar'
Idents(MS177_organoid_3_conditions) <- MS177_organoid_3_conditions$dataset
JQ1_motif_markers <- FindMarkers(MS177_organoid_3_conditions, ident.1 = "JQ1_24hr", ident.2 =  "control_24hr", min.pct = 0.25, logfc.threshold = 0.10, recorrect_umi=F, mean.fxn = rowMeans,fc.name = "avg_diff")

#To get all the upregulated genes after JQ1
DefaultAssay(MS177_organoid_3_conditions) <- 'SCT'
DEG_organoid_JQ1 <- FindMarkers(MS177_organoid_3_conditions, ident.1 = "JQ1_24hr", ident.2 = "control_24hr", min.pct = 0.25)
DEG_organoid_JQ1_all_up_regulated_genes <- DEG_organoid_JQ1[(DEG_organoid_JQ1$p_val_adj < 0.05) & (DEG_organoid_JQ1$avg_log2FC > 0),]
write.csv(DEG_organoid_JQ1_all_up_regulated_genes, '/home/yue1118/TNBC_scRNA_analysis_workflow/TNBC_paper_data_and_figure/all_upregulated_genes_after_JQ1_in_the_organoid_doublet_removed_pval_adj_lower_than_0.05_and_avg_log2FC_more_than_0_01082025.csv')
#The genes included in the dataset above were input into EnrichR database to look for potentially involved TFs
#The ChEA 2022 database was used to extract all the potentially involved TFs
#The result file is named "ChEA_2022_table_using_all_the_upregulated_genes_after_JQ1_doublet_removed_01082025.txt"


#To get all the upregulated genes in NFKB pathway after JQ1
DEG_organoid_JQ1_all_up_regulated_NFKB_genes <- DEG_organoid_JQ1[(DEG_organoid_JQ1$p_val_adj < 0.05) & (DEG_organoid_JQ1$avg_log2FC > 0) & (rownames(DEG_organoid_JQ1) %in% NFKB_gene_sets$gene),]
write.csv(DEG_organoid_JQ1_all_up_regulated_NFKB_genes, '/home/yue1118/TNBC_scRNA_analysis_workflow/TNBC_paper_data_and_figure/all_upregulated_NFKB_genes_after_JQ1_in_the_organoid_doublet_removed_pval_adj_lower_than_0.05_and_avg_log2FC_more_than_0_01082025.csv')
#The genes included in the dataset above were input into EnrichR database to look for potentially involved TFs
#The ChEA 2022 database was used to extract all the potentially involved TFs
#The result file is named "ChEA_2022_table_using_all_the_upregulated_NFKB_genes_after_JQ1_doublet_removed_01082025.txt"



pfms <- getMatrixSet(JASPAR2020,opts = list(species = 9606, all_versions = FALSE))

list_names <- c()
list_motifs <- rownames(MS177_organoid_3_conditions@assays$chromvar@data)
for (i in list_motifs) {
  list_names <- c(list_names, pfms@listData[[i]]@name)
}
length(list_names)
length(unique(list_names))
df_tf <- data.frame(Name = list_names, Motif = list_motifs)

df1 <- JQ1_motif_markers[(JQ1_motif_markers$p_val_adj < 0.05) & (JQ1_motif_markers$avg_diff > 0.25),]
df1$Motif <- rownames(df1)

df2 <- merge(df1, df_tf, by= "Motif")
df2 <- df2[,c("Name", "avg_diff", "Motif", "p_val_adj")]
colnames(df2) <- c("TF", "avg_diff", "Motif", "motif_p_val_adj")
View(df2)

list_tf_with_motf_changes <- df2$TF
list_NFKB_motifs <- intersect(NFKB_gene_sets$gene, list_tf_with_motf_changes)

#For all the up-regulated genes after MS177 in general  #To get the df_chea dataset as shown below, we get all the up-regulated
df_chea <- read.table('../TNBC_scRNA_analysis_workflow/TNBC_paper_data_and_figure/ChEA_2022_table_using_all_the_upregulated_genes_after_JQ1_doublet_removed_01082025.txt',header = TRUE, sep = "\t", fileEncoding = "UTF-8")
df_chea$TF <- sapply(str_split(df_chea$Term, " "), `[`, 1)
df_chea <- df_chea[df_chea$Adjusted.P.value < 0.05,]
df_chea <- df_chea[grep("Human", df_chea$Term),]
df_chea_2 <- df_chea[df_chea$TF %in% df2$TF,]
df_chea_2 <- df_chea_2[!duplicated(df_chea_2$TF), ]
df_chea_2 <- df_chea_2[, c("Term", "Adjusted.P.value", "TF")]
colnames(df_chea_2) <- c("global_DEG_TF_Term", "global_DEG_TF_enrichment_adjusted.P.value", "TF")
df3 <- merge(df2, df_chea_2, by= "TF", all.x=T)
View(df3)
View(df_chea)

df_chea_nfkb <- read.table('../TNBC_scRNA_analysis_workflow/TNBC_paper_data_and_figure/ChEA_2022_table_using_all_the_upregulated_NFKB_genes_after_JQ1_doublet_removed_01082025.txt',header = TRUE, sep = "\t", fileEncoding = "UTF-8")
df_chea_nfkb$TF <- sapply(str_split(df_chea_nfkb$Term, " "), `[`, 1)
df_chea_nfkb <- df_chea_nfkb[df_chea_nfkb$Adjusted.P.value < 0.05,]
df_chea_nfkb <- df_chea_nfkb[grep("Human", df_chea_nfkb$Term),]
df_chea_nfkb_2 <- df_chea_nfkb[df_chea_nfkb$TF %in% df2$TF,]
df_chea_nfkb_2 <- df_chea_nfkb_2[!duplicated(df_chea_nfkb_2$TF), ]
df_chea_nfkb_2 <- df_chea_nfkb_2[, c("Term", "Adjusted.P.value", "TF")]
colnames(df_chea_nfkb_2) <- c("NFKB_DEG_TF_Term", "NFKB_DEG_TF_enrichment_adjusted.P.value", "TF")
df4 <- merge(df3, df_chea_nfkb_2, by= "TF", all.x=T)
View(df4)

df4 <- df4[!(is.na(df4$global_DEG_TF_Term) & is.na(df4$NFKB_DEG_TF_Term)), ]
df4$NFKB_or_not <- ifelse(df4$TF %in% list_NFKB_motifs, "Yes", "No")
df4$NFKB_or_not <- factor(df4$NFKB_or_not, levels = c("Yes", "No"))
df4 <- df4[, c("TF", "Motif","avg_diff", "global_DEG_TF_enrichment_adjusted.P.value", "NFKB_DEG_TF_enrichment_adjusted.P.value", "NFKB_or_not")]
df4 <- df4[order(df4$NFKB_or_not), ]
df4$global_DEG_TF_enrichment_adjusted.P.value <- -log(df4$global_DEG_TF_enrichment_adjusted.P.value)
df4$NFKB_DEG_TF_enrichment_adjusted.P.value <- -log(df4$NFKB_DEG_TF_enrichment_adjusted.P.value)
rownames(df4) <- df4$TF


View(df4)
color_palette <- colorRampPalette(colors = c("white", "red"))(100)
color_palette_2 <- colorRampPalette(colors = c("white", "orange"))(100)
svglite('/home/yue1118/TNBC_scRNA_analysis_workflow/TNBC_paper_data_and_figure/JQ1_NFKB_enriched_TFs_doublet_removed_01082025.svg', width = 6, height = 12)
pheatmap(df4[, c("global_DEG_TF_enrichment_adjusted.P.value", "NFKB_DEG_TF_enrichment_adjusted.P.value")], cluster_rows = F, cluster_cols = F, color = color_palette)
dev.off()

svglite('/home/yue1118/TNBC_scRNA_analysis_workflow/TNBC_paper_data_and_figure/JQ1_upregulated_TF_motifs_doublet_removed_01082025.svg', width = 6, height = 12)
pheatmap(as.data.frame(df4[, c("avg_diff")]), cluster_rows = F, cluster_cols = F, color = color_palette_2)
dev.off()



#To plot all potentially involved TFs in JQ1 condition according to Fig 5I.
df1 <- read.csv('/home/yue1118/TNBC_scRNA_analysis_workflow/TNBC_paper_data_and_figure/TF_motifs_increased_after_JQ1_with_avg_difference_more_than_0.25_12282024.csv')
View(df1)
DefaultAssay(MS177_organoid_3_conditions) <- "SCT"
#MsigDB TNF signaling via NF-kappaB
list_module_NFkB <- c("ABCA1",	"ACKR3",	"AREG",	"ATF3","ATP2B1",	"B4GALT1",	"B4GALT5",	"BCL2A1",	"BCL3",	"BCL6",	"BHLHE40",	"BIRC2",	"BIRC3",	"BMP2",	"BTG1",	"BTG2",	"BTG3",	"CCL2",	"CCL20",	"CCL4","CCL5",
                      "CCN1",	"CCND1",	"CCNL1",	"CCRL2",	"CD44",	"CD69",	"CD80",	"CD83",	"CDKN1A",	"CEBPB",	"CEBPD",	"CFLAR",	"CLCF1",	"CSF1",	"CSF2",	"CXCL1",	"CXCL10",	"CXCL11","CXCL2",	"CXCL3",	"CXCL6",	"DENND5A",	"DNAJB4",	"DRAM1",	"DUSP1",	"DUSP2",	"DUSP4",	"DUSP5",	"EDN1",
                      "EFNA1",	"EGR1",	"EGR2",	"EGR3",	"EHD1",	"EIF1",	"ETS2",	"F2RL1",	"F3",	"FJX1",	"FOS",	"FOSB",	"FOSL1",	"FOSL2",	"FUT4",	"G0S2",	"GADD45A",	"GADD45B",	"GCH1",	"GEM",	"GFPT2","GPR183","HBEGF",	"HES1",	"ICAM1",	"ICOSLG",	"ID2",	"IER2",	"IER3",	"IER5",	"IFIH1",	"IFIT2",	"IFNGR2", "IL12B",	"IL15RA",	"IL18",	"IL1A", "IL1B",	 "IL23A",	"IL6",	"IL6ST",	"IL7R",	
                      "INHBA",	"IRF1",	"IRS2",	"JAG1",	"JUN",	"JUNB",	"KDM6B",	"KLF10",	"KLF2",	"KLF4",	"KLF6", "KLF9", "KYNU",	"LAMB3",	"LDLR", "LIF",	"LITAF", "MAFF",	"MAP2K3",	"MAP3K8",	"MARCKS",	"MCL1", "MSC", "MXD1",	"MYC",	"NAMPT",	"NFAT5", "NFE2L2",	"NFIL3",	"NFKB1",	"NFKB2",	"NFKBIA",	"NFKBIE",	"NINJ1",	"NR4A1", "NR4A2",	"NR4A3",	"OLR1", "PANX1",	"PDE4B", "PDLIM5",	
                      "PER1",	"PFKFB3",	"PHLDA1",	"PHLDA2",	"PLAU",	"PLAUR",	"PLEK",	"PLK2",	"PLPP3",	"PMEPA1",	"PNRC1",	"PPP1R15A", "PTGER4",	"PTGS2",	"PTPRE",	"PTX3",	"RCAN1",	"REL", "RELA",	"RELB",	"RHOB",	"RIGI",	"RIPK2",	"RNF19B",	"SAT1",	"SDC4",	"SERPINB2",	"SERPINB8", "SERPINE1",	"SGK1",	"SIK1",	"SLC16A6",	"SLC2A3",	"SLC2A6",	"SMAD3",	"SNN",	"SOCS3", "SOD2", "SPHK1",	"SPSB1",	"SQSTM1",	"STAT5A",	"TANK",	"TAP1",	"TGIF1",	"TIPARP",	"TLR2",	
                      "TNC",	"TNF",	"TNFAIP2",	"TNFAIP3",	"TNFAIP6",	"TNFAIP8",	"TNFRSF9",	"TNFSF9",	"TNIP1",	"TNIP2", 
                      "TRAF1",	"TRIB1", "TRIP10",	"TSC22D1",	"TUBB2A",	"VEGFA",	"YRDC",	"ZBTB10",	"ZC3H12A",	"ZFP36") 
Idents(MS177_organoid_3_conditions) <- MS177_organoid_3_conditions$dataset

for(i in c(16,17,28,32,41,47,48,49,50,52,91,109)){ #For the potentially involved TFs in JQ1 condition ELK1, ELK3, FLI1, FOXA2, FOXO3, GATA1, GATA2, GATA3, GATA4, GATA6, RUNX2, YY1
  h=df1[i,2]
  k=df1[i,1]
  row_vector <- as.vector(unlist(df1[i, 4:ncol(df1)], use.names = FALSE))
  row_vector <- row_vector[grepl("[A-Za-z]", row_vector)]
  #list_overlap <- list(intersect(row_vector, list_module_NFkB))
  list_overlap <- list(intersect(row_vector, list_module_NFkB))
  for (j in c("control_24hr", "JQ1_24hr","MS177_24hr")) {
    try({
      print(h)                                                                                                                                                  #"loop_NFKB_overlap_score", ""
      print(k)
      print(j)
      MS177_organoid_3_conditions <- AddModuleScore(MS177_organoid_3_conditions,features = list_overlap, pool = NULL,nbin = 10,ctrl = 5,k = FALSE,assay = NULL,name = "target_gene_score",seed = 1,search = FALSE,slot = "data")
      MS177_organoid_3_conditions_subset <- subset(MS177_organoid_3_conditions, idents = c(j))
      file_name_to_save <- paste0("/home/yue1118/TNBC_scRNA_analysis_workflow/TNBC_paper_data_and_figure/",j, "_",h, "_target_gene_correlation_doublet_removed_01082025.svg")
      svglite(file_name_to_save, width = 8, height = 6)
      par(mar = c(2, 2, 2, 2))
      plot(MS177_organoid_3_conditions_subset@assays$chromvar@data[k,], MS177_organoid_3_conditions_subset$target_gene_score1)
      chromvar_data <- MS177_organoid_3_conditions_subset@assays$chromvar@data
      independent_variable <- chromvar_data[k,]
      model <- lm(MS177_organoid_3_conditions_subset$target_gene_score1 ~ as.vector(MS177_organoid_3_conditions_subset@assays$chromvar@data[k,]))
      summary_model <- summary(model)
      summary_model$r.squared
      summary_model$coefficients
      abline(lm(MS177_organoid_3_conditions_subset$target_gene_score1 ~ as.vector(MS177_organoid_3_conditions_subset@assays$chromvar@data[k,])), col = "hotpink2", lwd = 20)
      dev.off()
      print(paste0("r_squared", summary_model$r.squared))
      print(paste0("pvalue", summary_model$coefficients[2,4]))
    })
  }
}

#For Fig 5L:
Idents(MS177_organoid_3_conditions) <- "dataset"
MS177_organoid_3_conditions$group <- factor(MS177_organoid_3_conditions$long_condition, levels = c("control_24hr_epithelial_like",  "control_24hr_stem-cell_like", "JQ1_24hr_epithelial_like", "JQ1_24hr_stem-cell_like", "MS177_24hr_epithelial_like",  "MS177_24hr_stem-cell_like"))

DefaultAssay(MS177_organoid_3_conditions) <- "SCT"
VlnPlot(MS177_organoid_3_conditions, features = c("MARCKS", "EDN1", "IRF1"), split.by = "group", pt.size = 0, stack = T, flip = T)
ggsave("/home/yue1118/TNBC_scRNA_analysis_workflow/TNBC_paper_data_and_figure/RELA_3_target_gene_expressions_with_peaks_changes_in_promoter_region_in_3_conditions_split_by_2_clusters_doublet_removed_01082025.svg", width = 10, height = 5)

#For the supplementary figure for 5L:
VlnPlot(MS177_organoid_3_conditions, features = c("CEBPD", "LAMB3", "AREG", "B4GALT1", "BTG1", "CXCL1", "EIF1", "FOSL2"), split.by = "group", pt.size = 0, stack = T, flip = T)
ggsave("/home/yue1118/TNBC_scRNA_analysis_workflow/TNBC_paper_data_and_figure/TNBC_paper_Fig5_supple_figure_part1_doublet_removed_01082025.svg", width = 15, height = 10)

VlnPlot(MS177_organoid_3_conditions, features = c("GADD45A", "HBEGF", "IER3", "INHBA", "IRS2", "JAG1","MYC"), split.by = "group", pt.size = 0, stack = T, flip = T)
ggsave("/home/yue1118/TNBC_scRNA_analysis_workflow/TNBC_paper_data_and_figure/TNBC_paper_Fig5_supple_figure_part2_doublet_removed_01082025.svg", width = 15, height = 10)

VlnPlot(MS177_organoid_3_conditions, features = c("NAMPT", "NFAT5", "NFE2L2", "NFKB1", "NFKBIA", "PLK2", "RCAN1"), split.by = "group", pt.size = 0, stack = T, flip = T)
ggsave("/home/yue1118/TNBC_scRNA_analysis_workflow/TNBC_paper_data_and_figure/TNBC_paper_Fig5_supple_figure_part3_doublet_removed_01082025.svg", width = 15, height = 10)


#For for Fig.5M
Idents(MS177_organoid_3_conditions) <- "dataset"
differential_peaks <- FindMarkers(
  object = MS177_organoid_3_conditions,
  ident.1 = "MS177_24hr",      # Group1, e.g., a cell type or condition
  ident.2 = "control_24hr",      # Group2, e.g., another cell type or condition
  assay = "peaks",          # Ensure you're analyzing the ATAC assay
  test.use = "LR",         # Use logistic regression, often suitable for sparse data
  min.pct = 0.1,           # Minimum percent of cells with the peak accessible
  logfc.threshold = 0.25   # Set log fold-change threshold for differential accessibility
)
#write.csv(differential_peaks, '/home/yue1118/TNBC_scRNA_analysis_workflow/TNBC_paper_data_and_figure/differential_peaks_between_MS177_and_control_in_the_organoid_doublet_removed_01082025.csv')
View(differential_peaks[grepl("chr5-13249", rownames(differential_peaks)), ]) #chr5-132495686-132498505 p_val_adj=4.44314e-07


Idents(MS177_organoid_3_conditions) <- "long_condition"
DefaultAssay(MS177_organoid_3_conditions) <- "peaks"
differential_peaks_2 <- FindMarkers(
  object = MS177_organoid_3_conditions,
  ident.1 = "MS177_24hr_stem-cell_like",      # Group1, e.g., a cell type or condition
  ident.2 = "control_24hr_stem-cell_like",      # Group2, e.g., another cell type or condition
  assay = "peaks",          # Ensure you're analyzing the ATAC assay
  test.use = "LR",         # Use logistic regression, often suitable for sparse data
  min.pct = 0.1,           # Minimum percent of cells with the peak accessible
  logfc.threshold = 0.25   # Set log fold-change threshold for differential accessibility
)
write.csv(differential_peaks_2, '/home/yue1118/TNBC_scRNA_analysis_workflow/TNBC_paper_data_and_figure/stem-cell_cluster_differential_peaks_between_MS177_and_control_in_the_organoid_doublet_removed_01082025.csv')
View(differential_peaks_2[grepl("chr5-13249", rownames(differential_peaks_2)), ]) #chr5-132495686-132498505 p_val_adj=3.986302e-05

differential_peaks_3 <- FindMarkers(
  object = MS177_organoid_3_conditions,
  ident.1 = "MS177_24hr_epithelial_like",      # Group1, e.g., a cell type or condition
  ident.2 = "control_24hr_epithelial_like",      # Group2, e.g., another cell type or condition
  assay = "peaks",          # Ensure you're analyzing the ATAC assay
  test.use = "LR",         # Use logistic regression, often suitable for sparse data
  min.pct = 0.1,           # Minimum percent of cells with the peak accessible
  logfc.threshold = 0.25   # Set log fold-change threshold for differential accessibility
)
write.csv(differential_peaks_3, '/home/yue1118/TNBC_scRNA_analysis_workflow/TNBC_paper_data_and_figure/epithelial-cell_cluster_differential_peaks_between_MS177_and_control_in_the_organoid_doublet_removed_01082025.csv')
View(differential_peaks_3[grepl("chr5-13249", rownames(differential_peaks_3)), ]) #chr5-132495686-132498505 p_val_adj=1

#IRF1
DefaultAssay(MS177_organoid_3_conditions) <- "peaks"
CoveragePlot(
  object = MS177_organoid_3_conditions,
  region = "chr5-132495686-132498505", # adjust to your region of interest
  group.by = "long_condition"
)
ggsave("/home/yue1118/TNBC_scRNA_analysis_workflow/TNBC_paper_data_and_figure/IRF1_chromatin_accessibility_in_promoter_region_in_3_conditions_split_by_2_clusters_doublet_removed_01082025.svg", width = 8, height = 5)


###To link differential peaks from the MS177 condition with genes using ChIPseeker 
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

#For the organoid as a whole after MS177 treatment
differential_peaks <-  read.csv('/home/yue1118/TNBC_scRNA_analysis_workflow/TNBC_paper_data_and_figure/differential_peaks_between_MS177_and_control_in_the_organoid_doublet_removed_01082025.csv')
differential_peaks <- differential_peaks[differential_peaks$p_val_adj < 0.05,]
differential_peaks$peaks <- differential_peaks$X
peaks_split <- do.call(rbind, strsplit(differential_peaks$X, "[-]"))
colnames(peaks_split) <- c("chr", "start", "end")
peaks_df <- data.frame( chr = peaks_split[, 1],start = as.numeric(peaks_split[, 2]),end = as.numeric(peaks_split[, 3]))
granges_obj <- GRanges(seqnames = peaks_df$chr,ranges = IRanges(start = peaks_df$start, end = peaks_df$end))
peakAnno_MS177 <- annotatePeak(granges_obj, tssRegion = c(-3000, 3000), TxDb = txdb, annoDb = "org.Hs.eg.db") 
df_peakAnno_MS177 <- as.data.frame(peakAnno_MS177)
df_peakAnno_MS177 <- df_peakAnno_MS177 %>% mutate(peaks = paste(seqnames, start, end, sep = "-"))
df_peak_merged <- merge(df_peakAnno_MS177, differential_peaks, by="peaks", all.y=TRUE)
write.csv(df_peak_merged, "/home/yue1118/TNBC_scRNA_analysis_workflow/TNBC_paper_data_and_figure/peaks_annotated_with_genes_after_MS177_doublet_removed_01082025.csv") #Supplementary Table 3

#For the epithelial-like cluster after MS177 treatment
differential_peaks_2 <-  read.csv('/home/yue1118/TNBC_scRNA_analysis_workflow/TNBC_paper_data_and_figure/epithelial-cell_cluster_differential_peaks_between_MS177_and_control_in_the_organoid_doublet_removed_01082025.csv')
differential_peaks_2 <- differential_peaks_2[differential_peaks_2$p_val_adj < 0.05,]
differential_peaks_2$peaks <- differential_peaks_2$X
peaks_split <- do.call(rbind, strsplit(differential_peaks_2$X, "[-]"))
colnames(peaks_split) <- c("chr", "start", "end")
peaks_df <- data.frame( chr = peaks_split[, 1],start = as.numeric(peaks_split[, 2]),end = as.numeric(peaks_split[, 3]))
granges_obj <- GRanges(seqnames = peaks_df$chr,ranges = IRanges(start = peaks_df$start, end = peaks_df$end))
peakAnno_epithelial <- annotatePeak(granges_obj, tssRegion = c(-3000, 3000), TxDb = txdb, annoDb = "org.Hs.eg.db") 
df_peakAnno_epithelial <- as.data.frame(peakAnno_epithelial)
df_peakAnno_epithelial <- df_peakAnno_epithelial %>% mutate(peaks = paste(seqnames, start, end, sep = "-"))
df_peak_merged <- merge(df_peakAnno_epithelial, differential_peaks_2, by="peaks", all.y=TRUE)
write.csv(df_peak_merged, "/home/yue1118/TNBC_scRNA_analysis_workflow/TNBC_paper_data_and_figure/epithelial_cluster_peaks_annotated_with_genes_after_MS177_doublet_removed_01082025.csv") #Supplementary Table 4

#For the stem-cell-like cluster after MS177 treatment
differential_peaks_3 <-  read.csv('/home/yue1118/TNBC_scRNA_analysis_workflow/TNBC_paper_data_and_figure/stem-cell_cluster_differential_peaks_between_MS177_and_control_in_the_organoid_doublet_removed_01082025.csv') 
differential_peaks_3 <- differential_peaks_3[differential_peaks_3$p_val_adj < 0.05,]
differential_peaks_3$peaks <- differential_peaks_3$X
peaks_split <- do.call(rbind, strsplit(differential_peaks_3$X, "[-]"))
colnames(peaks_split) <- c("chr", "start", "end")
peaks_df <- data.frame( chr = peaks_split[, 1],start = as.numeric(peaks_split[, 2]),end = as.numeric(peaks_split[, 3]))
granges_obj <- GRanges(seqnames = peaks_df$chr,ranges = IRanges(start = peaks_df$start, end = peaks_df$end))
peakAnno_stem <- annotatePeak(granges_obj, tssRegion = c(-3000, 3000), TxDb = txdb, annoDb = "org.Hs.eg.db") 
df_peakAnno_stem <- as.data.frame(peakAnno_stem)
df_peakAnno_stem <- df_peakAnno_stem %>% mutate(peaks = paste(seqnames, start, end, sep = "-"))
df_peak_merged <- merge(df_peakAnno_stem, differential_peaks_3, by="peaks", all.y=TRUE)
write.csv(df_peak_merged, "/home/yue1118/TNBC_scRNA_analysis_workflow/TNBC_paper_data_and_figure/stem-cell_cluster_peaks_annotated_with_genes_after_MS177_doublet_removed_01082025.csv") #Supplementary Table 5


###To link differential ATAC peaks with genes after JQ1 treatment

#To find differential ATAC peaks first
Idents(MS177_organoid_3_conditions) <- "dataset"
differential_peaks <- FindMarkers(
  object = MS177_organoid_3_conditions,
  ident.1 = "JQ1_24hr",      # Group1, e.g., a cell type or condition
  ident.2 = "control_24hr",      # Group2, e.g., another cell type or condition
  assay = "peaks",          # Ensure you're analyzing the ATAC assay
  test.use = "LR",         # Use logistic regression, often suitable for sparse data
  min.pct = 0.1,           # Minimum percent of cells with the peak accessible
  logfc.threshold = 0.25   # Set log fold-change threshold for differential accessibility
)
write.csv(differential_peaks, '/home/yue1118/TNBC_scRNA_analysis_workflow/TNBC_paper_data_and_figure/differential_peaks_between_JQ1_and_control_in_the_organoid_doublet_removed_01082025.csv')

Idents(MS177_organoid_3_conditions) <- "long_condition"
DefaultAssay(MS177_organoid_3_conditions) <- "peaks"
differential_peaks_2 <- FindMarkers(
  object = MS177_organoid_3_conditions,
  ident.1 = "JQ1_24hr_stem-cell_like",      # Group1, e.g., a cell type or condition
  ident.2 = "control_24hr_stem-cell_like",      # Group2, e.g., another cell type or condition
  assay = "peaks",          # Ensure you're analyzing the ATAC assay
  test.use = "LR",         # Use logistic regression, often suitable for sparse data
  min.pct = 0.1,           # Minimum percent of cells with the peak accessible
  logfc.threshold = 0.25   # Set log fold-change threshold for differential accessibility
)
write.csv(differential_peaks_2, '/home/yue1118/TNBC_scRNA_analysis_workflow/TNBC_paper_data_and_figure/stem-cell_cluster_differential_peaks_between_JQ1_and_control_in_the_organoid_doublet_removed_01082025.csv')

differential_peaks_3 <- FindMarkers(
  object = MS177_organoid_3_conditions,
  ident.1 = "JQ1_24hr_epithelial_like",      # Group1, e.g., a cell type or condition
  ident.2 = "control_24hr_epithelial_like",      # Group2, e.g., another cell type or condition
  assay = "peaks",          # Ensure you're analyzing the ATAC assay
  test.use = "LR",         # Use logistic regression, often suitable for sparse data
  min.pct = 0.1,           # Minimum percent of cells with the peak accessible
  logfc.threshold = 0.25   # Set log fold-change threshold for differential accessibility
)
write.csv(differential_peaks_3, '/home/yue1118/TNBC_scRNA_analysis_workflow/TNBC_paper_data_and_figure/epithelial-cell_cluster_differential_peaks_between_JQ1_and_control_in_the_organoid_doublet_removed_01082025.csv')

#For the organoid as a whole after JQ1 treatment
differential_peaks <-  read.csv('/home/yue1118/TNBC_scRNA_analysis_workflow/TNBC_paper_data_and_figure/differential_peaks_between_JQ1_and_control_in_the_organoid_doublet_removed_01082025.csv')
differential_peaks <- differential_peaks[differential_peaks$p_val_adj < 0.05,]
differential_peaks$peaks <- differential_peaks$X
peaks_split <- do.call(rbind, strsplit(differential_peaks$X, "[-]"))
colnames(peaks_split) <- c("chr", "start", "end")
peaks_df <- data.frame(chr = peaks_split[, 1],start = as.numeric(peaks_split[, 2]),end = as.numeric(peaks_split[, 3]))
granges_obj <- GRanges(seqnames = peaks_df$chr,ranges = IRanges(start = peaks_df$start, end = peaks_df$end))
peakAnno_JQ1 <- annotatePeak(granges_obj, tssRegion = c(-3000, 3000), TxDb = txdb, annoDb = "org.Hs.eg.db") 
df_peakAnno_JQ1 <- as.data.frame(peakAnno_JQ1)
df_peakAnno_JQ1 <- df_peakAnno_JQ1 %>% mutate(peaks = paste(seqnames, start, end, sep = "-"))
df_peak_merged <- merge(df_peakAnno_JQ1, differential_peaks, by="peaks", all.y=TRUE)
write.csv(df_peak_merged, "/home/yue1118/TNBC_scRNA_analysis_workflow/TNBC_paper_data_and_figure/peaks_annotated_with_genes_after_JQ1_doublet_removed_01082025.csv") #Supplementary Table 7


#For the epithelial-like cluster after JQ1 treatment
differential_peaks_2 <-  read.csv('/home/yue1118/TNBC_scRNA_analysis_workflow/TNBC_paper_data_and_figure/epithelial-cell_cluster_differential_peaks_between_JQ1_and_control_in_the_organoid_doublet_removed_01082025.csv')
differential_peaks_2 <- differential_peaks_2[differential_peaks_2$p_val_adj < 0.05,]
differential_peaks_2$peaks <- differential_peaks_2$X
peaks_split <- do.call(rbind, strsplit(differential_peaks_2$X, "[-]"))
colnames(peaks_split) <- c("chr", "start", "end")
peaks_df <- data.frame(chr = peaks_split[, 1],start = as.numeric(peaks_split[, 2]),end = as.numeric(peaks_split[, 3]))
granges_obj <- GRanges(seqnames = peaks_df$chr,ranges = IRanges(start = peaks_df$start, end = peaks_df$end))
peakAnno_epithelial <- annotatePeak(granges_obj, tssRegion = c(-3000, 3000), TxDb = txdb, annoDb = "org.Hs.eg.db") 
df_peakAnno_epithelial <- as.data.frame(peakAnno_epithelial)
df_peakAnno_epithelial <- df_peakAnno_epithelial %>% mutate(peaks = paste(seqnames, start, end, sep = "-"))
df_peak_merged <- merge(df_peakAnno_epithelial, differential_peaks_2, by="peaks", all.y=TRUE)
write.csv(df_peak_merged, "/home/yue1118/TNBC_scRNA_analysis_workflow/TNBC_paper_data_and_figure/epithelial_cluster_peaks_annotated_with_genes_after_JQ1_doublet_removed_01082025.csv") #Supplementary Table 8

#For the stem-cell-like cluster after JQ1 treatment
differential_peaks_3 <-  read.csv('/home/yue1118/TNBC_scRNA_analysis_workflow/TNBC_paper_data_and_figure/stem-cell_cluster_differential_peaks_between_JQ1_and_control_in_the_organoid_doublet_removed_01082025.csv') 
differential_peaks_3 <- differential_peaks_3[differential_peaks_3$p_val_adj < 0.05,]
differential_peaks_3$peaks <- differential_peaks_3$X
peaks_split <- do.call(rbind, strsplit(differential_peaks_3$X, "[-]"))
colnames(peaks_split) <- c("chr", "start", "end")
peaks_df <- data.frame(chr = peaks_split[, 1],start = as.numeric(peaks_split[, 2]),end = as.numeric(peaks_split[, 3]))
granges_obj <- GRanges(seqnames = peaks_df$chr,ranges = IRanges(start = peaks_df$start, end = peaks_df$end))
peakAnno_stem <- annotatePeak(granges_obj, tssRegion = c(-3000, 3000), TxDb = txdb, annoDb = "org.Hs.eg.db") 
df_peakAnno_stem <- as.data.frame(peakAnno_stem)
df_peakAnno_stem <- df_peakAnno_stem %>% mutate(peaks = paste(seqnames, start, end, sep = "-"))
df_peak_merged <- merge(df_peakAnno_stem, differential_peaks_3, by="peaks", all.y=TRUE)
write.csv(df_peak_merged, "/home/yue1118/TNBC_scRNA_analysis_workflow/TNBC_paper_data_and_figure/stem-cell_cluster_peaks_annotated_with_genes_after_JQ1_doublet_removed_01082025.csv") #Supplementary Table 9

df1 <- read.csv("/home/yue1118/TNBC_scRNA_analysis_workflow/TNBC_paper_data_and_figure/stem-cell_cluster_peaks_annotated_with_genes_after_JQ1_doublet_removed_01082025.csv")
View(df1[grepl("chr5-13249", df1$peaks), ])

df2 <- read.csv("/home/yue1118/TNBC_scRNA_analysis_workflow/TNBC_paper_data_and_figure/stem-cell_cluster_peaks_annotated_with_genes_after_MS177_doublet_removed_01082025.csv")
View(df2[grepl("chr5-13249", df2$peaks), ])

df3 <- read.csv("/home/yue1118/TNBC_scRNA_analysis_workflow/TNBC_paper_data_and_figure/peaks_annotated_with_genes_after_JQ1_doublet_removed_01082025.csv")
View(df3[grepl("chr5-13249", df3$peaks), ])


df4 <- read.csv("/home/yue1118/TNBC_scRNA_analysis_workflow/TNBC_paper_data_and_figure/epithelial_cluster_peaks_annotated_with_genes_after_JQ1_doublet_removed_01082025.csv")
View(df4[grepl("chr5-13249", df4$peaks), ])
