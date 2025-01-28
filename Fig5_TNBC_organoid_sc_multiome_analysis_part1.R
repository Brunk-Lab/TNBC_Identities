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
library(DoubletFinder)

#Codes from here to line 386 are for quality control and building the organoid Seurat object.

# get gene annotations for hg38
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevels(annotation) <- paste0('chr', seqlevels(annotation))

#STEP1 ---to merge scATAC datasets from 3 treatment conditions This is derived from https://stuartlab.org/signac/articles/merging
# read in peak sets

organoid_atac_control_24hr <- read.table(
  file = "/datastore/lbcfs/labs/brunk_lab/private/TNBC/cell_ranger_arc_results/TNBC_organoid_control_24hr/Control_24hr/outs/atac_peaks.bed",
  col.names = c("chr", "start", "end")
)


organoid_atac_JQ1_24hr <- read.table(
  file = "/datastore/lbcfs/labs/brunk_lab/private/TNBC/cell_ranger_arc_results/TNBC_organoid_JQ1_24hr/JQ1_24hr/outs/atac_peaks.bed",
  col.names = c("chr", "start", "end")
)


organoid_atac_MS177_24hr <- read.table(
  file = "/datastore/lbcfs/labs/brunk_lab/private/TNBC/cell_ranger_arc_results/TNBC_organoid_MSL77_24hr/MSL77_24hr/outs/atac_peaks.bed",
  col.names = c("chr", "start", "end")
)

organoid_atac_control_24hr <- makeGRangesFromDataFrame(organoid_atac_control_24hr)
organoid_atac_JQ1_24hr <- makeGRangesFromDataFrame(organoid_atac_JQ1_24hr)
organoid_atac_MS177_24hr <- makeGRangesFromDataFrame(organoid_atac_MS177_24hr)


# Create a unified set of peaks to quantify in each dataset
combined.peaks <- reduce(x = c(organoid_atac_control_24hr, organoid_atac_JQ1_24hr, organoid_atac_MS177_24hr))

# Filter out bad peaks based on length
peakwidths <- width(combined.peaks)
combined.peaks <- combined.peaks[peakwidths  < 10000 & peakwidths > 20]


# load metadata

md.control_24hr <- read.table(
  file = "/datastore/lbcfs/labs/brunk_lab/private/TNBC/cell_ranger_arc_results/TNBC_organoid_control_24hr/Control_24hr/outs/per_barcode_metrics.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)


md.JQ1_24hr <- read.table(
  file = "/datastore/lbcfs/labs/brunk_lab/private/TNBC/cell_ranger_arc_results/TNBC_organoid_JQ1_24hr/JQ1_24hr/outs/per_barcode_metrics.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)

md.MS177_24hr <- read.table(
  file = "/datastore/lbcfs/labs/brunk_lab/private/TNBC/cell_ranger_arc_results/TNBC_organoid_MSL77_24hr/MSL77_24hr/outs/per_barcode_metrics.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)

# perform an initial filtering of low count cells

md.control_24hr <- md.control_24hr[md.control_24hr$atac_fragments > 500, ]
md.JQ1_24hr <- md.JQ1_24hr[md.JQ1_24hr$atac_fragments > 500, ]
md.MS177_24hr <- md.MS177_24hr[md.MS177_24hr$atac_fragments > 500, ]

# create fragment objects

frags_control_24hr <- CreateFragmentObject(
  path = "/datastore/lbcfs/labs/brunk_lab/private/TNBC/cell_ranger_arc_results/TNBC_organoid_control_24hr/Control_24hr/outs/atac_fragments.tsv.gz",
  cells = rownames(md.control_24hr)
)


frags_JQ1_24hr <- CreateFragmentObject(
  path = "/datastore/lbcfs/labs/brunk_lab/private/TNBC/cell_ranger_arc_results/TNBC_organoid_JQ1_24hr/JQ1_24hr/outs/atac_fragments.tsv.gz",
  cells = rownames(md.JQ1_24hr)
)


frags_MS177_24hr <- CreateFragmentObject(
  path = "/datastore/lbcfs/labs/brunk_lab/private/TNBC/cell_ranger_arc_results/TNBC_organoid_MSL77_24hr/MSL77_24hr/outs/atac_fragments.tsv.gz",
  cells = rownames(md.MS177_24hr)
)

#quantify peaks in each dataset

control_24hr.counts <- FeatureMatrix(
  fragments = frags_control_24hr,
  features = combined.peaks,
  cells = rownames(md.control_24hr)
)


JQ1_24hr.counts <- FeatureMatrix(
  fragments = frags_JQ1_24hr,
  features = combined.peaks,
  cells = rownames(md.JQ1_24hr)
)


MS177_24hr.counts <- FeatureMatrix(
  fragments = frags_MS177_24hr,
  features = combined.peaks,
  cells = rownames(md.MS177_24hr)
)


control_24hr_assay <- CreateChromatinAssay(control_24hr.counts, fragments = frags_control_24hr, annotation = annotation)
control_24hr <- CreateSeuratObject(control_24hr_assay, assay = "ATAC", meta.data=md.control_24hr, annotation = annotation) #the annotation parameter was likely not utilized or stored by Seurat during the execution of CreateSeuratObject because the codes ran through successfully without any issue.


JQ1_24hr_assay <- CreateChromatinAssay(JQ1_24hr.counts, fragments = frags_JQ1_24hr, annotation = annotation)
JQ1_24hr <- CreateSeuratObject(JQ1_24hr_assay, assay = "ATAC", meta.data=md.JQ1_24hr, annotation = annotation)


MS177_24hr_assay <- CreateChromatinAssay(MS177_24hr.counts, fragments = frags_MS177_24hr, annotation = annotation)
MS177_24hr <- CreateSeuratObject(MS177_24hr_assay, assay = "ATAC", meta.data=md.MS177_24hr, annotation = annotation)



control_24hr$dataset <- 'control_24hr'
JQ1_24hr$dataset <- 'JQ1_24hr'
MS177_24hr$dataset <- 'MS177_24hr'


# merge all datasets, adding a cell ID to make sure cell names are unique
atac_combined <- merge(
  x = control_24hr,
  y = list(JQ1_24hr, MS177_24hr),
  add.cell.ids = c("control_24hr", "JQ1_24hr", "MS177_24hr")
)

atac_combined[["ATAC"]]


#STEP2 --- to merge scRNA data

control_24hr.multiome <- Read10X_h5("/datastore/lbcfs/labs/brunk_lab/private/TNBC/cell_ranger_arc_results/TNBC_organoid_control_24hr/Control_24hr/outs/filtered_feature_bc_matrix.h5")
JQ1_24hr.multiome <- Read10X_h5("/datastore/lbcfs/labs/brunk_lab/private/TNBC/cell_ranger_arc_results/TNBC_organoid_JQ1_24hr/JQ1_24hr/outs/filtered_feature_bc_matrix.h5")
MS177_24hr.multiome <- Read10X_h5("/datastore/lbcfs/labs/brunk_lab/private/TNBC/cell_ranger_arc_results/TNBC_organoid_MSL77_24hr/MSL77_24hr/outs/filtered_feature_bc_matrix.h5")


control_24hr_rna_counts <- control_24hr.multiome$`Gene Expression`
JQ1_24hr_rna_counts <- JQ1_24hr.multiome$`Gene Expression`
MS177_24hr_rna_counts <- MS177_24hr.multiome$`Gene Expression`


control_24hr_rna <- CreateSeuratObject(counts = control_24hr_rna_counts, project = "TNBC_organoid")
JQ1_24hr_rna <- CreateSeuratObject(counts = JQ1_24hr_rna_counts, project = "TNBC_organoid")
MS177_24hr_rna <- CreateSeuratObject(counts = MS177_24hr_rna_counts, project = "TNBC_organoid")

rna.combined <- merge(control_24hr_rna, y = list(JQ1_24hr_rna,  MS177_24hr_rna), add.cell.ids = c("control_24hr",  "JQ1_24hr","MS177_24hr"), project = "TNBC_organoid")

atac_combined_2 <- atac_combined[,colnames(atac_combined) %in% colnames(rna.combined)]
rna.combined_2 <- rna.combined[, colnames(rna.combined) %in% colnames(atac_combined_2)] #It is checked that atac_combined_2 and rna.combined_2 have column names(cell names) in the same order.
atac_combined_2[['RNA']] <- CreateAssayObject(counts = GetAssayData(rna.combined_2[['RNA']]))
DefaultAssay(atac_combined_2) <- 'RNA'
atac_combined_2[["percent.mt"]] <- PercentageFeatureSet(atac_combined_2, pattern = "^MT-")

DefaultAssay(atac_combined_2) <- "ATAC"

atac_combined_2 <- NucleosomeSignal(atac_combined_2)
atac_combined_2 <- TSSEnrichment(atac_combined_2)

VlnPlot(atac_combined_2, features = c("nCount_ATAC", "nCount_RNA", "nFeature_RNA", "TSS.enrichment", "nucleosome_signal", "percent.mt"), ncol = 3,
        log = TRUE, pt.size = 0) + NoLegend()

MS177_organoid_3_conditions <- subset(
  x = atac_combined_2,
  subset = nCount_ATAC < 7e4 &
    nCount_ATAC > 1000 &
    nCount_RNA < 25000 &
    nCount_RNA > 1000 &
    nFeature_RNA > 500  &
    nucleosome_signal < 2 &
    TSS.enrichment > 1 &
    percent.mt < 20
)

MS177_organoid_3_conditions

DefaultAssay(MS177_organoid_3_conditions) <- "ATAC"
atac_combined_peaks <- CallPeaks(MS177_organoid_3_conditions)
atac_combined_peaks <- keepStandardChromosomes(atac_combined_peaks, pruning.mode = "coarse")
atac_combined_peaks <- subsetByOverlaps(x = atac_combined_peaks, ranges = blacklist_hg38_unified, invert = TRUE)

DefaultAssay(MS177_organoid_3_conditions) <- "ATAC"
# quantify counts in each peak
organoid_macs2_counts <- FeatureMatrix(
  fragments = Fragments(MS177_organoid_3_conditions),
  features = atac_combined_peaks,
  cells = colnames(MS177_organoid_3_conditions)
)

# create a new assay using the MACS2 peak set and add it to the Seurat object
MS177_organoid_3_conditions[["peaks"]] <- CreateChromatinAssay(
  counts = organoid_macs2_counts,
  fragments = Fragments(MS177_organoid_3_conditions),
  annotation = annotation
)

# ATAC analysis
# We exclude the first dimension as this is typically correlated with sequencing depth
DefaultAssay(MS177_organoid_3_conditions) <- "peaks"
MS177_organoid_3_conditions <- FindTopFeatures(MS177_organoid_3_conditions, min.cutoff = 'q0')
MS177_organoid_3_conditions <- RunTFIDF(MS177_organoid_3_conditions)
MS177_organoid_3_conditions <- RunSVD(MS177_organoid_3_conditions)
MS177_organoid_3_conditions <- RunUMAP(MS177_organoid_3_conditions, reduction = 'lsi', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_", seed.use = 7)

# RNA analysis
DefaultAssay(MS177_organoid_3_conditions) <- "RNA"
MS177_organoid_3_conditions <- SCTransform(MS177_organoid_3_conditions, verbose = FALSE) %>% RunPCA() %>% RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_', seed.use = 7)

MS177_organoid_3_conditions <- FindMultiModalNeighbors(MS177_organoid_3_conditions, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50))
MS177_organoid_3_conditions <- RunUMAP(MS177_organoid_3_conditions, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_", seed.use = 7)
MS177_organoid_3_conditions <- FindClusters(MS177_organoid_3_conditions, graph.name = "wsnn", algorithm = 3, verbose = FALSE)
MS177_organoid_3_conditions$wnn_cluster_identity <- Idents(MS177_organoid_3_conditions)

DefaultAssay(MS177_organoid_3_conditions) <- "ATAC"
main.chroms <- standardChromosomes(BSgenome.Hsapiens.UCSC.hg38)
keep.peaks <- which(as.character(seqnames(granges(MS177_organoid_3_conditions))) %in% main.chroms)
MS177_organoid_3_conditions[["ATAC"]] <- subset(MS177_organoid_3_conditions[["ATAC"]], features = rownames(MS177_organoid_3_conditions[["ATAC"]])[keep.peaks])

DefaultAssay(MS177_organoid_3_conditions) <- "ATAC"
pwm_set <- getMatrixSet(x = JASPAR2020, opts = list(species = 9606, all_versions = FALSE))
motif.matrix <- CreateMotifMatrix(features = granges(MS177_organoid_3_conditions), pwm = pwm_set, genome = 'hg38', use.counts = FALSE)
motif.object <- CreateMotifObject(data = motif.matrix, pwm = pwm_set)
MS177_organoid_3_conditions <- SetAssayData(MS177_organoid_3_conditions, assay = 'ATAC', slot = 'motifs', new.data = motif.object)

# Note that this step can take 30-60 minutes 
MS177_organoid_3_conditions<- RunChromVAR(
  object = MS177_organoid_3_conditions,
  genome = BSgenome.Hsapiens.UCSC.hg38
)

#To compute gene activity for each gene based on their chromatin accessibility
DefaultAssay(MS177_organoid_3_conditions) <- "peaks"
organoid.gene.activities <- GeneActivity(MS177_organoid_3_conditions, biotypes = NULL) # to include not only protein coding genes but also long non coding genes
MS177_organoid_3_conditions[['gene_activity']] <- CreateAssayObject(counts = organoid.gene.activities)
MS177_organoid_3_conditions <- NormalizeData(
  object = MS177_organoid_3_conditions,
  assay = 'gene_activity',
  normalization.method = 'LogNormalize',
  scale.factor = median(MS177_organoid_3_conditions$nCount_RNA)
)

##To differentiate epithelial cluster from stem-cell like cluster
DefaultAssay(MS177_organoid_3_conditions) <- "SCT"
DimPlot(MS177_organoid_3_conditions, label = T, reduction = "wnn.umap")
Idents(MS177_organoid_3_conditions) <- MS177_organoid_3_conditions$wnn_cluster_identity
FeaturePlot(MS177_organoid_3_conditions, features = c("CD24", "CD44", "VIM", "FN1"), label = TRUE, reduction = "wnn.umap")#epithelial-like cluster: 0, 1, 5, 8, 9, 11, 12, 13, 14, 15
#stem-cell-like cluster: 2, 3, 4, 6, 7, 10
Idents(MS177_organoid_3_conditions) <- MS177_organoid_3_conditions$dataset
DimPlot(MS177_organoid_3_conditions, reduction = "wnn.umap", group.by = "dataset", label = FALSE, label.size = 2.5, repel = TRUE)


func_epithelial_or_stem <- function(x){
  if (x %in% c(0, 1, 5, 8, 9, 11, 12, 13, 14, 15)) {
    "epithelial_like"
  } else {
    "stem-cell_like"
  }
}

df_wnn_cluster_identity <- as.data.frame(MS177_organoid_3_conditions@meta.data[, c("wnn_cluster_identity")])
rownames(df_wnn_cluster_identity) <- colnames(MS177_organoid_3_conditions)
View(df_wnn_cluster_identity)
colnames(df_wnn_cluster_identity) <- c("wnn_cluster_identity")

df_wnn_cluster_identity$epithelial_or_stem <- lapply(df_wnn_cluster_identity$wnn_cluster_identity, func_epithelial_or_stem)
df_wnn_cluster_identity$epithelial_or_stem <- as.character(df_wnn_cluster_identity$epithelial_or_stem)
MS177_organoid_3_conditions$epithelial_or_stem <- df_wnn_cluster_identity$epithelial_or_stem
table(MS177_organoid_3_conditions$epithelial_or_stem)
MS177_organoid_3_conditions$long_condition <- paste(MS177_organoid_3_conditions$dataset, MS177_organoid_3_conditions$epithelial_or_stem, sep = "_")

table(MS177_organoid_3_conditions$long_condition)

#saveRDS(MS177_organoid_3_conditions, './Organoid_TNBC_JQ1_MS177_24hr_3_conditions_with_gene_activity_with_cluster_identities_09272024.rds')
MS177_organoid_3_conditions <- readRDS('/home/yue1118/TNBC_scRNA_analysis_workflow/Organoid_TNBC_JQ1_MS177_24hr_3_conditions_with_gene_activity_with_cluster_identities_09272024.rds')

###Now run DoubletFinder to remove potential doublets from the dataset

suppressMessages(require(DoubletFinder))
Idents(MS177_organoid_3_conditions) <- "dataset"
#To run doubletfinder on different conditions separately
#For the control condition
MS177_organoid_3_conditions_subset_control <- subset(MS177_organoid_3_conditions, idents = "control_24hr")
#To look for the best pK value
sweep.res.list_control <- paramSweep(MS177_organoid_3_conditions_subset_control, PCs = 1:20, sct =TRUE)
sweep.stats_control <- summarizeSweep(sweep.res.list_control, GT = FALSE)
bcmvn_control <- find.pK(sweep.stats_control)
pK <- bcmvn_control %>% filter(BCmetric == max(BCmetric)) %>% select(pK) #pK=0.17
#Homotypic Doublet Proportion Estimate
cluster_annotation <- MS177_organoid_3_conditions_subset_control@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(cluster_annotation)
#To estimate the number of doublets
nExp_poi <- round(0.05*nrow(MS177_organoid_3_conditions_subset_control@meta.data)) #assuming there are 5% doublets
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

MS177_organoid_3_conditions_subset_control <- DoubletFinder::doubletFinder(MS177_organoid_3_conditions_subset_control, pN = 0.25, pK = 0.17, nExp = nExp_poi.adj, PCs = 1:10, sct = TRUE)
df_control_doublet <- MS177_organoid_3_conditions_subset_control@meta.data[,c("pANN_0.25_0.17_272", "DF.classifications_0.25_0.17_272")]
colnames(df_control_doublet) <- c("pANN", "DF.classification")

#For the MS177 condition
MS177_organoid_3_conditions_subset_MS177 <- subset(MS177_organoid_3_conditions, idents = "MS177_24hr")
#To look for the best pK value
sweep.res.list_MS177 <- paramSweep(MS177_organoid_3_conditions_subset_MS177, PCs = 1:20, sct =TRUE)
sweep.stats_MS177 <- summarizeSweep(sweep.res.list_MS177, GT = FALSE)
bcmvn_MS177 <- find.pK(sweep.stats_MS177)
pK <- bcmvn_MS177 %>% filter(BCmetric == max(BCmetric)) %>% select(pK) #pK=0.07
#Homotypic Doublet Proportion Estimate
cluster_annotation <- MS177_organoid_3_conditions_subset_MS177@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(cluster_annotation)
#To estimate the number of doublets
nExp_poi <- round(0.05*nrow(MS177_organoid_3_conditions_subset_MS177@meta.data)) #assuming there are 5% doublets
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
MS177_organoid_3_conditions_subset_MS177 <- DoubletFinder::doubletFinder(MS177_organoid_3_conditions_subset_MS177, pN = 0.25, pK = 0.07 , nExp = nExp_poi.adj, PCs = 1:10, sct = TRUE)
df_MS177_doublet <- MS177_organoid_3_conditions_subset_MS177@meta.data[,c("pANN_0.25_0.07_203", "DF.classifications_0.25_0.07_203")]
colnames(df_MS177_doublet) <- c("pANN", "DF.classification")

#For the JQ1 condition
MS177_organoid_3_conditions_subset_JQ1 <- subset(MS177_organoid_3_conditions, idents = "JQ1_24hr")
#To look for the best pK value
sweep.res.list_JQ1 <- paramSweep(MS177_organoid_3_conditions_subset_JQ1, PCs = 1:20, sct =TRUE)
sweep.stats_JQ1 <- summarizeSweep(sweep.res.list_JQ1, GT = FALSE)
bcmvn_JQ1 <- find.pK(sweep.stats_JQ1)
pK <- bcmvn_JQ1 %>% filter(BCmetric == max(BCmetric)) %>% select(pK) #pK=0.09
#Homotypic Doublet Proportion Estimate
cluster_annotation <- MS177_organoid_3_conditions_subset_JQ1@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(cluster_annotation)
#To estimate the number of doublets
nExp_poi <- round(0.05*nrow(MS177_organoid_3_conditions_subset_JQ1@meta.data)) #assuming there are 5% doublets
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
MS177_organoid_3_conditions_subset_JQ1 <- DoubletFinder::doubletFinder(MS177_organoid_3_conditions_subset_JQ1, pN = 0.25, pK = 0.09 , nExp = nExp_poi.adj, PCs = 1:10, sct = TRUE)
df_JQ1_doublet <- MS177_organoid_3_conditions_subset_JQ1@meta.data[,c("pANN_0.25_0.09_170", "DF.classifications_0.25_0.09_170")]
colnames(df_JQ1_doublet) <- c("pANN", "DF.classification")

df_doublet <- rbind(df_control_doublet, df_JQ1_doublet, df_MS177_doublet)
df1 <- as.data.frame(rownames(df_doublet) == rownames(MS177_organoid_3_conditions@meta.data))
table(df1$`rownames(df_doublet) == rownames(MS177_organoid_3_conditions@meta.data)`)

MS177_organoid_3_conditions$DF.classification <- df_doublet$DF.classification

saveRDS(MS177_organoid_3_conditions, '/home/yue1118/TNBC_scRNA_analysis_workflow/Organoid_TNBC_JQ1_MS177_24hr_3_conditions_with_gene_activity_with_cluster_identities_with_doublets_classification_01082025.rds')
MS177_organoid_3_conditions <- readRDS('/home/yue1118/TNBC_scRNA_analysis_workflow/Organoid_TNBC_JQ1_MS177_24hr_3_conditions_with_gene_activity_with_cluster_identities_with_doublets_classification_01082025.rds')

Idents(MS177_organoid_3_conditions) <- "DF.classification"
MS177_organoid_3_conditions <- subset(MS177_organoid_3_conditions, idents = "Singlet")

dim(MS177_organoid_3_conditions)
#For Fig 5A
DimPlot(MS177_organoid_3_conditions, reduction = "wnn.umap", group.by = "dataset", label = FALSE, label.size = 2.5, repel = TRUE)      
ggsave('/home/yue1118/TNBC_scRNA_analysis_workflow/TNBC_paper_data_and_figure/DimPlot_tnbc_organoid_treatment_condition_doublets_removed_01082025.svg', dpi = 300, height = 10, width = 10)

custom_colors <- c("epithelial_like" = "#078992", "stem-cell_like" = "#9E6300")
DimPlot(MS177_organoid_3_conditions, reduction = "wnn.umap", group.by = "epithelial_or_stem", label = FALSE, label.size = 2.5, repel = TRUE,cols = custom_colors)      
ggsave('/home/yue1118/TNBC_scRNA_analysis_workflow/TNBC_paper_data_and_figure/DimPlot_tnbc_organoid_epithelial_or_stem_doublets_removed_01082025.svg', dpi = 300, height = 10, width = 10)


#For Fig 5B
table(MS177_organoid_3_conditions$long_condition)
#Then numbers of each cluster are used to make stacked barplot in adobe illustrator


#For Fig 5C
w <- read.delim("/home/yue1118/TNBC_paper_data/CCLE_TNBC_model/NMF_using_524_genes_selected_using_pan_cancer_cell_lines_including_MYC_on_30_updated_TNBC_cell_lines_normalized_for_samples_12222022/7/w.tsv")
w_matrix_genes <- w$X
empty_df <- data.frame(row.names = w$X)
empty_df$ID <- w$X
w_numeric <- as.matrix(sapply(w[, -1], as.numeric))

#to plot pie charts for the organoid as a whole for 3 treatment conditions
DefaultAssay(MS177_organoid_3_conditions) <- "SCT"
Idents(MS177_organoid_3_conditions) <- MS177_organoid_3_conditions$dataset


for (i in c("control_24hr", "JQ1_24hr", "MS177_24hr")) {
  df_state <- data.frame(matrix(ncol = 7, nrow = 0))
  organoid_individual <- subset(MS177_organoid_3_conditions, idents = i)
  sct_data <- GetAssayData(organoid_individual, assay = "SCT", slot = "counts")
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
  rownames(H_matrix) <- colnames(organoid_individual)
  df_state <- rbind(df_state, H_matrix)
  zero_rows <- df_state[apply(df_state, 1, function(x) all(x == 0)), ]
  df_non_zero_rows <- df_state[apply(df_state, 1, function(x) any(x != 0)), ]#0 rows with all 0 values
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
  ggsave(paste0("/home/yue1118/TNBC_scRNA_analysis_workflow/TNBC_paper_data_and_figure/Organoid_", i, "_pie_chart_doublets_removed_01082025.svg"), height = 5, width = 5) #11242024.svg
}

#For Fig 5G
#to plot pie charts for the 2 clusters in the organoid for 3 treatment conditions
DefaultAssay(MS177_organoid_3_conditions) <- "SCT"
Idents(MS177_organoid_3_conditions) <- MS177_organoid_3_conditions$long_condition

for (i in c("control_24hr_epithelial_like", "control_24hr_stem-cell_like", "JQ1_24hr_epithelial_like", "JQ1_24hr_stem-cell_like","MS177_24hr_epithelial_like", "MS177_24hr_stem-cell_like")) {
  df_state <- data.frame(matrix(ncol = 7, nrow = 0))
  organoid_individual <- subset(MS177_organoid_3_conditions, idents = i)
  sct_data <- GetAssayData(organoid_individual, assay = "SCT", slot = "counts")
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
  rownames(H_matrix) <- colnames(organoid_individual)
  df_state <- rbind(df_state, H_matrix)
  zero_rows <- df_state[apply(df_state, 1, function(x) all(x == 0)), ]
  df_non_zero_rows <- df_state[apply(df_state, 1, function(x) any(x != 0)), ]#0 rows with all 0 values
  max_col_index <- apply(df_non_zero_rows, 1, which.max)
  max_col_names <- colnames(df_non_zero_rows)[max_col_index]
  df_max_state <- data.frame(rownames(df_non_zero_rows), max_col_names)
  df1 <- as.data.frame(table(df_max_state$max_col_names))
  df_state_percentage <- data.frame(category = df1$Var1, count = df1$Freq)
  df_state_percentage$percentage <- (df_state_percentage$count / sum(df_state_percentage$count)) *100
  print(df_state_percentage)
  ggplot(df_state_percentage, aes(x = "", y = count, fill = category)) +
    geom_bar(width = 1, stat = "identity", color = "white") +
    coord_polar("y", start = 0) +
    theme_void() +
    scale_fill_manual(values = c("V1" = "lightblue", "V2" = "blue",  "V3" = "chocolate1", "V4" = "darkred","V5"="brown2", "V6" = "darksalmon",  "V7"="aquamarine"))# + 
  #labs(title = paste("single-cell state map"), fill = "Category")
  ggsave(paste0("/home/yue1118/TNBC_scRNA_analysis_workflow/TNBC_paper_data_and_figure/Organoid_", i, "_pie_chart_doublets_removed_01082025.svg"), height = 5, width = 5) ##11242024.svg
}

#For Fig 5D

# Extract cluster information
cluster1 <- "control_24hr_epithelial_like" 
cluster2 <- "control_24hr_stem-cell_like"
cluster3 <- "JQ1_24hr_epithelial_like"
cluster4 <- "JQ1_24hr_stem-cell_like"
cluster5 <- "MS177_24hr_epithelial_like"
cluster6 <- "MS177_24hr_stem-cell_like"

genes_of_interest <- c("SFRP1","SGCD","DST","CALD1","TPM1","MGP","FN1","SLIT3","CD59","SAT1","CD44")
# Fetch data for genes of interest
expression_data <- FetchData(MS177_organoid_3_conditions, vars =genes_of_interest, layer = "SCT")
expression_data$cluster <- MS177_organoid_3_conditions@meta.data$long_condition

# Subset data by clusters
cluster1_data <- expression_data[expression_data$cluster == cluster1, ]
cluster2_data <- expression_data[expression_data$cluster == cluster2, ]
cluster3_data <- expression_data[expression_data$cluster == cluster3, ]
cluster4_data <- expression_data[expression_data$cluster == cluster4, ]
cluster5_data <- expression_data[expression_data$cluster == cluster5, ]
cluster6_data <- expression_data[expression_data$cluster == cluster6, ]

# Calculate mean expression for each cluster
cluster1_mean <- colMeans(cluster1_data[, genes_of_interest])
cluster2_mean <- colMeans(cluster2_data[, genes_of_interest])
cluster3_mean <- colMeans(cluster3_data[, genes_of_interest])
cluster4_mean <- colMeans(cluster4_data[, genes_of_interest])
cluster5_mean <- colMeans(cluster5_data[, genes_of_interest])
cluster6_mean <- colMeans(cluster6_data[, genes_of_interest])

# Combine the results
mean_expression <- data.frame(
  Gene = genes_of_interest,
  Cluster1 = cluster1_mean,
  Cluster2 = cluster2_mean,
  Cluster3 = cluster3_mean,
  Cluster4 = cluster4_mean,
  Cluster5 = cluster5_mean,
  Cluster6 = cluster6_mean
)
colnames(mean_expression) <- c("Gene","control_24hr_epithelial_like","control_24hr_stem-cell_like", "JQ1_24hr_epithelial_like","JQ1_24hr_stem-cell_like","MS177_24hr_epithelial_like","MS177_24hr_stem-cell_like")
mean_expression_heatmap <- mean_expression[, c("control_24hr_stem-cell_like", "JQ1_24hr_stem-cell_like", "MS177_24hr_stem-cell_like", "control_24hr_epithelial_like", "JQ1_24hr_epithelial_like", "MS177_24hr_epithelial_like")]

svg('/home/yue1118/TNBC_scRNA_analysis_workflow/TNBC_paper_data_and_figure/EMT_markers_DEG_heatmap_between_two_clusters_in_the_organoid_doublet_removed_01082025.svg',width = 6, height = 6)
pheatmap(as.matrix(mean_expression_heatmap), color = colorRampPalette(c("white", "red"))(256), cluster_rows = FALSE)
dev.off()

#Fig 5E.

DefaultAssay(MS177_organoid_3_conditions) <- "SCT"
#MsigDB hallmark gene set TNF signaling via NF-kappaB
list_module_NFkB <- c("ABCA1",	"ACKR3",	"AREG",	"ATF3","ATP2B1",	"B4GALT1",	"B4GALT5",	"BCL2A1",	"BCL3",	"BCL6",	"BHLHE40",	"BIRC2",	"BIRC3",	"BMP2",	"BTG1",	"BTG2",	"BTG3",	"CCL2",	"CCL20",	"CCL4","CCL5",
                      "CCN1",	"CCND1",	"CCNL1",	"CCRL2",	"CD44",	"CD69",	"CD80",	"CD83",	"CDKN1A",	"CEBPB",	"CEBPD",	"CFLAR",	"CLCF1",	"CSF1",	"CSF2",	"CXCL1",	"CXCL10",	"CXCL11","CXCL2",	"CXCL3",	"CXCL6",	"DENND5A",	"DNAJB4",	"DRAM1",	"DUSP1",	"DUSP2",	"DUSP4",	"DUSP5",	"EDN1",
                      "EFNA1",	"EGR1",	"EGR2",	"EGR3",	"EHD1",	"EIF1",	"ETS2",	"F2RL1",	"F3",	"FJX1",	"FOS",	"FOSB",	"FOSL1",	"FOSL2",	"FUT4",	"G0S2",	"GADD45A",	"GADD45B",	"GCH1",	"GEM",	"GFPT2","GPR183","HBEGF",	"HES1",	"ICAM1",	"ICOSLG",	"ID2",	"IER2",	"IER3",	"IER5",	"IFIH1",	"IFIT2",	"IFNGR2", "IL12B",	"IL15RA",	"IL18",	"IL1A", "IL1B",	 "IL23A",	"IL6",	"IL6ST",	"IL7R",	
                      "INHBA",	"IRF1",	"IRS2",	"JAG1",	"JUN",	"JUNB",	"KDM6B",	"KLF10",	"KLF2",	"KLF4",	"KLF6", "KLF9", "KYNU",	"LAMB3",	"LDLR", "LIF",	"LITAF", "MAFF",	"MAP2K3",	"MAP3K8",	"MARCKS",	"MCL1", "MSC", "MXD1",	"MYC",	"NAMPT",	"NFAT5", "NFE2L2",	"NFIL3",	"NFKB1",	"NFKB2",	"NFKBIA",	"NFKBIE",	"NINJ1",	"NR4A1", "NR4A2",	"NR4A3",	"OLR1", "PANX1",	"PDE4B", "PDLIM5",	
                      "PER1",	"PFKFB3",	"PHLDA1",	"PHLDA2",	"PLAU",	"PLAUR",	"PLEK",	"PLK2",	"PLPP3",	"PMEPA1",	"PNRC1",	"PPP1R15A", "PTGER4",	"PTGS2",	"PTPRE",	"PTX3",	"RCAN1",	"REL", "RELA",	"RELB",	"RHOB",	"RIGI",	"RIPK2",	"RNF19B",	"SAT1",	"SDC4",	"SERPINB2",	"SERPINB8", "SERPINE1",	"SGK1",	"SIK1",	"SLC16A6",	"SLC2A3",	"SLC2A6",	"SMAD3",	"SNN",	"SOCS3", "SOD2", "SPHK1",	"SPSB1",	"SQSTM1",	"STAT5A",	"TANK",	"TAP1",	"TGIF1",	"TIPARP",	"TLR2",	
                      "TNC",	"TNF",	"TNFAIP2",	"TNFAIP3",	"TNFAIP6",	"TNFAIP8",	"TNFRSF9",	"TNFSF9",	"TNIP1",	"TNIP2", 
                      "TRAF1",	"TRIB1", "TRIP10",	"TSC22D1",	"TUBB2A",	"VEGFA",	"YRDC",	"ZBTB10",	"ZC3H12A",	"ZFP36") 

list_module_NFkB <- list(list_module_NFkB)
MS177_organoid_3_conditions <- AddModuleScore(MS177_organoid_3_conditions,features = list_module_NFkB, pool = NULL,nbin = 24,ctrl = 5,k = FALSE,assay = NULL,name = "NFkB_score_new",seed = 1,search = FALSE,slot = "data")
MS177_organoid_3_conditions$group <- factor(MS177_organoid_3_conditions$long_condition, levels = c("control_24hr_epithelial_like",  "control_24hr_stem-cell_like", "JQ1_24hr_epithelial_like", "JQ1_24hr_stem-cell_like", "MS177_24hr_epithelial_like",  "MS177_24hr_stem-cell_like"))
Idents(MS177_organoid_3_conditions) <- MS177_organoid_3_conditions$dataset
VlnPlot(MS177_organoid_3_conditions, features = c("NFkB_score_new1"), split.by = "group", pt.size = 0)
ggsave("./TNBC_paper_data_and_figure/Organoid_3_conditions_split_by_2_clusters_NFKB_score_doublets_removed_01082025.svg", width = 10, height = 5, dpi = 300)

Idents(MS177_organoid_3_conditions) <- "long_condition"
MS177_organoid_3_conditions_subset_1 <- subset(MS177_organoid_3_conditions, idents = c("control_24hr_epithelial_like"))
MS177_organoid_3_conditions_subset_2 <- subset(MS177_organoid_3_conditions, idents = c("control_24hr_stem-cell_like"))
MS177_organoid_3_conditions_subset_3 <- subset(MS177_organoid_3_conditions, idents = c("JQ1_24hr_epithelial_like"))
MS177_organoid_3_conditions_subset_4 <- subset(MS177_organoid_3_conditions, idents = c("JQ1_24hr_stem-cell_like"))
MS177_organoid_3_conditions_subset_5 <- subset(MS177_organoid_3_conditions, idents = c("MS177_24hr_epithelial_like"))
MS177_organoid_3_conditions_subset_6 <- subset(MS177_organoid_3_conditions, idents = c("MS177_24hr_stem-cell_like"))


median(MS177_organoid_3_conditions_subset_1$NFkB_score_new1)
median(MS177_organoid_3_conditions_subset_2$NFkB_score_new1)
median(MS177_organoid_3_conditions_subset_3$NFkB_score_new1)
median(MS177_organoid_3_conditions_subset_4$NFkB_score_new1)
median(MS177_organoid_3_conditions_subset_5$NFkB_score_new1)
median(MS177_organoid_3_conditions_subset_6$NFkB_score_new1)

wilcox.test(MS177_organoid_3_conditions_subset_1$NFkB_score_new1, MS177_organoid_3_conditions_subset_3$NFkB_score_new1, paired = F) #p-value < 2.2e-16
wilcox.test(MS177_organoid_3_conditions_subset_1$NFkB_score_new1, MS177_organoid_3_conditions_subset_5$NFkB_score_new1, paired = F)#p-value < 2.2e-16
wilcox.test(MS177_organoid_3_conditions_subset_2$NFkB_score_new1, MS177_organoid_3_conditions_subset_4$NFkB_score_new1, paired = F)# p-value < 2.2e-16
wilcox.test(MS177_organoid_3_conditions_subset_2$NFkB_score_new1, MS177_organoid_3_conditions_subset_6$NFkB_score_new1, paired = F)#p-value < 2.2e-16
